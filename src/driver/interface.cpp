#include <unordered_map>
#include <memory>
#include <cstring>
#include <fstream>
#include <mutex>
#include <thread>
#include <anydsl_runtime.hpp>

#if defined(__x86_64__) || defined(__amd64__) || defined(_M_X64)
#include <x86intrin.h>
#endif

#include "interface.h"
#include "bvh.h"
#include "obj.h"
#include "image.h"
#include "buffer.h"

void send_rays(float *, size_t, size_t, bool);
int  recv_rays(float **, size_t, bool, int, bool);
void master_save_ray_batches(float*, size_t, size_t, size_t);
int32_t thread_buffer_size(); 

template <typename Node, typename Tri>
struct Bvh {
    anydsl::Array<Node> nodes;
    anydsl::Array<Tri>  tris;
};

using Bvh2Tri1 = Bvh<Node2, Tri1>;
using Bvh4Tri4 = Bvh<Node4, Tri4>;
using Bvh8Tri4 = Bvh<Node8, Tri4>;

#ifdef ENABLE_EMBREE_DEVICE
#include <embree3/rtcore.h>

template <int N>
struct EmbreeIntersect {};

template <>
struct EmbreeIntersect<8> {
    static void intersect(const int* valid, RTCScene scene, RTCRayHitNt<8>& ray_hit, bool coherent) {
        RTCIntersectContext context;
        context.flags = coherent ? RTC_INTERSECT_CONTEXT_FLAG_COHERENT : RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT;
        context.filter = NULL;
        context.instID[0] = RTC_INVALID_GEOMETRY_ID;
        rtcIntersect8(valid, scene, &context, reinterpret_cast<RTCRayHit8*>(&ray_hit));
    }

    static void occluded(const int* valid, RTCScene scene, RTCRayNt<8>& ray) {
        RTCIntersectContext context;
        context.flags = RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT;
        context.filter = NULL;
        context.instID[0] = RTC_INVALID_GEOMETRY_ID;
        rtcOccluded8(valid, scene, &context, reinterpret_cast<RTCRay8*>(&ray));
    }
};

template <>
struct EmbreeIntersect<4> {
    static void intersect(const int* valid, RTCScene scene, RTCRayHitNt<4>& ray_hit, bool coherent) {
        RTCIntersectContext context;
        context.flags = coherent ? RTC_INTERSECT_CONTEXT_FLAG_COHERENT : RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT;
        context.filter = NULL;
        context.instID[0] = RTC_INVALID_GEOMETRY_ID;
        rtcIntersect4(valid, scene, &context, reinterpret_cast<RTCRayHit4*>(&ray_hit));
    }

    static void occluded(const int* valid, RTCScene scene, RTCRayNt<4>& ray) {
        RTCIntersectContext context;
        context.flags = RTC_INTERSECT_CONTEXT_FLAG_INCOHERENT;
        context.filter = NULL;
        context.instID[0] = RTC_INVALID_GEOMETRY_ID;
        rtcOccluded4(valid, scene, &context, reinterpret_cast<RTCRay4*>(&ray));
    }
};

struct EmbreeDevice {
    std::vector<uint32_t> indices;
    std::vector<float3>   vertices;
    RTCDevice device;
    RTCScene scene;

    static void error_handler(void*, const RTCError code, const char* str) {
        if (code == RTC_ERROR_NONE)
            return;
        std::cerr << "Embree error ";
        switch (code) {
            case RTC_ERROR_UNKNOWN:             std::cerr << "RTC_ERROR_UNKNOWN"; break;
            case RTC_ERROR_INVALID_ARGUMENT:    std::cerr << "RTC_ERROR_INVALID_ARGUMENT"; break;
            case RTC_ERROR_INVALID_OPERATION:   std::cerr << "RTC_ERROR_INVALID_OPERATION"; break;
            case RTC_ERROR_OUT_OF_MEMORY:       std::cerr << "RTC_ERROR_OUT_OF_MEMORY"; break;
            case RTC_ERROR_UNSUPPORTED_CPU:     std::cerr << "RTC_ERROR_UNSUPPORTED_CPU"; break;
            case RTC_ERROR_CANCELLED:           std::cerr << "RTC_ERROR_CANCELLED"; break;
            default: break;
        }
        if (str) std::cerr << ": " << str;
        std::cerr << std::endl;
        abort();
    }

    EmbreeDevice() {
        read_buffer("data/vertices.bin", vertices);
        read_buffer("data/indices.bin", indices);

        if (vertices.empty() || indices.empty())
            error("Cannot build scene due to missing data files");

        device = rtcNewDevice(nullptr);
        if (!device)
            error("Cannot initialize Embree");
        rtcSetDeviceErrorFunction(device, error_handler, nullptr);
        auto mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);
        scene = rtcNewScene(device);

        auto vertex_ptr = (float4*)rtcSetNewGeometryBuffer(mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, sizeof(float) * 4, vertices.size());
        for (auto& v : vertices)
            *(vertex_ptr++) = float4(v, 1.0f);
        auto index_ptr = (uint32_t*)rtcSetNewGeometryBuffer(mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, sizeof(uint32_t) * 4, indices.size() / 4);
        for (auto& i : indices)
            *(index_ptr++) = i;

        rtcCommitGeometry(mesh);
        rtcAttachGeometry(scene, mesh);
        rtcCommitScene(scene);

        info("Embree device initialized successfully");
    }

    ~EmbreeDevice() {
        rtcReleaseScene(scene);
        rtcReleaseDevice(device);
    }

    template <int N>
    void load_ray(RTCRayHitNt<N>& ray_hit, const RayStream& rays, size_t i, size_t size) {
        size_t n = std::min(size, i + N);
        for (size_t j = i, k = 0; j < n; ++j, ++k) {
            ray_hit.ray.org_x[k] = rays.org_x[j];
            ray_hit.ray.org_y[k] = rays.org_y[j];
            ray_hit.ray.org_z[k] = rays.org_z[j];

            ray_hit.ray.dir_x[k] = rays.dir_x[j];
            ray_hit.ray.dir_y[k] = rays.dir_y[j];
            ray_hit.ray.dir_z[k] = rays.dir_z[j];

            ray_hit.ray.tnear[k] = rays.tmin[j];
            ray_hit.ray.tfar[k]  = rays.tmax[j];
            ray_hit.ray.mask[k]  = 0xFFFFFFFF;
            ray_hit.ray.id[k]    = i;
            ray_hit.ray.flags[k] = 0;

            ray_hit.hit.primID[k] = -1;
            ray_hit.hit.geomID[k] = -1;
        }
    }

    template <int N>
    void load_ray(RTCRayNt<N>& ray, const RayStream& rays, size_t i, size_t size) {
        size_t n = std::min(size, i + N);
        for (size_t j = i, k = 0; j < n; ++j, ++k) {
            ray.org_x[k] = rays.org_x[j];
            ray.org_y[k] = rays.org_y[j];
            ray.org_z[k] = rays.org_z[j];

            ray.dir_x[k] = rays.dir_x[j];
            ray.dir_y[k] = rays.dir_y[j];
            ray.dir_z[k] = rays.dir_z[j];

            ray.tnear[k] = rays.tmin[j];
            ray.tfar[k]  = rays.tmax[j];
            ray.mask[k]  = 0xFFFFFFFF;
            ray.id[k]    = i;
            ray.flags[k] = 0;
        }
    }

    template <int N>
    void save_hit(const RTCRayHitNt<N>& ray_hit, PrimaryStream& primary, size_t i, int32_t invalid_id) {
        size_t n = std::min(primary.size, int32_t(i + N));
        for (size_t j = i, k = 0; j < n; ++j, ++k) {
            auto prim_id = ray_hit.hit.primID[k];
            primary.geom_id[j] = prim_id == RTC_INVALID_GEOMETRY_ID ? invalid_id : indices[prim_id * 4 + 3];
            primary.prim_id[j] = prim_id;
            primary.t[j]       = ray_hit.ray.tfar[k];
            primary.u[j]       = ray_hit.hit.u[k];
            primary.v[j]       = ray_hit.hit.v[k];
        }
    }

    template <int N>
    void save_hit(const RTCRayNt<N>& ray, SecondaryStream& secondary, size_t i) {
        size_t n = std::min(secondary.size, int32_t(i + N));
        for (size_t j = i, k = 0; j < n; ++j, ++k) {
            secondary.prim_id[j] = ray.tfar[k] < 0 ? 0 : -1;
        }
    }

    template <int N>
    void intersect(PrimaryStream& primary, int32_t invalid_id, bool coherent) {
        alignas(32) int valid[N];
        alignas(32) RTCRayHitNt<N> ray_hit;
        memset(valid, 0xFF, sizeof(int) * N);
        for (size_t i = 0; i < primary.size; i += N) {
            load_ray(ray_hit, primary.rays, i, primary.size);
            EmbreeIntersect<N>::intersect(valid, scene, ray_hit, coherent);
            save_hit(ray_hit, primary, i, invalid_id);
        }
    }

    template <int N>
    void intersect(SecondaryStream& secondary) {
        alignas(32) int valid[N];
        alignas(32) RTCRayNt<N> ray;
        memset(valid, 0xFF, sizeof(int) * N);
        for (size_t i = 0; i < secondary.size; i += N) {
            load_ray(ray, secondary.rays, i, secondary.size);
            EmbreeIntersect<N>::occluded(valid, scene, ray);
            save_hit(ray, secondary, i);
        }
    }
};
#endif

struct Interface {
    using DeviceImage = std::tuple<anydsl::Array<uint8_t>, int32_t, int32_t>;

    struct DeviceData {
        std::unordered_map<std::string, Bvh2Tri1> bvh2_tri1;
        std::unordered_map<std::string, Bvh4Tri4> bvh4_tri4;
        std::unordered_map<std::string, Bvh8Tri4> bvh8_tri4;
        std::unordered_map<std::string, anydsl::Array<uint8_t>> buffers;
        std::unordered_map<std::string, DeviceImage> images;
        anydsl::Array<int32_t> tmp_buffer;
        
        anydsl::Array<float> first_primary;
        anydsl::Array<float> second_primary;
        anydsl::Array<float> outgoing_primary;
        
        anydsl::Array<float> first_secondary;
        anydsl::Array<float> second_secondary;
        anydsl::Array<float> outgoing_secondary;
        anydsl::Array<float> film_pixels;
    };
    std::unordered_map<int32_t, DeviceData> devices;
    
    static thread_local anydsl::Array<float> cpu_primary;
    static thread_local anydsl::Array<float> cpu_secondary;
    static thread_local anydsl::Array<float> cpu_primary_outgoing;
    static thread_local anydsl::Array<float> cpu_secondary_outgoing;
    static thread_local anydsl::Array<float> cpu_film;

#ifdef ENABLE_EMBREE_DEVICE
    EmbreeDevice embree_device;
#endif

    anydsl::Array<float> host_pixels;
    anydsl::Array<float> tmp_pixels;
    anydsl::Array<float> host_primary;
    anydsl::Array<float> host_secondary;
    anydsl::Array<float> host_outgoing_primary;
    anydsl::Array<float> host_outgoing_secondary;
    anydsl::Array<int32_t> host_primary_num;
    anydsl::Array<int32_t> host_prerender;

    size_t film_width;
    size_t film_height;
    
    Interface(size_t width, size_t height)
        : film_width(width)
        , film_height(height)
    {
        
        host_pixels = anydsl::Array<float>(width * height * 3);
        tmp_pixels  = anydsl::Array<float>(width * height * 3);
        
        host_primary_num      = anydsl::Array<int32_t>(128);
        host_prerender        = anydsl::Array<int32_t>(width * height * 8/*maxlength*/);
    }

    template <typename T>
    anydsl::Array<T>& resize_array(int32_t dev, anydsl::Array<T>& array, size_t size, size_t multiplier) {
        auto capacity = (size & ~((1 << 5) - 1)) + 32; // round to 32
        if (array.size() < capacity) {
            auto n = capacity * multiplier;
            array = std::move(anydsl::Array<T>(dev, reinterpret_cast<T*>(anydsl_alloc(dev, sizeof(T) * n)), n));
        }
        return array;
    }

    //cpu 
    anydsl::Array<float>& cpu_primary_stream(size_t size) {
        return resize_array(0, cpu_primary, size, 21);
    }

    anydsl::Array<float>& cpu_primary_outgoing_stream(size_t size) {
        return resize_array(0, cpu_primary_outgoing, size, 21);
    }
    
    anydsl::Array<float>& cpu_secondary_stream(size_t size) {
        return resize_array(0, cpu_secondary, size, 14);
    }
    
    anydsl::Array<float>& cpu_secondary_outgoing_stream(size_t size) {
        return resize_array(0, cpu_secondary_outgoing, size, 14);
    }
    //gpu
    void initial_gpu_host_data(size_t size) {
        auto capacity = (size & ~((1 << 5) - 1)) + 32;
        resize_array(0, host_primary,  capacity, 21); 
        resize_array(0, host_outgoing_primary, capacity, 21); 
        resize_array(0, host_secondary,  capacity, 14); 
        resize_array(0, host_outgoing_secondary, capacity, 14); 
    } 
    
    anydsl::Array<float>& gpu_first_primary_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].first_primary, size, 21);
    }

    anydsl::Array<float>& gpu_second_primary_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].second_primary, size, 21);
    }
    
    anydsl::Array<float>& gpu_primary_buffer_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].outgoing_primary, size, 21);
    }
    anydsl::Array<float>& gpu_outgoing_primary_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].outgoing_primary, size, 21);
    }

    anydsl::Array<float>& gpu_first_secondary_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].first_secondary, size, 14);
    }

    anydsl::Array<float>& gpu_second_secondary_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].second_secondary, size, 14);
    }

    anydsl::Array<float>& gpu_secondary_buffer_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].outgoing_secondary, size, 14);
    }
    
    anydsl::Array<float>& gpu_outgoing_secondary_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].outgoing_secondary, size, 14);
    }
    
    anydsl::Array<int32_t>& gpu_tmp_buffer(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].tmp_buffer, size, 1);
    }

    //bvh
    const Bvh2Tri1& load_bvh2_tri1(int32_t dev, const std::string& filename) {
        auto& bvh2_tri1 = devices[dev].bvh2_tri1;
        auto it = bvh2_tri1.find(filename);
        if (it != bvh2_tri1.end())
            return it->second;
        return bvh2_tri1[filename] = std::move(load_bvh<Node2, Tri1>(dev, filename));
    }

    const Bvh4Tri4& load_bvh4_tri4(int32_t dev, const std::string& filename) {
        auto& bvh4_tri4 = devices[dev].bvh4_tri4;
        auto it = bvh4_tri4.find(filename);
        if (it != bvh4_tri4.end())
            return it->second;
        return bvh4_tri4[filename] = std::move(load_bvh<Node4, Tri4>(dev, filename));
    }

    const Bvh8Tri4& load_bvh8_tri4(int32_t dev, const std::string& filename) {
        printf("before bvh device\n");
        auto& bvh8_tri4 = devices[dev].bvh8_tri4;
        printf("after bvh device\n");
        auto it = bvh8_tri4.find(filename);
        printf("after bvh find file\n");
        if (it != bvh8_tri4.end())
            return it->second;
        printf("before bvh std::move\n");
        return bvh8_tri4[filename] = std::move(load_bvh<Node8, Tri4>(dev, filename));
    }

    template <typename T>
    anydsl::Array<T> copy_to_device(int32_t dev, const T* data, size_t n) {
        anydsl::Array<T> array(dev, reinterpret_cast<T*>(anydsl_alloc(dev, n * sizeof(T))), n);
        anydsl_copy(0, data, 0, dev, array.data(), 0, sizeof(T) * n);
        return array;
    }

    template <typename T>
    anydsl::Array<T> copy_to_device(int32_t dev, const std::vector<T>& host) {
        return copy_to_device(dev, host.data(), host.size());
    }

    DeviceImage copy_to_device(int32_t dev, const ImageRgba32& img) {
        return DeviceImage(copy_to_device(dev, img.pixels.get(), img.width * img.height * 4), img.width, img.height);
    }

    template <typename Node, typename Tri>
    Bvh<Node, Tri> load_bvh(int32_t dev, const std::string& filename) {
        std::ifstream is(filename, std::ios::binary);
        if (!is)
            error("Cannot open BVH '", filename, "'");
        do {
            size_t node_size = 0, tri_size = 0;
            is.read((char*)&node_size, sizeof(uint32_t));
            is.read((char*)&tri_size,  sizeof(uint32_t));
            if (node_size == sizeof(Node) &&
                tri_size  == sizeof(Tri)) {
                info("Loaded BVH file '", filename, "'");
                std::vector<Node> nodes;
                std::vector<Tri>  tris;
                printf("before read nodes\n");
                read_buffer(is, nodes);
                printf("before read tris\n");
                read_buffer(is, tris);
                printf("before move copy to device\n");
                return Bvh<Node, Tri> { std::move(copy_to_device(dev, nodes)), std::move(copy_to_device(dev, tris)) };
            }
            skip_buffer(is);
            skip_buffer(is);
        } while (!is.eof() && is);
        error("Invalid BVH file");
    }

    const anydsl::Array<uint8_t>& load_buffer(int32_t dev, const std::string& filename) {
        info("Prepare Load buffer '", filename, "'");
        auto& buffers = devices[dev].buffers;
        auto it = buffers.find(filename);
        if (it != buffers.end())
            return it->second;
        std::ifstream is(filename, std::ios::binary);
        if (!is)
            error("Cannot open buffer '", filename, "'");
        std::vector<uint8_t> vector;
        read_buffer(is, vector);
        info("Loaded buffer '", filename, "'");
        buffers[filename] = std::move(copy_to_device(dev, vector));
        info("after buffer copy'", filename, "'");
        return buffers[filename];
    }

    anydsl::Array<float>& load_gpu_primary_stream(int dev, int32_t size, const std::string& filename) {
        auto& primary = devices[dev].first_primary; 
        //    gpu_first_primary_stream(dev, size); 
        std::ifstream is(filename, std::ios::binary);
        if (!is)
           error("Cannot open primary buffer");
        read_buffer(is, host_primary);
        info("loaded buffer primary buffer '", filename, "'"); 
        anydsl::copy(host_primary, primary);
        return primary;
    }
    
    const DeviceImage& load_png(int32_t dev, const std::string& filename) {
        auto& images = devices[dev].images;
        auto it = images.find(filename);
        if (it != images.end())
            return it->second;
        ImageRgba32 img;
        if (!::load_png(filename, img))
            error("Cannot load PNG file '", filename, "'");
        info("Loaded PNG file '", filename, "'");
        return images[filename] = std::move(copy_to_device(dev, img));
    }

    const DeviceImage& load_jpg(int32_t dev, const std::string& filename) {
        printf("before load jpg\n");
        auto& images = devices[dev].images;
        printf("image \n");
        auto it = images.find(filename);
        printf("it \n");
        if (it != images.end())
            return it->second;
        ImageRgba32 img;
        printf("img \n");
        if (!::load_jpg(filename, img))
            error("Cannot load JPG file '", filename, "'");
        info("Loaded JPG file '", filename, "'");
        return images[filename] = std::move(copy_to_device(dev, img));
    }

    void present(int32_t dev) {
        anydsl::copy(devices[dev].film_pixels, host_pixels);
     //   int size = film_width * film_height * 3;
     //   for (int i = 0; i < size; i++){
     //       host_pixels[i] += tmp_pixels[i];
     //   }
    }
    
    int save_first_primary(int32_t dev) {
        anydsl::copy(devices[dev].first_primary, host_primary);
        return devices[dev].first_primary.size();
    }
    int save_second_primary(int32_t dev) {
        anydsl::copy(devices[dev].second_primary, host_primary);
        return devices[dev].second_primary.size();
    }
    void load_first_primary(int32_t dev) {
        anydsl::copy(host_primary, devices[dev].first_primary);
    }
    void load_second_primary(int32_t dev) {
        anydsl::copy(host_primary, devices[dev].second_primary);
    }

    int save_first_secondary(int32_t dev) {
        anydsl::copy(devices[dev].first_secondary, host_secondary);
        return devices[dev].first_secondary.size();
    }
    int save_second_secondary(int32_t dev) {
        anydsl::copy(devices[dev].second_secondary, host_secondary);
        return devices[dev].second_secondary.size();
    }
    void load_first_secondary(int32_t dev) {
        anydsl::copy(host_secondary, devices[dev].first_secondary);
    }
    void load_second_secondary(int32_t dev) {
        anydsl::copy(host_secondary, devices[dev].second_secondary);
    }
    void save_first_primary_num(int32_t dev) {
        anydsl::copy(devices[dev].tmp_buffer, host_primary_num);
    }
    int save_outgoing_primary(int32_t dev) {
        anydsl::copy(devices[dev].outgoing_primary, host_outgoing_primary);
        return devices[dev].outgoing_primary.size();
    }
    int save_buffer_secondary(int32_t dev) {
        anydsl::copy(devices[dev].outgoing_secondary, host_outgoing_secondary);
        return devices[dev].outgoing_secondary.size();
    }
    void load_buffer_primary(int32_t dev) {
        anydsl::copy(host_outgoing_primary, devices[dev].outgoing_primary);
    }
    void clear() {
        std::fill(host_pixels.begin(), host_pixels.end(), 0.0f);
        //clear cpu thread film
        for (auto& pair : devices) {
            auto& device_pixels = devices[pair.first].film_pixels;
            if (device_pixels.size())
                anydsl::copy(host_pixels, device_pixels);
        }
    }
};
thread_local anydsl::Array<float> Interface::cpu_primary;
thread_local anydsl::Array<float> Interface::cpu_secondary;
thread_local anydsl::Array<float> Interface::cpu_primary_outgoing;
thread_local anydsl::Array<float> Interface::cpu_secondary_outgoing;
thread_local anydsl::Array<float> Interface::cpu_film;

static std::unique_ptr<Interface> interface;

void setup_interface(size_t width, size_t height) {
    interface.reset(new Interface(width, height));
}

void cleanup_interface() {
    interface.reset();
}

float* get_pixels() {
    return interface->host_pixels.data();
}

int32_t* get_prerender_result() {
    return interface->host_prerender.data();
}

float* get_first_primary() {
    return interface->host_primary.data();    
}

float* get_buffer_primary() {
    return interface->host_outgoing_primary.data();
}

void clear_pixels() {
    return interface->clear();
}

inline void get_ray_stream(RayStream& rays, float* ptr, size_t capacity) {
    rays.id    = (int*)ptr;
    rays.org_x = ptr + 1 * capacity;
    rays.org_y = ptr + 2 * capacity;
    rays.org_z = ptr + 3 * capacity;
    rays.dir_x = ptr + 4 * capacity;
    rays.dir_y = ptr + 5 * capacity;
    rays.dir_z = ptr + 6 * capacity;
    rays.tmin  = ptr + 7 * capacity;
    rays.tmax  = ptr + 8 * capacity;
}

inline void get_primary_stream(PrimaryStream& primary, float* ptr, size_t capacity) {
    get_ray_stream(primary.rays, ptr, capacity);
    primary.grid_id   = (int*)ptr + 9 * capacity;
    primary.geom_id   = (int*)ptr + 10 * capacity;
    primary.prim_id   = (int*)ptr + 11 * capacity;
    primary.t         = ptr + 12 * capacity;
    primary.u         = ptr + 13 * capacity;
    primary.v         = ptr + 14 * capacity;
    primary.rnd       = (unsigned int*)ptr + 15 * capacity;
    primary.mis       = ptr + 16 * capacity;
    primary.contrib_r = ptr + 17 * capacity;
    primary.contrib_g = ptr + 18 * capacity;
    primary.contrib_b = ptr + 19 * capacity;
    primary.depth     = (int*)ptr + 20 * capacity;
    primary.size = 0;
}

inline void get_secondary_stream(SecondaryStream& secondary, float* ptr, size_t capacity) {
    get_ray_stream(secondary.rays, ptr, capacity);
    secondary.grid_id = (int*)ptr + 9 * capacity;
    secondary.prim_id = (int*)ptr + 10 * capacity;
    secondary.color_r = ptr + 11 * capacity;
    secondary.color_g = ptr + 12 * capacity;
    secondary.color_b = ptr + 13 * capacity;
    secondary.size = 0;
}

inline void get_out_ray_stream(OutRayStream& out, float* ptr, size_t capacity, size_t width) {
    out.fptr  = ptr;
    out.iptr  = (int*)ptr;
    out.uiptr = (unsigned int*) ptr;
    out.capacity = capacity;
    out.width = width;
    out.size = 0;
}

inline void copy_primary_ray(PrimaryStream a, PrimaryStream b, int src_id, int dst_id, bool keep_hit) {
   b.rays.id[dst_id]    = a.rays.id[src_id];
   b.rays.org_x[dst_id] = a.rays.org_x[src_id];
   b.rays.org_y[dst_id] = a.rays.org_y[src_id];
   b.rays.org_z[dst_id] = a.rays.org_z[src_id];
   b.rays.dir_x[dst_id] = a.rays.dir_x[src_id];
   b.rays.dir_y[dst_id] = a.rays.dir_y[src_id];
   b.rays.dir_z[dst_id] = a.rays.dir_z[src_id];
   b.rays.tmin[dst_id]  = a.rays.tmin[src_id];
   b.rays.tmax[dst_id]  = a.rays.tmax[src_id];
   b.grid_id[dst_id]    = a.grid_id[src_id];
   if (keep_hit) {
       b.geom_id[dst_id] = a.geom_id[src_id];
       b.prim_id[dst_id] = a.prim_id[src_id];
       b.t[dst_id]       = a.t[src_id];
       b.u[dst_id]       = a.u[src_id];
       b.v[dst_id]       = a.v[src_id];
   }
   b.rnd[dst_id]        = a.rnd[src_id];
   b.mis[dst_id]        = a.mis[src_id];
   b.contrib_r[dst_id]  = a.contrib_r[src_id];
   b.contrib_g[dst_id]  = a.contrib_g[src_id];
   b.contrib_b[dst_id]  = a.contrib_b[src_id];
   b.depth[dst_id]      = a.depth[src_id];
}

extern "C" {

//void rodent_get_film_data(int32_t dev, float** pixels, int32_t* width, int32_t* height) {
//    if (dev != 0) {
//        auto& device = interface->devices[dev];
//        if (!device.film_pixels.size()) {
//            auto film_size = interface->film_width * interface->film_height * 3;
//            auto film_data = reinterpret_cast<float*>(anydsl_alloc(dev, sizeof(float) * film_size));
//            device.film_pixels = std::move(anydsl::Array<float>(dev, film_data, film_size));
//            anydsl::copy(interface->host_pixels, device.film_pixels);
//        }
//        *pixels = device.film_pixels.data();
//    } else {
//        *pixels = interface->host_pixels.data();
//    }
//    *width  = interface->film_width;
//    *height = interface->film_height;
//}

void rodent_get_film_data(int32_t dev, float** pixels, int32_t* width, int32_t* height, bool host) {
    if (dev != 0) {
        auto& device = interface->devices[dev];
        if (!device.film_pixels.size()) {
            auto film_size = interface->film_width * interface->film_height * 3;
            auto film_data = reinterpret_cast<float*>(anydsl_alloc(dev, sizeof(float) * film_size));
            device.film_pixels = std::move(anydsl::Array<float>(dev, film_data, film_size));
            anydsl::copy(interface->host_pixels, device.film_pixels);
        }
        *pixels = device.film_pixels.data();
    } else {
        *pixels = interface->host_pixels.data();
    }
    *width  = interface->film_width;
    *height = interface->film_height;

//    if(host) {
//       *pixels = interface->host_pixels.data();
//    } else {
//        auto film_size = interface->film_width * interface->film_height * 3;
//        auto film_data = reinterpret_cast<float*>(anydsl_alloc(dev, sizeof(float) * film_size));
//        // thread film 
//        printf("get film data dev %d\n", dev);
//        if(dev == 0) { 
//            interface->cpu_film = std::move(anydsl::Array<float>(dev, film_data, film_size));
//            anydsl::copy(interface->host_pixels, interface->cpu_film);
//            *pixels = interface->cpu_film.data();
//        } else {
//            auto& device = interface->devices[dev];
//            if (!device.film_pixels.size()) {
//                device.film_pixels = std::move(anydsl::Array<float>(dev, film_data, film_size));
//                anydsl::copy(interface->host_pixels, device.film_pixels);
//            }
//            *pixels = device.film_pixels.data();
//        } 
//    }
//    *width  = interface->film_width;
//    *height = interface->film_height;
}

void rodent_get_prerender_data(int** prerender) {
    std::fill(interface->host_prerender.begin(), interface->host_prerender.end(), -1);
    *prerender = interface->host_prerender.data();
}

void rodent_load_png(int32_t dev, const char* file, uint8_t** pixels, int32_t* width, int32_t* height) {
    auto& img = interface->load_png(dev, file);
    *pixels = const_cast<uint8_t*>(std::get<0>(img).data());
    *width  = std::get<1>(img);
    *height = std::get<2>(img);
}

void rodent_load_jpg(int32_t dev, const char* file, uint8_t** pixels, int32_t* width, int32_t* height) {
    auto& img = interface->load_jpg(dev, file);
    *pixels = const_cast<uint8_t*>(std::get<0>(img).data());
    *width  = std::get<1>(img);
    *height = std::get<2>(img);
}

uint8_t* rodent_load_buffer(int32_t dev, const char* file) {
    auto& array = interface->load_buffer(dev, file);
    return const_cast<uint8_t*>(array.data());
}

void rodent_load_bvh2_tri1(int32_t dev, const char* file, Node2** nodes, Tri1** tris) {
    auto& bvh = interface->load_bvh2_tri1(dev, file);
    *nodes = const_cast<Node2*>(bvh.nodes.data());
    *tris  = const_cast<Tri1*>(bvh.tris.data());
}

void rodent_load_bvh4_tri4(int32_t dev, const char* file, Node4** nodes, Tri4** tris) {
    auto& bvh = interface->load_bvh4_tri4(dev, file);
    *nodes = const_cast<Node4*>(bvh.nodes.data());
    *tris  = const_cast<Tri4*>(bvh.tris.data());
}

void rodent_load_bvh8_tri4(int32_t dev, const char* file, Node8** nodes, Tri4** tris) {
    auto& bvh = interface->load_bvh8_tri4(dev, file);
    printf("load bvh 8 4\n");
    *nodes = const_cast<Node8*>(bvh.nodes.data());
    *tris  = const_cast<Tri4*>(bvh.tris.data());
}

//cpu
void rodent_cpu_get_primary_stream(PrimaryStream* primary, int32_t size) {
    auto& array = interface->cpu_primary_stream(size);
    get_primary_stream(*primary, array.data(), array.size() / 21);
}

void rodent_cpu_get_out_ray_stream(OutRayStream * stream, int32_t capacity, bool primary) {
    if(primary) {
        auto& array = interface->cpu_primary_outgoing_stream(capacity);
        get_out_ray_stream(*stream, array.data(), capacity, 21);
    } else {
        auto& array = interface->cpu_secondary_outgoing_stream(capacity);
        get_out_ray_stream(*stream, array.data(), capacity, 14);
    }
}

void rodent_cpu_get_secondary_stream(SecondaryStream* secondary, int32_t size) {
    auto& array = interface->cpu_secondary_stream(size);
    get_secondary_stream(*secondary, array.data(), array.size() / 14);
}

void rodent_cpu_get_secondary_outgoing_stream(SecondaryStream* buffer, int32_t size) {
    auto& array = interface->cpu_secondary_outgoing_stream(size);
    get_secondary_stream(*buffer, array.data(), array.size() / 14);
}

int32_t rodent_cpu_get_thread_num() {
    printf("hardware concurrency %d\n", std::thread::hardware_concurrency());
    return 1;//std::thread::hardware_concurrency(); 
}
//gpu
//
void rodent_initial_gpu_host_data(int32_t size) {
    interface->initial_gpu_host_data(size); 
}

void rodent_gpu_get_tmp_buffer(int32_t dev, int32_t** buf, int32_t size) {
    *buf = interface->gpu_tmp_buffer(dev, size).data();
}

void rodent_gpu_get_first_primary_stream(int32_t dev, PrimaryStream* primary, int32_t size) {
    auto& array = interface->gpu_first_primary_stream(dev, size);
    get_primary_stream(*primary, array.data(), array.size() / 21);
}

void rodent_gpu_get_second_primary_stream(int32_t dev, PrimaryStream* primary, int32_t size) {
    auto& array = interface->gpu_second_primary_stream(dev, size);
    get_primary_stream(*primary, array.data(), array.size() / 21);
}

void rodent_gpu_get_primary_buffer_stream(int32_t dev, PrimaryStream* primary, int32_t size) {
    auto& array = interface->gpu_primary_buffer_stream(dev, size);
    get_primary_stream(*primary, array.data(), array.size() / 21);
}


void rodent_gpu_get_first_secondary_stream(int32_t dev, SecondaryStream* secondary, int32_t size) {
    auto& array = interface->gpu_first_secondary_stream(dev, size);
    get_secondary_stream(*secondary, array.data(), array.size() / 14);
}

void rodent_gpu_get_second_secondary_stream(int32_t dev, SecondaryStream* secondary, int32_t size) {
    auto& array = interface->gpu_second_secondary_stream(dev, size);
    get_secondary_stream(*secondary, array.data(), array.size() / 14);
}

void rodent_gpu_get_secondary_buffer_stream(int32_t dev, SecondaryStream* secondary, int32_t size) {
    auto& array = interface->gpu_secondary_buffer_stream(dev, size);
    get_secondary_stream(*secondary, array.data(), array.size() / 14);
}

void rodent_gpu_get_out_ray_stream(int32_t dev, OutRayStream * stream, int32_t capacity, bool primary) {
    if(primary) {
        auto& array = interface->gpu_outgoing_primary_stream(dev, capacity);
        get_out_ray_stream(*stream, array.data(), capacity, 21);
    } else {
        auto& array = interface->gpu_outgoing_secondary_stream(dev, capacity);
        get_out_ray_stream(*stream, array.data(), capacity, 14);
    }
}

//embree
#ifdef ENABLE_EMBREE_DEVICE
void rodent_cpu_intersect_primary_embree(PrimaryStream* primary, int32_t invalid_id, int32_t coherent) {
    interface->embree_device.intersect<8>(*primary, invalid_id, coherent != 0);
}

void rodent_cpu_intersect_secondary_embree(SecondaryStream* secondary) {
    interface->embree_device.intersect<8>(*secondary);
}
#endif
 
//debug
void rodent_secondary_check(int32_t dev, int primary_size, int32_t print_mark, bool is_first_primary) {
    SecondaryStream secondary;
    printf("\n%d secondary size  %d\n", print_mark, primary_size);
    if (dev == -1) {
        auto& array = interface->cpu_secondary;
        get_secondary_stream(secondary, array.data(), array.size() / 14);
    } 
    for(int i = 0; i < 5; i++){
        printf("ray id %d chunk %d tmin %f ray org %f %f %f dir %f %f %f\n", 
                    secondary.rays.id[i],    secondary.grid_id[i],    secondary.rays.tmin[i],
                    secondary.rays.org_x[i], secondary.rays.org_y[i], secondary.rays.org_z[i],
                    secondary.rays.dir_x[i], secondary.rays.dir_y[i], secondary.rays.dir_z[i]);
    }
    SecondaryStream buffer;
    auto& array = interface->cpu_secondary_outgoing;
    get_secondary_stream(buffer, array.data(), array.size() / 14);
    for(int i = 0; i < 3; i++){
        printf("buffer id %d chunk %d tmin %f ray org %f %f %f dir %f %f %f\n", 
                    buffer.rays.id[i], buffer.grid_id[i], buffer.rays.tmin[i],
                    buffer.rays.org_x[i], buffer.rays.org_y[i], buffer.rays.org_z[i],
                    buffer.rays.dir_x[i], buffer.rays.dir_y[i], buffer.rays.dir_z[i]);
    }
}

void rodent_first_primary_check(int32_t dev, int primary_size, int32_t print_mark, bool is_first_primary) {
    PrimaryStream primary;
    printf("%d primary size %d\n", print_mark, primary_size);
    if (dev == -1) {
        auto& array = interface->cpu_primary;
        get_primary_stream(primary, array.data(), array.size() / 21);
    } else {
        int capacity = is_first_primary? interface->save_first_primary(dev) / 21 : interface->save_second_primary(dev) / 21;
        auto&  array = interface->host_primary;
        get_primary_stream(primary, array.data(), capacity);
    }
    for(int i = 0; i < 5; i++) {
        printf("ray id %d chunk %d tmin %f ray org %f %f %f dir %f %f %f geom %d prim %d\n", 
                    primary.rays.id[i+0], primary.grid_id[i+0], primary.rays.tmin[i],
                    primary.rays.org_x[i], primary.rays.org_y[i], primary.rays.org_z[i],
                    primary.rays.dir_x[i], primary.rays.dir_y[i], primary.rays.dir_z[i],
                    primary.geom_id[i], primary.prim_id[i]);
    }
}

int rodent_first_primary_save(int32_t dev, int primary_size, int32_t chunk_num, bool is_first_primary) {
    int capacity = is_first_primary? interface->save_first_primary(dev) / 21 : interface->save_second_primary(dev) / 21;
    auto&  array = interface->host_primary;
    int size_new = 0; //rays_transfer(array.data(), primary_size, capacity);
    is_first_primary ? interface->load_first_primary(dev) : interface->load_second_primary(dev);
    return size_new;
}

//worker

int32_t rodent_thread_buffer_size() {
    return thread_buffer_size(); 
}

void rodent_worker_primary_send(int32_t dev, int buffer_size) {
    if(dev == -1) {
        auto& array = interface->cpu_primary_outgoing;
        printf("interface before send primary\n");
        send_rays(array.data(), buffer_size, array.size() / 21, true);
        printf("interface after send primary\n");
    } else {
        int capacity = interface->save_outgoing_primary(dev) / 21;
        auto& array  = interface->host_outgoing_primary;
        printf("interface before send secondary\n");
        send_rays(array.data(), buffer_size, capacity, true);
        printf("interface after send secondary\n");
    }
}

void rodent_worker_secondary_send(int32_t dev, int buffer_size) {
    if(dev == -1) {
        auto& array = interface->cpu_secondary_outgoing;
//            printf("array.size %d buffer size %d\n ", array.size() / 14, buffer_size);
        send_rays(array.data(), buffer_size, array.size() / 14, false);
    } else {
        int capacity = interface->save_buffer_secondary(dev) / 14;
        auto& array  = interface->host_outgoing_secondary;
        send_rays(array.data(), buffer_size, capacity, false);
    }
}

int32_t rodent_worker_primary_recv(int32_t dev, int32_t rays_size, bool isFirst, int32_t thread_id, bool thread_wait) {
    int size_new = 0;
    if(dev == -1){
        float* array = interface->cpu_primary.data();
        size_new = recv_rays(&array, rays_size, true, thread_id, thread_wait); 
    } else {
        float* array = interface->host_primary.data();
        size_new = recv_rays(&array, rays_size, true, thread_id, thread_wait);
        if(size_new > 0)
            isFirst ? interface->load_first_primary(dev) : interface->load_second_primary(dev);
    }
    return size_new;
}

int32_t rodent_worker_secondary_recv(int32_t dev, int32_t rays_size, bool isFirst, int32_t thread_id, bool thread_wait) {
    int size_new = 0;
    if(dev == -1){
        float* array = interface->cpu_secondary.data();
        size_new = recv_rays(&array, rays_size, false, thread_id, thread_wait); 
    } else {
        float* array = interface->host_secondary.data();
        size_new = recv_rays(&array, rays_size, false, thread_id, thread_wait);
        if(size_new > 0) {
            isFirst ? interface->load_first_secondary(dev) : interface->load_second_secondary(dev);
        }
    }
    return size_new;
}

//master
void rodent_master_save_ray_batches(int32_t dev, int32_t thread_id, int32_t rays_size) {
    if(dev == -1) {
//        printf("interface tid%d size %d\n", thread_id, rays_size);
        auto& array = interface->cpu_primary;
        master_save_ray_batches(array.data(), rays_size, array.size() / 21, thread_id);
    } else {
        printf("master only cpu\n");
//        int capacity = interface->save_buffer_primary(dev) / 21;
//        auto& array  = interface->host_buffer_primary;
//        master_save_rays(array.data(), buffer_size, capacity, true, send_all);
    }
}

void rodent_gpu_first_primary_load(int32_t dev, PrimaryStream* primary, int32_t size) {
    auto& array = interface->load_gpu_primary_stream(dev, size, "data/primary.bin");
    size = array.size() / 21;
    get_primary_stream(*primary, array.data(), size);
}

void rodent_present(int32_t dev) {
    if (dev != 0)
        interface->present(dev);
}

int64_t clock_us() {
#if defined(__x86_64__) || defined(__amd64__) || defined(_M_X64)
#define CPU_FREQ 4e9
    return __rdtsc() * int64_t(1000000) / int64_t(CPU_FREQ);
#else
    using namespace chrono;
    return duration_cast<microseconds>(high_resolution_clock::now().time_since_epoch());
#endif
}

} // extern "C"

