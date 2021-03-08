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
int  recv_rays(float **, bool, int);
void master_save_ray_batches(float*, size_t, size_t, size_t);

int dfw_thread_num();
int dfw_chunk_size();
int dfw_mpi_rank();
int dfw_stream_logic_capacity(); 
int dfw_stream_store_capacity();
int dfw_out_stream_capacity();
void dfw_time_start(std::string);
void dfw_time_end(std::string);

const int primary_width   = PRIMARY_WIDTH;
const int secondary_width = SECONDARY_WIDTH;
const int chunk_hit_res = CHUNK_HIT_RES;

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
        
        anydsl::Array<int32_t> cpu_render_chunk_hit;
    };
    std::unordered_map<int32_t, DeviceData> devices;
    
    static thread_local anydsl::Array<float> cpu_primary;
    static thread_local anydsl::Array<float> cpu_secondary;
    anydsl::Array<int32_t> *cpu_record_chunk_hit;
    anydsl::Array<float> *cpu_primary_outgoing;
    anydsl::Array<float> *cpu_secondary_outgoing;
#ifdef ENABLE_EMBREE_DEVICE
    EmbreeDevice embree_device;
#endif

    anydsl::Array<float> host_pixels;
    anydsl::Array<float> tmp_pixels;
    anydsl::Array<float> host_primary;
    anydsl::Array<float> host_secondary;
    anydsl::Array<float> host_outgoing_primary;
    anydsl::Array<float> host_outgoing_secondary;

    //gpu & cpu
    anydsl::Array<int32_t> host_primary_num;
    anydsl::Array<int32_t> host_record_chunk_hit;
    anydsl::Array<int32_t> host_render_chunk_hit;
    
    anydsl::Array<int32_t> host_pass_record;

    size_t film_width;
    size_t film_height;

    std::mutex mtx;    
    bool first_write; 
    std::ofstream os;

    Interface(size_t width, size_t height)
        : film_width(width)
        , film_height(height)
    {
        
        host_pixels = anydsl::Array<float>(width * height * 3);
        tmp_pixels  = anydsl::Array<float>(width * height * 3);
       
        int chunk_size = dfw_chunk_size();
        host_primary_num = anydsl::Array<int32_t>(128);
        printf("light field size %d\n", chunk_hit_res * chunk_hit_res * 6 * chunk_size);
        
        int thread_num = dfw_thread_num();
        cpu_record_chunk_hit = new anydsl::Array<int32_t>[thread_num];
        cpu_primary_outgoing = new anydsl::Array<float>[thread_num];
        cpu_secondary_outgoing = new anydsl::Array<float>[thread_num];

        int out_size = dfw_out_stream_capacity(); 
        int out_capacity = (out_size & ~((1 << 5) - 1)) + 32; // round to 32
        int chunk_hit_size = chunk_hit_res * chunk_hit_res * 6 * chunk_size;
        int pass_record_size = chunk_size * (chunk_size + 3); // chunk_size * (send size of dst chunk, recv and no_hit)

        for(int i = 0; i < thread_num; i++) {
            cpu_record_chunk_hit[i] = anydsl::Array<int32_t>(chunk_hit_size);
            std::fill(cpu_record_chunk_hit[i].begin(), cpu_record_chunk_hit[i].end(), 0);
            cpu_primary_outgoing[i] = anydsl::Array<float>(primary_width * out_capacity);
            cpu_secondary_outgoing[i] = anydsl::Array<float>(secondary_width * out_capacity);
        }
        host_record_chunk_hit = anydsl::Array<int32_t>(chunk_hit_size);
        std::fill(host_record_chunk_hit.begin(), host_record_chunk_hit.end(), 0);
        
        host_render_chunk_hit = anydsl::Array<int32_t>(chunk_hit_size);
        std::fill(host_render_chunk_hit.begin(), host_render_chunk_hit.end(), 0);

        host_pass_record = anydsl::Array<int32_t>(pass_record_size);
        std::fill(host_pass_record.begin(), host_pass_record.end(), 0);

        os = std::ofstream("interface.log");
    }

    template <typename T>
    anydsl::Array<T>& resize_array(int32_t dev, anydsl::Array<T>& array, size_t size, size_t multiplier) {
        auto capacity = size; //(size & ~((1 << 5) - 1)) + 32; // round to 32
        if (array.size() < capacity) {
            auto n = capacity * multiplier;
            array = std::move(anydsl::Array<T>(dev, reinterpret_cast<T*>(anydsl_alloc(dev, sizeof(T) * n)), n));
        }
        return array;
    }

    anydsl::Array<int32_t>& render_chunk_hit_list(int32_t dev) {
        int size = chunk_hit_res * chunk_hit_res * 6 * dfw_chunk_size() * 2;
        return resize_array(dev, devices[dev].cpu_render_chunk_hit, size, 1);
    }
    

    //cpu 
    anydsl::Array<float>& cpu_primary_stream(size_t size) {
        return resize_array(0, cpu_primary, size, primary_width);
    }
    
    anydsl::Array<float>& cpu_secondary_stream(size_t size) {
        return resize_array(0, cpu_secondary, size, secondary_width);
    }
    
    //gpu
    void initial_gpu_host_data() {
        auto capacity = dfw_stream_store_capacity();
        resize_array(0, host_primary,  capacity, primary_width); 
        resize_array(0, host_outgoing_primary, capacity, primary_width); 
        resize_array(0, host_secondary,  capacity, secondary_width); 
        resize_array(0, host_outgoing_secondary, capacity, secondary_width); 
    } 
    
    anydsl::Array<float>& gpu_first_primary_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].first_primary, size, primary_width);
    }

    anydsl::Array<float>& gpu_second_primary_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].second_primary, size, primary_width);
    }
    
    anydsl::Array<float>& gpu_outgoing_primary_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].outgoing_primary, size, primary_width);
    }

    anydsl::Array<float>& gpu_first_secondary_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].first_secondary, size, secondary_width);
    }

    anydsl::Array<float>& gpu_second_secondary_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].second_secondary, size, secondary_width);
    }

    anydsl::Array<float>& gpu_outgoing_secondary_stream(int32_t dev, size_t size) {
        return resize_array(dev, devices[dev].outgoing_secondary, size, secondary_width);
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
        auto& bvh8_tri4 = devices[dev].bvh8_tri4;
        auto it = bvh8_tri4.find(filename);
        if (it != bvh8_tri4.end())
            return it->second;
        return bvh8_tri4[filename] = std::move(load_bvh<Node8, Tri4>(dev, filename));
    }
   
    void unload_chunk_data(int32_t chunk, int32_t dev) {
        printf("unload bvh\n");
        std::string chunk_path = "data/";
        chunk_path  += (chunk > 9 ? "0" : "00") + std::to_string(chunk) + "/";
        std::cout<<"unload chunk path "<<chunk_path<<"\n";
        devices[dev].bvh4_tri4.erase(chunk_path + "bvh_sse.bin");
        devices[dev].bvh8_tri4.erase(chunk_path + "bvh_avx.bin");
        devices[dev].bvh2_tri1.erase(chunk_path + "bvh_nvvm.bin");
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
                read_buffer(is, nodes);
                read_buffer(is, tris);
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
        auto& images = devices[dev].images;
        auto it = images.find(filename);
        if (it != images.end())
            return it->second;
        ImageRgba32 img;
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
    int gpu_outgoing_primary(int32_t dev) {
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
    primary.next_chk   = (int*)ptr + 9 * capacity;
    primary.org_chk   = (int*)ptr + 10 * capacity;
    primary.geom_id   = (int*)ptr + 11 * capacity;
    primary.prim_id   = (int*)ptr + 12 * capacity;
    primary.t         = ptr + 13 * capacity;
    primary.u         = ptr + 14 * capacity;
    primary.v         = ptr + 15 * capacity;
    primary.rnd       = (unsigned int*)ptr + 16 * capacity;
    primary.mis       = ptr + 17 * capacity;
    primary.contrib_r = ptr + 18 * capacity;
    primary.contrib_g = ptr + 19 * capacity;
    primary.contrib_b = ptr + 20 * capacity;
    primary.depth     = (int*)ptr + 21 * capacity;
    primary.size = 0;
    primary.capacity = dfw_stream_logic_capacity();
}

inline void get_secondary_stream(SecondaryStream& secondary, float* ptr, size_t capacity) {
    get_ray_stream(secondary.rays, ptr, capacity);
    secondary.next_chk = (int*)ptr + 9 * capacity;
    secondary.org_chk = (int*)ptr + 10 * capacity;
    secondary.pri_chk = (int*)ptr + 11 * capacity;
    secondary.prim_id = (int*)ptr + 12 * capacity;
    secondary.color_r = ptr + 13 * capacity;
    secondary.color_g = ptr + 14 * capacity;
    secondary.color_b = ptr + 15 * capacity;
    secondary.depth   = (int*)ptr + 16 * capacity;
    secondary.size = 0;
    secondary.capacity = dfw_stream_logic_capacity();
}

inline void get_out_ray_stream(OutRayStream& out, float* ptr, size_t capacity, size_t width) {
    out.fptr  = ptr;
    out.iptr  = (int*)ptr;
    out.uiptr = (unsigned int*) ptr;
    out.capacity = capacity;
    out.width = width;
    out.size = 0;
}

inline void get_chunk_hit(ChunkHit& chunk_hit, int* ptr) {
    chunk_hit.ctrb     = ptr;
    chunk_hit.res      = chunk_hit_res;
    chunk_hit.chk_cap  = chunk_hit_res * chunk_hit_res * 6;
    chunk_hit.chk_size = dfw_chunk_size();
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
    b.next_chk[dst_id]    = a.next_chk[src_id];
    b.org_chk[dst_id]    = a.org_chk[src_id];
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

//scene data
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
    dfw_time_start("interface => load bvh");
    auto& bvh = interface->load_bvh2_tri1(dev, file);
    *nodes = const_cast<Node2*>(bvh.nodes.data());
    *tris  = const_cast<Tri1*>(bvh.tris.data());
    dfw_time_end("interface => load bvh");
}

void rodent_load_bvh4_tri4(int32_t dev, const char* file, Node4** nodes, Tri4** tris) {
    dfw_time_start("interface => load bvh");
    auto& bvh = interface->load_bvh4_tri4(dev, file);
    *nodes = const_cast<Node4*>(bvh.nodes.data());
    *tris  = const_cast<Tri4*>(bvh.tris.data());
    dfw_time_end("interface => load bvh");
}

void rodent_load_bvh8_tri4(int32_t dev, const char* file, Node8** nodes, Tri4** tris) {
    dfw_time_start("interface => load bvh");
    auto& bvh = interface->load_bvh8_tri4(dev, file);
    printf("load bvh 8 4\n");
    *nodes = const_cast<Node8*>(bvh.nodes.data());
    *tris  = const_cast<Tri4*>(bvh.tris.data());
    dfw_time_end("interface => load bvh");
}

//cpu

void rodent_cpu_get_thread_data( PrimaryStream* primary, SecondaryStream* secondary
                               , OutRayStream * out_primary, OutRayStream * out_secondary
                               , ChunkHit *render_chunk_hit, ChunkHit *record_chunk_hit
                               , int dev, int tid, bool new_camera
                               ) 
{
    dfw_time_start("interface => prepare data");
    
    int rays_capacity = dfw_stream_store_capacity(); 
    auto& array_primary = interface->cpu_primary_stream(rays_capacity);
    get_primary_stream(*primary, array_primary.data(), array_primary.size() / primary_width);
    auto& array_secondary = interface->cpu_secondary_stream(rays_capacity);
    get_secondary_stream(*secondary, array_secondary.data(), array_secondary.size() / secondary_width);

    int out_size = dfw_out_stream_capacity(); 
    int out_capacity = (out_size & ~((1 << 5) - 1)) + 32; // round to 32
    
    auto& array_out_primary = interface->cpu_primary_outgoing[tid];
    get_out_ray_stream(*out_primary, array_out_primary.data(), out_size, primary_width);
    auto& array_out_secondary = interface->cpu_secondary_outgoing[tid];
    get_out_ray_stream(*out_secondary, array_out_secondary.data(), out_size, secondary_width);

    auto& array_record_chunk_hit = interface->cpu_record_chunk_hit[tid];
    get_chunk_hit(*record_chunk_hit, array_record_chunk_hit.data()); 
    // all thread share
    get_chunk_hit(*render_chunk_hit, interface->host_render_chunk_hit.data()); 

    if(new_camera) {
        interface->first_write = true;
        std::fill(array_record_chunk_hit.begin(), array_record_chunk_hit.end(), 0);
    }
    
    dfw_time_end("interface => prepare data");
}

//gpu
void rodent_gpu_get_data( PrimaryStream* primary, PrimaryStream* other_primary, SecondaryStream* secondary
                        , OutRayStream * out_primary, OutRayStream * out_secondary
                        , ChunkHit *render_chunk_hit, ChunkHit *record_chunk_hit
                        , int dev, int tid, bool new_camera) 
{
    int rays_capacity = dfw_stream_store_capacity(); 
    interface->initial_gpu_host_data(); 
    
    auto& array_first_primary = interface->gpu_first_primary_stream(dev, rays_capacity);
    get_primary_stream(*primary, array_first_primary.data(), array_first_primary.size() / primary_width);
    
    auto& array_second_primary = interface->gpu_second_primary_stream(dev, rays_capacity);
    get_primary_stream(*primary, array_second_primary.data(), array_second_primary.size() / primary_width);
        
    auto& array_first_secondary = interface->gpu_first_secondary_stream(dev, rays_capacity);
    get_secondary_stream(*secondary, array_first_secondary.data(), array_first_secondary.size() / secondary_width);

//    auto& array_second_secondary = interface->gpu__secondary_stream(dev, size);
//    get_secondary_stream(*secondary, array_first_secondary.data(), array_first_secondary.size() / secondary_width);

    int out_size = dfw_out_stream_capacity(); 
    int out_capacity = (out_size & ~((1 << 5) - 1)) + 32; // round to 32
    auto& array_out_primary = interface->gpu_outgoing_primary_stream(dev, out_capacity);
    get_out_ray_stream(*out_primary, array_out_primary.data(), out_size, primary_width);
    auto& array_out_secondary = interface->gpu_outgoing_secondary_stream(dev, out_capacity);
    get_out_ray_stream(*out_secondary, array_out_secondary.data(), out_size, secondary_width);
    
    auto& array_record_chunk_hit = interface->cpu_record_chunk_hit[0];
    if(new_camera) {
        interface->first_write = true;
        std::fill(array_record_chunk_hit.begin(), array_record_chunk_hit.end(), 0);
    }
    get_chunk_hit(*record_chunk_hit, array_record_chunk_hit.data()); 

    dfw_time_end("render prepare data");
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
void rodent_secondary_check(int32_t dev, int32_t tid, int primary_size, int buffer_size, int32_t print_mark, bool is_first_primary) {
    SecondaryStream secondary;
    printf("\n rank %d |%d rthread secondary size  %d\n", dfw_mpi_rank(), print_mark, primary_size);
    if (dev == -1) {
        auto& array = interface->cpu_secondary;
        get_secondary_stream(secondary, array.data(), array.size() / secondary_width);
    } 
    for(int i = 0; i < std::min(5, primary_size); i++){
        printf("ray id %d chunk %d tmin %f ray org %f %f %f dir %f %f %f\n", 
                    secondary.rays.id[i],    secondary.next_chk[i],    secondary.rays.tmin[i],
                    secondary.rays.org_x[i], secondary.rays.org_y[i], secondary.rays.org_z[i],
                    secondary.rays.dir_x[i], secondary.rays.dir_y[i], secondary.rays.dir_z[i]);
    }
    auto& array = interface->cpu_secondary_outgoing[tid];
    int *data = (int*) array.data();
    float *fptr = (float*) array.data();
    printf("buffer size%d\n", buffer_size);
    for(int i = 0; i < std::min(3, buffer_size); i++){
        printf("buffer id %d chunk %d tmin %f ray org %f %f %f dir %f %f %f\n", 
                    data[i*secondary_width], data[i*secondary_width + 9], data[i*secondary_width + 7],
                    fptr[i*secondary_width + 1], fptr[i*secondary_width + 2], fptr[i*secondary_width + 3], 
                    fptr[i*secondary_width + 4], fptr[i*secondary_width + 5], fptr[i*secondary_width + 6]) ;
    }
}

void rodent_first_primary_check(int32_t dev, int primary_size, int32_t print_mark, bool is_first_primary) {
    PrimaryStream primary;
    printf("rank %d %d rthread primary size %d\n", dfw_mpi_rank(), print_mark, primary_size);
    if (dev == -1) {
        auto& array = interface->cpu_primary;
        get_primary_stream(primary, array.data(), array.size() / primary_width);
    } else {
        int capacity = is_first_primary? interface->save_first_primary(dev) / primary_width : interface->save_second_primary(dev) / primary_width;
        auto& array = interface->host_primary;
        get_primary_stream(primary, array.data(), capacity);
    }
    for(int i = 0; i < 5; i++) {
        printf("ray id %d chunk %d tmin %f ray org %f %f %f dir %f %f %f geom %d prim %d depth %d\n", 
                    primary.rays.id[i+0], primary.next_chk[i+0], primary.rays.tmin[i],
                    primary.rays.org_x[i], primary.rays.org_y[i], primary.rays.org_z[i],
                    primary.rays.dir_x[i], primary.rays.dir_y[i], primary.rays.dir_z[i],
                    primary.geom_id[i], primary.prim_id[i], primary.depth[i]);
    }
}

//worker

int32_t rodent_cpu_thread_num() {
    return dfw_thread_num();
}

int32_t rodent_chunk_hit_resolution() {
    return CHUNK_HIT_RES;
}

inline int ctrb_cmp(const int &a, const int &b) {
    int res = 0;
    for(int i = 0; i < 8; i++) {
        int bit = i * 4;
        //res += (std::min(((a >> bit) & 0xF), ((b >> bit) & 0xF)) << bit);  
        res += (std::max(((a >> bit) & 0xF), ((b >> bit) & 0xF)) << bit);  
        //if((a >> bit) & 0xF == 0 )
        //    res += (((b >> bit) & 0xF) << bit);  
    }
    return res;
}

void rodent_save_chunk_hit(int32_t dev, int32_t tid) {
    interface->os<<"chunk hit from thread "<<tid<<"\n";
    if(dev == -1) {
        std::lock_guard <std::mutex> lock(interface->mtx);
        int* array = interface->cpu_record_chunk_hit[tid].data();
        int* host  = interface->host_record_chunk_hit.data();
        int  size  = chunk_hit_res * chunk_hit_res * 6 * dfw_chunk_size();
        if(interface->first_write) {
            memcpy(host, array, size * sizeof(int));
            interface->first_write = false;
            return;
        }
        for (int i = 0; i < size; i++){
            //interface->os<<host[i]<<" "<<array[i]<<" ";
            if(array[i] == 0) continue;
            if(host[i] == 0)
                host[i] = array[i];
            else 
                host[i] = ctrb_cmp(host[i], array[i]);
            //interface->os<<host[i]<<"\n";
        }
    } else {
        error("only cpu");
    }
}

int* rodent_get_chunk_hit() {
    //only save cpu data!!
    int thread_num = dfw_thread_num();
    for(int i = 0; i < thread_num; i++) { 
        rodent_save_chunk_hit(-1, i); 
        std::fill(interface->cpu_record_chunk_hit[i].begin(), interface->cpu_record_chunk_hit[i].end(), 0);
    }
    return interface->host_record_chunk_hit.data();
}

void rodent_worker_primary_send(int32_t dev, int32_t tid, int buffer_size) {
    if(buffer_size == 0) return;
    if(dev == -1) {
        auto& array = interface->cpu_primary_outgoing[tid];
        send_rays(array.data(), buffer_size, array.size() / primary_width, true);
    } else {
        int capacity = interface->gpu_outgoing_primary(dev) / primary_width;
        auto& array  = interface->host_outgoing_primary;
        send_rays(array.data(), buffer_size, capacity, true);
    }
}

void rodent_worker_secondary_send(int32_t dev, int32_t tid, int buffer_size) {
    if(buffer_size == 0) return;
    if(dev == -1) {
        auto& array = interface->cpu_secondary_outgoing[tid];
//            printf("array.size %d buffer size %d\n ", array.size() / secondary_width, buffer_size);
        send_rays(array.data(), buffer_size, array.size() / secondary_width, false);
    } else {
        int capacity = interface->save_buffer_secondary(dev) / secondary_width;
        auto& array  = interface->host_outgoing_secondary;
        send_rays(array.data(), buffer_size, capacity, false);
    }
}
int rodent_rays_export(int32_t dev, int32_t out_primary_size, int32_t out_secondary_size, bool isFirst, bool primary, int32_t tid) {
    
    rodent_worker_primary_send(dev, tid, out_primary_size);
    rodent_worker_secondary_send(dev, tid,  out_secondary_size);
    
    int size_new = 0;
    
    if(dev == -1) {
        if(primary) {
            float* array = interface->cpu_primary.data();
            size_new = recv_rays(&array, true, tid); 
        } else {
            float* array = interface->cpu_secondary.data();
            size_new = recv_rays(&array, false, tid); 
        }
    } else {
        if(primary) {
            float* array = interface->host_primary.data();
            size_new = recv_rays(&array, true, tid);
            if(size_new > 0)
                isFirst ? interface->load_first_primary(dev) : interface->load_second_primary(dev);
        } else {
            float* array = interface->host_secondary.data();
            size_new = recv_rays(&array, false, tid);
            if(size_new > 0) {
                isFirst ? interface->load_first_secondary(dev) : interface->load_second_secondary(dev);
            }
        }
    }
    return size_new;
} 

void rodent_memory_check(int32_t tag) {
    interface->os<<tag<<" "<<physical_memory_used_by_process()<<"mb \n";
}

void rodent_update_render_chunk_hit(int32_t* new_chunk_hit, int32_t size) {
    memcpy(interface->host_render_chunk_hit.data(), new_chunk_hit, size * sizeof(int));
}

void rodent_unload_chunk_data(int32_t chunk, int32_t dev) {   
    interface->unload_chunk_data(chunk, dev);
} 

int* rodent_get_render_chunk_hit_data() {
    interface->host_render_chunk_hit.data();
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
