#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <cstring>
#include <limits>
#include <mpi.h>
#include <algorithm>

#include "bvh.h"
#ifdef ENABLE_EMBREE_BVH
#include "embree_bvh.h"
#endif
#include "interface.h"
#include "../distributed/MeshChunk.h"
#include "obj.h"
#include "buffer.h"
#include "simplify.h"
#ifdef WIN32
#include <direct.h>
#define create_directory(d) _mkdir(d)
#else
#include <sys/stat.h>
#define create_directory(d) { umask(0); mkdir(d, 0777); }
#endif
#include <thread>
#include <mutex>
#include <condition_variable>
#include "../distributed/statistic.h"


enum class Target : uint32_t {
    GENERIC = 0,
    AVX2,
    AVX2_EMBREE,
    AVX,
    SSE42,
    ASIMD,
    NVVM_STREAMING,
    NVVM_MEGAKERNEL,
    AMDGPU_STREAMING,
    AMDGPU_MEGAKERNEL,
    INVALID
};

inline Target cpuid() {
    std::ifstream info(CPUINFO_PATH);
    if (!info) return Target(0);

    std::string line;
    std::vector<std::string> isa_list{
        "asimd", "sse4_2", "avx", "avx2"
    };
    std::unordered_set<std::string> detected;
    while (std::getline(info, line)) {
        for (auto isa : isa_list) {
            if (line.find(isa) != std::string::npos)
                detected.insert(isa);
        }
    }
    if (detected.count("avx2"))   return Target::AVX2;
    if (detected.count("avx"))    return Target::AVX;
    if (detected.count("sse4_2")) return Target::SSE42;
    if (detected.count("asimd"))  return Target::ASIMD;
    return Target::GENERIC;
}

void copy_file(const std::string& src, const std::string& dst) {
    for (size_t i = 0; i < dst.size();) {
        auto j = dst.find_first_of("/\\", i);
        if (j == std::string::npos)
            break;
        create_directory(dst.substr(0, j).c_str());
        i = j + 1;
    }
    info("Copying '", src, "' to '", dst, "'");
    std::ifstream is(src, std::ios::binary);
    std::ofstream os(dst, std::ios::binary);
    assert(is && os);
    os << is.rdbuf();
}

inline std::string make_id(const std::string& str) {
    auto id = str;
    std::transform(id.begin(), id.end(), id.begin(), [] (char c) {
        if (std::isspace(c) || !std::isalnum(c)) return '_';
        return c;
    });
    return id;
}

inline std::string fix_file(const std::string& str) {
    auto id = str;
    std::transform(id.begin(), id.end(), id.begin(), [] (char c) {
        if (c == '\\') return '/';
        return c;
    });
    return id;
}

inline bool ends_with(const std::string& str, const std::string& ext) {
    return str.rfind(ext) == str.length() - ext.length();
}

template <size_t N, size_t M>
struct BvhNTriM {};

template <>
struct BvhNTriM<8, 4> {
    using Node = Node8;
    using Tri  = Tri4;
};

template <>
struct BvhNTriM<4, 4> {
    using Node = Node4;
    using Tri  = Tri4;
};

template <>
struct BvhNTriM<2, 1> {
    using Node = Node2;
    using Tri  = Tri1;
};

template <size_t N, size_t M>
class BvhNTriMAdapter {
    struct CostFn {
        static float leaf_cost(int count, float area) {
            return count * area;
        }
        static float traversal_cost(float area) {
            return area;
        }
    };

    using BvhBuilder = SplitBvhBuilder<N, CostFn>;
    using Adapter    = BvhNTriMAdapter;
    using Node       = typename BvhNTriM<N, M>::Node;
    using Tri        = typename BvhNTriM<N, M>::Tri;

    std::vector<Node>& nodes_;
    std::vector<Tri>&  tris_;
    BvhBuilder         builder_;

public:
    BvhNTriMAdapter(std::vector<Node>& nodes, std::vector<Tri>& tris)
        : nodes_(nodes), tris_(tris)
    {}

    void build(const std::vector<uint32_t>& indices, const std::vector<::Tri>& tris) {
        builder_.build(tris, NodeWriter(*this), LeafWriter(*this, tris, indices), M / 2);
    }

#ifdef STATISTICS
    void print_stats() const override { builder_.print_stats(); }
#endif

private:
    struct NodeWriter {
        Adapter& adapter;

        NodeWriter(Adapter& adapter)
            : adapter(adapter)
        {}

        template <typename BBoxFn>
        int operator() (int parent, int child, const BBox& /*parent_bb*/, size_t count, BBoxFn bboxes) {
            auto& nodes = adapter.nodes_;

            size_t i = nodes.size();
            nodes.emplace_back();

            if (parent >= 0 && child >= 0) {
                assert(parent >= 0 && parent < nodes.size());
                assert(child >= 0 && child < N);
                nodes[parent].child[child] = i + 1;
            }

            assert(count >= 2 && count <= N);

            for (size_t j = 0; j < count; j++) {
                const BBox& bbox = bboxes(j);
                nodes[i].bounds[0][j] = bbox.min.x;
                nodes[i].bounds[2][j] = bbox.min.y;
                nodes[i].bounds[4][j] = bbox.min.z;

                nodes[i].bounds[1][j] = bbox.max.x;
                nodes[i].bounds[3][j] = bbox.max.y;
                nodes[i].bounds[5][j] = bbox.max.z;
            }

            for (size_t j = count; j < N; ++j) {
                nodes[i].bounds[0][j] = std::numeric_limits<float>::infinity();
                nodes[i].bounds[2][j] = std::numeric_limits<float>::infinity();
                nodes[i].bounds[4][j] = std::numeric_limits<float>::infinity();

                nodes[i].bounds[1][j] = -std::numeric_limits<float>::infinity();
                nodes[i].bounds[3][j] = -std::numeric_limits<float>::infinity();
                nodes[i].bounds[5][j] = -std::numeric_limits<float>::infinity();

                nodes[i].child[j] = 0;
            }

            return i;
        }
    };

    struct LeafWriter {
        Adapter& adapter;
        const std::vector<::Tri>& in_tris;
        const std::vector<uint32_t>& indices;

        LeafWriter(Adapter& adapter, const std::vector<::Tri>& in_tris, const std::vector<uint32_t>& indices)
            : adapter(adapter)
            , in_tris(in_tris)
            , indices(indices)
        {}

        template <typename RefFn>
        void operator() (int parent, int child, const BBox& /*leaf_bb*/, size_t ref_count, RefFn refs) {
            auto& nodes   = adapter.nodes_;
            auto& tris    = adapter.tris_;

            nodes[parent].child[child] = ~tris.size();

            // Group triangles by packets of M
            for (size_t i = 0; i < ref_count; i += M) {
                const size_t c = i + M <= ref_count ? M : ref_count - i;

                Tri tri;
                std::memset(&tri, 0, sizeof(Tri));
                for (size_t j = 0; j < c; j++) {
                    const int id = refs(i + j);
                    auto& in_tri = in_tris[id];
                    const float3 e1 = in_tri.v0 - in_tri.v1;
                    const float3 e2 = in_tri.v2 - in_tri.v0;
                    const float3 n = cross(e1, e2);
                    tri.v0[0][j] = in_tri.v0.x;
                    tri.v0[1][j] = in_tri.v0.y;
                    tri.v0[2][j] = in_tri.v0.z;

                    tri.e1[0][j] = e1.x;
                    tri.e1[1][j] = e1.y;
                    tri.e1[2][j] = e1.z;

                    tri.e2[0][j] = e2.x;
                    tri.e2[1][j] = e2.y;
                    tri.e2[2][j] = e2.z;

                    tri.n[0][j] = n.x;
                    tri.n[1][j] = n.y;
                    tri.n[2][j] = n.z;

                    tri.prim_id[j] = id;
                    tri.geom_id[j] = indices[id * 4 + 3];
                }

                for (size_t j = c; j < 4; j++)
                    tri.prim_id[j] = 0xFFFFFFFF;

                tris.emplace_back(tri);
            }
            assert(ref_count > 0);
            tris.back().prim_id[M - 1] |= 0x80000000;
        }
    };
};

template <>
class BvhNTriMAdapter<2, 1> {
    struct CostFn {
        static float leaf_cost(int count, float area) {
            return count * area;
        }
        static float traversal_cost(float area) {
            return area;
        }
    };

    using BvhBuilder = SplitBvhBuilder<2, CostFn>;
    using Adapter    = BvhNTriMAdapter;
    using Node       = Node2;
    using Tri        = Tri1;

    std::vector<Node>& nodes_;
    std::vector<Tri>&  tris_;
    BvhBuilder         builder_;

public:
    BvhNTriMAdapter(std::vector<Node>& nodes, std::vector<Tri>& tris)
        : nodes_(nodes), tris_(tris)
    {}

    void build(const std::vector<uint32_t>& indices, const std::vector<::Tri>& tris) {
        builder_.build(tris, NodeWriter(*this), LeafWriter(*this, tris, indices), 2);
    }

#ifdef STATISTICS
    void print_stats() const override { builder_.print_stats(); }
#endif

private:
    struct NodeWriter {
        Adapter& adapter;

        NodeWriter(Adapter& adapter)
            : adapter(adapter)
        {}

        template <typename BBoxFn>
        int operator() (int parent, int child, const BBox& /*parent_bb*/, size_t count, BBoxFn bboxes) {
            auto& nodes = adapter.nodes_;

            size_t i = nodes.size();
            nodes.emplace_back();

            if (parent >= 0 && child >= 0) {
                assert(parent >= 0 && parent < nodes.size());
                assert(child >= 0 && child < 2);
                nodes[parent].child[child] = i + 1;
            }

            assert(count >= 1 && count <= 2);

            const BBox& bbox1 = bboxes(0);
            nodes[i].bounds[0] = bbox1.min.x;
            nodes[i].bounds[2] = bbox1.min.y;
            nodes[i].bounds[4] = bbox1.min.z;
            nodes[i].bounds[1] = bbox1.max.x;
            nodes[i].bounds[3] = bbox1.max.y;
            nodes[i].bounds[5] = bbox1.max.z;

            if (count == 2) {
                const BBox& bbox2 = bboxes(1);
                nodes[i].bounds[ 6] = bbox2.min.x;
                nodes[i].bounds[ 8] = bbox2.min.y;
                nodes[i].bounds[10] = bbox2.min.z;
                nodes[i].bounds[ 7] = bbox2.max.x;
                nodes[i].bounds[ 9] = bbox2.max.y;
                nodes[i].bounds[11] = bbox2.max.z;
            } else {
                nodes[i].bounds[ 6] =  std::numeric_limits<float>::infinity();
                nodes[i].bounds[ 8] =  std::numeric_limits<float>::infinity();
                nodes[i].bounds[10] =  std::numeric_limits<float>::infinity();
                nodes[i].bounds[ 7] = -std::numeric_limits<float>::infinity();
                nodes[i].bounds[ 9] = -std::numeric_limits<float>::infinity();
                nodes[i].bounds[11] = -std::numeric_limits<float>::infinity();
            }

            return i;
        }
    };

    struct LeafWriter {
        Adapter& adapter;
        const std::vector<::Tri>& in_tris;
        const std::vector<uint32_t>& indices;

        LeafWriter(Adapter& adapter, const std::vector<::Tri>& in_tris, const std::vector<uint32_t>& indices)
            : adapter(adapter)
            , in_tris(in_tris)
            , indices(indices)
        {}

        template <typename RefFn>
        void operator() (int parent, int child, const BBox& /*leaf_bb*/, size_t ref_count, RefFn refs) {
            auto& nodes   = adapter.nodes_;
            auto& tris    = adapter.tris_;

            nodes[parent].child[child] = ~tris.size();

            for (int i = 0; i < ref_count; i++) {
                const int ref = refs(i);
                auto& tri = in_tris[ref];
                auto e1 = tri.v0 - tri.v1;
                auto e2 = tri.v2 - tri.v0;
                auto n = cross(e1, e2);
                int geom_id = indices[ref * 4 + 3];
                tris.emplace_back(Tri1 {
                    { tri.v0.x, tri.v0.y, tri.v0.z}, 0,
                    { e1.x, e1.y, e1.z}, geom_id,
                    { e2.x, e2.y, e2.z}, ref
                });
            }

            // Add sentinel
            tris.back().prim_id |= 0x80000000;
        }
    };
};

template <typename T>
static std::vector<uint8_t> pad_buffer(const std::vector<T>& elems, bool enable, size_t size) {
    std::vector<uint8_t> new_elems;
    if (!enable) {
        new_elems.resize(sizeof(T) * elems.size());
        memcpy(new_elems.data(), elems.data(), sizeof(T) * elems.size());
        return new_elems;
    }
    assert(size >= sizeof(T));
    new_elems.resize(size * elems.size(), 0);
    uint8_t* ptr = new_elems.data();
    for (auto& elem : elems) {
        memcpy(ptr, &elem, sizeof(T));
        ptr += size;
    }
    return new_elems;
}

template <typename Array>
static void write_buffer_hetero(std::string path, std::string file_name, const Array& array, unsigned short mask, bool padding) {
    if(mask & 1){
        create_directory((path + "cpu/").c_str());
        if(padding)
            write_buffer(path + "cpu/" + file_name, pad_buffer(array, false, sizeof(float) * 4));
        else
            write_buffer(path + "cpu/" + file_name, array);

    }
    if(mask & 2){
        create_directory((path + "gpu/").c_str());
        if(padding)
            write_buffer(path + "gpu/" + file_name, pad_buffer(array, true, sizeof(float) * 4));
        else 
            write_buffer(path + "gpu/" + file_name, array);
    }
}

static void write_tri_mesh(std::string path, const obj::TriMesh& tri_mesh, unsigned short target ) {
    printf("write tri mesh %s\n", path.c_str());
    write_buffer_hetero(path,  "normals.bin", tri_mesh.normals,      target, true);
    write_buffer_hetero(path,  "face_normals.bin", tri_mesh.face_normals, target, true);
    write_buffer_hetero(path,  "textcoords.bin", tri_mesh.texcoords,    target, true);
    write_buffer(path + "indices.bin",      tri_mesh.indices);
}


template <size_t N, size_t M>
static void build_bvh(const obj::TriMesh& tri_mesh,
                      std::vector<typename BvhNTriM<N, M>::Node>& nodes,
                      std::vector<typename BvhNTriM<N, M>::Tri>& tris) {
    BvhNTriMAdapter<N, M> adapter(nodes, tris);
    auto num_tris = tri_mesh.indices.size() / 4;
    std::vector<::Tri> in_tris(num_tris);
    for (size_t i = 0; i < num_tris; i++) {
        auto& v0 = tri_mesh.vertices[tri_mesh.indices[i * 4 + 0]];
        auto& v1 = tri_mesh.vertices[tri_mesh.indices[i * 4 + 1]];
        auto& v2 = tri_mesh.vertices[tri_mesh.indices[i * 4 + 2]];
        in_tris[i] = Tri(v0, v1, v2);
    }
    adapter.build(tri_mesh.indices, in_tris);
}

template <size_t N, size_t M>
static void build_muti_bvh(const obj::TriMesh& tri_mesh, obj::File& file, 
                      std::vector<typename BvhNTriM<N, M>::Node>& nodes,
                      std::vector<typename BvhNTriM<N, M>::Tri>& tris) {
    BvhNTriMAdapter<N, M> adapter(nodes, tris);
    
    auto num_tris = tri_mesh.indices.size() / 4;
    std::vector<::Tri> in_tris; //(num_tris);
    std::vector<uint32_t> indices; //(tri_mesh.indices.size());
    for (size_t i = 0; i < num_tris; i++) {
        auto& v0 = tri_mesh.vertices[tri_mesh.indices[i * 4 + 0]];
        auto& v1 = tri_mesh.vertices[tri_mesh.indices[i * 4 + 1]];
        auto& v2 = tri_mesh.vertices[tri_mesh.indices[i * 4 + 2]];
        in_tris.push_back(Tri(v0, v1, v2));
        for(size_t j = 0; j < 4; j++){
            uint32_t id = tri_mesh.indices.at(i * 4 + j); 
            indices.push_back(id);
        }
    }
    adapter.build(indices, in_tris);
}

template <typename Node, typename Tri>
static void write_bvh_buffer(std::vector<Node>& nodes, std::vector<Tri>& tris, std::string path) {
    std::ofstream of(path, std::ios::binary);
//    std::ofstream of("data/bvh.bin", std::ios::app | std::ios::binary);
    size_t node_size = sizeof(Node);
    size_t tri_size  = sizeof(Tri);
    of.write((char*)&node_size, sizeof(uint32_t));
    of.write((char*)&tri_size,  sizeof(uint32_t));
    printf("write nvh nodes %ld\n", nodes.size());
    
   
//    std::ofstream of_nodes("data/nodes", std::ios::binary);
//    size_t in_size = nodes.size() * sizeof(nodes[0]);
//    printf("nodes dat %d\n", in_size);
//    of_nodes.write((char*)&in_size,  sizeof(size_t));
//    of_nodes.write((char*)nodes.data(), in_size);
//
//    std::ofstream of_tris("data/tris", std::ios::binary);
//    in_size = tris.size() * sizeof(tris[0]);
//    printf("tris dat size %ld sizeof tris[0] %ld\n", tris.size(), sizeof(tris[0]));
//    of_tris.write((char*)&in_size,  sizeof(size_t));
//    of_tris.write((char*)tris.data(), in_size);

    write_buffer(of, nodes);
    printf("write bvh tri %ld\n", tris.size());
    write_buffer(of, tris );
    info("BVH with ", nodes.size(), " node(s), ", tris.size(), " tri(s)");
}

static bool operator == (const obj::Material& a, const obj::Material& b) {
    return a.ka == b.ka &&
           a.kd == b.kd &&
           a.ks == b.ks &&
           a.ke == b.ke &&
           a.ns == b.ns &&
           a.ni == b.ni &&
           a.tf == b.tf &&
           // ignored: a.tr == b.tr &&
           // ignored: a.d == b.d &&
           a.illum == b.illum &&
           // ignored: a.map_ka == b.map_ka &&
           a.map_kd == b.map_kd &&
           a.map_ks == b.map_ks &&
           a.map_ke == b.map_ke &&
           // ignored: a.map_bump == b.map_bump &&
           // ignored: a.map_d == b.map_d;
           true;
}

inline bool is_simple(const obj::Material& mat) {
    return mat.illum != 5 && mat.illum != 7 &&              // Must be diffuse
           mat.ke == rgb(0.0f) && mat.map_ke == "" &&       // Must not be emitting
           mat.map_kd == "" && mat.map_ks == "" &&          // Must not contain any texture
           (mat.kd != rgb(0.0f) || mat.ks != rgb(0.0f));    // Must not be completely black
}

static size_t cleanup_mtl(obj::MaterialLib& mtl_lib) {
    // Create a dummy material
    auto& dummy_mat = mtl_lib.map[""];
    dummy_mat.ka = rgb(0.0f);
    dummy_mat.kd = rgb(0.0f, 1.0f, 1.0f);
    dummy_mat.ks = rgb(0.0f);
    dummy_mat.ke = rgb(0.0f);
    dummy_mat.ns = 1.0f;
    dummy_mat.ni = 1.0f;
    dummy_mat.tf = rgb(0.0f);
    dummy_mat.tr = 1.0f;
    dummy_mat.d  = 1.0f;
    dummy_mat.illum = 2;
    dummy_mat.map_ka = "";
    dummy_mat.map_kd = "";
    dummy_mat.map_ks = "";
    dummy_mat.map_ke = "";
    dummy_mat.map_bump = "";
    dummy_mat.map_d = "";

    //set mtl_lib.list 
    std::cout<<"mtl_lib list"<<"\n ";
    std::unordered_map<std::string, obj::Material>::iterator iter ;   
    for(iter = mtl_lib.map.begin(); iter != mtl_lib.map.end(); iter++) {
        mtl_lib.list.emplace_back(iter->first);
    }
    std::cout<<mtl_lib.list.size()<<"\n";

    
    // Put simple materials at the end
    std::vector<std::string> new_materials = mtl_lib.list;
    size_t num_complex = mtl_lib.list.size();
    num_complex = std::partition(new_materials.begin(), new_materials.end(), [&] (auto& mtl_name) {
        return !is_simple(mtl_lib.map[mtl_name]);
    }) - new_materials.begin();
    std::swap(mtl_lib.list, new_materials);
    
    info("total ", new_materials.size(), "material(s)");
   

    //Delete duplicate
    new_materials.clear();
    for(int i = 0; i < mtl_lib.list.size(); i++) { 
        auto& mtl1_name = mtl_lib.list[i];
        auto& mtl1 = mtl_lib.map[mtl1_name]; 
        
        int find = 0;
        for (find = 0; find < new_materials.size(); find++) {
            auto& mtl2_name = new_materials[find];
            auto& mtl2 = mtl_lib.map[mtl2_name];
            if (mtl1 == mtl2) 
                break;
        }
        mtl_lib.ids.insert(std::make_pair(mtl1_name, find));
        if(find == new_materials.size()) 
            new_materials.emplace_back(mtl1_name);
    }
    std::swap(mtl_lib.list, new_materials);


    for (size_t i = 0; i < mtl_lib.list.size(); ++i) {
        auto& mtl1_name = mtl_lib.list[i];
        std::cout<<"mtl1 name "<<i<<" "<<mtl1_name<<" "<<"\n";
    }
    std::unordered_map<std::string, int>::iterator iter2 ;   
    for(iter2 = mtl_lib.ids.begin(); iter2 != mtl_lib.ids.end(); iter2++) {
        std::cout<<"mtl lib id "<<iter2->first<<" "<<iter2->second<<"\n";
    }

    return mtl_lib.list.size();
}

void write_bvh(obj::TriMesh &tri_mesh, Target target, unsigned short &bvh_export, std::string &data_path, std::string &chunk_path) 
{
    printf("write bvh\n");
    if (target == Target::NVVM_STREAMING   || target == Target::NVVM_MEGAKERNEL ||
        target == Target::AMDGPU_STREAMING || target == Target::AMDGPU_MEGAKERNEL) {
        printf("nvvm\n");
        if(!(bvh_export&1)){
            create_directory((data_path + "gpu/").c_str());
            std::vector<typename BvhNTriM<2, 1>::Node> nodes;
            std::vector<typename BvhNTriM<2, 1>::Tri> tris;
            build_bvh<2, 1>(tri_mesh, nodes, tris);
            printf("bvh built\n");
            write_bvh_buffer(nodes, tris, chunk_path + "bvh_nvvm.bin");
            bvh_export ++;
        }
    } else if (target == Target::GENERIC || target == Target::ASIMD || target == Target::SSE42) {
        printf("sse42 %d bvh_export %d \n", bvh_export, bvh_export&2);
        if(!(bvh_export&2)){
            create_directory((data_path + "cpu/").c_str());
            std::vector<typename BvhNTriM<4, 4>::Node> nodes;
            std::vector<typename BvhNTriM<4, 4>::Tri> tris;
            build_bvh<4, 4>(tri_mesh, nodes, tris);
            write_bvh_buffer(nodes, tris, chunk_path + "bvh_sse.bin");
            bvh_export += 2;
        }
    } else {
        printf("avx \n");
        if(!(bvh_export&4)){
            std::cout<<"path :" + data_path + "cpu/" + " " + chunk_path + "bvh_avx.bin";
            create_directory((data_path + "cpu/").c_str());
            std::vector<typename BvhNTriM<8, 4>::Node> nodes;
            std::vector<typename BvhNTriM<8, 4>::Tri> tris;
            statistics.start("build bvh");
            build_bvh<8, 4>(tri_mesh, nodes, tris);
            statistics.end("build bvh");
            statistics.start("write bvh");
            write_bvh_buffer(nodes, tris, chunk_path + "bvh_avx.bin");
            statistics.end("write bvh");
            bvh_export += 4;
        }
    }
}
template <typename T>
void mpi_gather(std::vector<T> &data, int *out_counts, int in_size, int global_size, int proc_rank, int proc_size) {
    
    printf("%d mpi gather global size %d\n", proc_rank, global_size);
    std::vector<T> recv_buffer(data.size());
    
    std::vector<int> rcounts(proc_size);
    std::vector<int> displs(proc_size);
   
    if(proc_rank == 0) { 
        recv_buffer.resize(global_size * in_size);
        for(int i = 0; i < proc_size; i++) {
            rcounts[i] = out_counts[i] * in_size; 
            printf("rcounts %d out counts %d", rcounts[i], out_counts[i]);
        }
        displs[0] = 0;
        for(int i = 1; i < proc_size; i++) {
            displs[i] = displs[i - 1] + rcounts[i - 1];
            printf("displs %d ", displs[i]);
        }
    }
    printf(" %d before gather %d \n", proc_rank, out_counts[proc_rank]);
    MPI_Gatherv((float*)data.data(), out_counts[proc_rank] * in_size, MPI_FLOAT, (float*)recv_buffer.data(), rcounts.data(), displs.data(), MPI_FLOAT, 0, MPI_COMM_WORLD);
    if(proc_rank == 0) {
        std::swap(recv_buffer, data);
    }
}

void mpi_gather_light( std::vector<int>   &light_ids
                    , std::vector<int>    &num_lights
                    , std::vector<int>    &num_tris
                    , std::vector<float3> &light_verts
                    , std::vector<float3> &light_norms
                    , std::vector<float>  &light_areas
                    , std::vector<rgb>    &light_colors
                    , int proc_rank
                    , int proc_size)
{
    if(proc_size > 1) {
        std::vector<int> reduce_num_lights(proc_size);
        MPI_Reduce(num_lights.data(), reduce_num_lights.data(), proc_size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
        
        std::vector<int> reduce_num_tris(proc_size);
        MPI_Reduce(num_tris.data(), reduce_num_tris.data(), proc_size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
        
        if(proc_rank == 0) {
            std::swap(reduce_num_lights, num_lights);
            std::swap(reduce_num_tris, num_tris);
        }

    }
    
    int global_light_size = 0;
    //Gather light data
    if(proc_rank == 0) { 
        printf("%d proc reduce light %d %d %d %d\n", proc_rank, num_lights[0], num_lights[1], num_lights[2], num_lights[3]);
        printf("%d proc reduce tri %d %d %d %d\n", proc_rank, num_tris[0], num_tris[1], num_tris[2], num_tris[3]);
        for(int i = 0; i < num_lights.size(); i++) {
            printf("check light ids chunk %d tris %d lights %d\n\n", i, num_tris[i], num_lights[i]);
        }
        printf("light sizes :");
        for(int i = 0; i < num_lights.size(); i++) {
            printf(" %d ", num_lights[i]); 
            global_light_size += num_lights[i];
        }
        printf("global light size %d\n", global_light_size);
    }
    
    if(proc_size > 1) {
        mpi_gather(light_ids, num_lights.data(), 1, global_light_size, proc_rank, proc_size); 
        mpi_gather(light_verts, num_lights.data(), 9, global_light_size, proc_rank, proc_size); 
        mpi_gather(light_norms, num_lights.data(), 3, global_light_size, proc_rank, proc_size); 
        mpi_gather(light_colors, num_lights.data(), 3, global_light_size, proc_rank, proc_size); 
        mpi_gather(light_areas, num_lights.data(), 1, global_light_size, proc_rank, proc_size); 
        
        if(proc_rank != 0) return ;
        
        for(int i = 0 ; i < light_ids.size(); i ++) 
            printf(" after light ids %d\n", light_ids[i]);

    }

    std::vector<float3> new_light_verts;
    std::vector<float3> new_light_norms;
    std::vector<rgb>    new_light_colors;
    std::vector<float>  new_light_areas;

    std::vector<int> light_ids_map(global_light_size * 2);
    for(int i = 0 /*light_verts.size() / 3*/; i < global_light_size; i++) {
        auto& v0 = light_verts[i * 3 + 0];
        auto& v1 = light_verts[i * 3 + 1];
        auto& v2 = light_verts[i * 3 + 2];
        int find_light = 0;
        for(find_light; find_light < new_light_areas.size(); find_light++) {
            int vert_id = find_light * 3;
            if(v0 == new_light_verts[vert_id] && 
               v1 == new_light_verts[vert_id + 1] && 
               v2 == new_light_verts[vert_id + 2])
                break;
        } 
        light_ids_map[i * 2] = light_ids[i];
        light_ids_map[i * 2 + 1] = find_light;
        if(find_light != new_light_areas.size()) {
            continue;
        }

        new_light_verts.emplace_back(v0);
        new_light_verts.emplace_back(v1);
        new_light_verts.emplace_back(v2);
        new_light_norms.emplace_back(light_norms[i]);
        new_light_areas.emplace_back(light_areas[i]);
        new_light_colors.emplace_back(light_colors[i]);
    }
    std::swap(new_light_verts,  light_verts);
    std::swap(new_light_norms,  light_norms);
    std::swap(new_light_colors, light_colors);
    std::swap(new_light_areas,  light_areas);

    int n = 0;
    printf("num lights %ld\n", num_lights.size());
    for(int i = 0; i < num_lights.size(); i++) {
        std::string light_ids_path = "data/";
        if (i < 100) light_ids_path += "0";
        if (i < 10) light_ids_path += "0";
        light_ids_path += std::to_string(i) + "/light_ids.bin";

        printf("light ids chunk %d tris %d lights %d\n\n", i, num_tris[i], num_lights[i]);
        std::vector<int> tri_light_ids(num_tris[i]);

        //set all light id 0
        for(int j = 0; j < num_lights[i]; j++) {
            int dst = (n + j) * 2;
            printf("light ids map dst %d %d %d", dst, light_ids_map[dst], light_ids_map[dst + 1]);
            tri_light_ids[light_ids_map[dst]] = light_ids_map[dst + 1];
        }
        n += num_lights[i];
        
        write_buffer(light_ids_path, tri_light_ids);
        printf("end light ids chunk %d %ld\n", i, tri_light_ids.size());
    }

    printf("%ld %ld %ld %ld\n", light_verts.size(), light_norms.size(), light_colors.size(), light_areas.size());
    printf("end mpi light gather\n");
}

bool convert_simple_mesh(const obj::ScenePath& obj_file_paths, obj::MaterialLib &mtl_lib, obj::Light &light, size_t dev_num, 
                        Target* target_list, size_t* dev_list, size_t chunk_size, int padding_flag) 
{
    //set simple mesh 
    obj::TriMesh simple_mesh; 
    //get xml/obj file path 
    BBox global_bbox;
    auto ticks = std::chrono::high_resolution_clock::now();
    if(obj_file_paths.simple.empty()) {
        printf("No simple mesh given, simplify origin mesh!\n");
        // Read all mtls 
        for( int j = 0; j < obj_file_paths.mesh.size(); j ++) {
            obj::File obj_file;
            std::string file_name =  obj_file_paths.mesh[j] + ".obj";
            if (!obj::load_obj(file_name, obj_file, mtl_lib)) {
                error("Invalid OBJ file '", file_name, "'");
                return false;
            }
            if(j == 0) {
                simple_mesh = Simplify::simplify(obj_file, mtl_lib, 7, false);
            } else {
                auto sub_mesh = Simplify::simplify(obj_file, mtl_lib, 7, false);
                obj::mesh_add(simple_mesh, sub_mesh);
            }
        }
    } else {
        printf("Load simple mesh\n");
        for(int j = 0; j < obj_file_paths.simple.size(); j ++) {
            obj::File obj_file;
            std::string file_name =  obj_file_paths.simple[j] + ".obj";
            if (!obj::load_obj(file_name, obj_file, mtl_lib)) {
                error("Invalid OBJ file '", file_name, "'");
                return false;
            }
            if(j == 0) {
                global_bbox = obj_file.bbox;
                simple_mesh = compute_tri_mesh(obj_file, mtl_lib, 0, global_bbox, false);
            } else {
                auto sub_mesh = compute_tri_mesh(obj_file, mtl_lib, 0, global_bbox, false);
                obj::mesh_add(simple_mesh, sub_mesh);
            }
        }
    }

    std::cout<<obj_file_paths.mesh[0]<<"\n";

    float elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - ticks).count();
    // only proc 0 process the left part
    printf("simplify time %f\n", elapsed_ms);

    ticks = std::chrono::high_resolution_clock::now();

    std::string data_path  = "data/"; 
    std::string simple_mesh_path = "data/";
    if (chunk_size < 100) simple_mesh_path += "0";
    if (chunk_size < 10) simple_mesh_path += "0";
    simple_mesh_path += std::to_string(chunk_size) + "/";

    create_directory(simple_mesh_path.c_str());
    unsigned short t = 0; 
    write_tri_mesh(simple_mesh_path, simple_mesh, padding_flag);
    write_bvh(simple_mesh, Target::AVX2, t, data_path, simple_mesh_path);
    std::vector<int> light_ids(simple_mesh.indices.size() / 4, 0);
    for (size_t i = 0; i < simple_mesh.indices.size(); i += 4) {
        // Do not leave this array undefined, even if this triangle is not a light
        light_ids[i / 4] = 0;
        auto& mtl_name = mtl_lib.list[simple_mesh.indices[i + 3]];
        if (mtl_name == "")
            continue;
        auto& mat = mtl_lib.map.find(mtl_name)->second;
        if (mat.ke == rgb(0.0f) && mat.map_ke == "")
            continue;

        auto& v0 = simple_mesh.vertices[simple_mesh.indices[i + 0]];
        auto& v1 = simple_mesh.vertices[simple_mesh.indices[i + 1]];
        auto& v2 = simple_mesh.vertices[simple_mesh.indices[i + 2]];

        int find_light = 0;
        for(find_light; find_light < light.areas.size(); find_light++){
            int vert_id = find_light * 3;
            if(v0 == light.verts[vert_id] && v1 == light.verts[vert_id + 1] && v2 == light.verts[vert_id + 2])
                break;
        }
        light_ids[i / 4] = find_light;
    }
    write_buffer(simple_mesh_path + "light_ids.bin", light_ids);
    
    elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - ticks).count();
    // only proc 0 process the left part
    printf("convert simplify time %f\n", elapsed_ms);

    //obj::write_obj(simple_mesh, mtl_lib, chunk_size); 
    return true;
}

static bool convert_obj(const std::string& file_name, size_t dev_num, Target* target_list, size_t* dev_list, 
                        size_t max_path_len, size_t spp, bool embree_bvh, std::ostream& os, size_t chunk_size) 
{
    info("Converting OBJ file '", file_name, "'");
   
    FilePath scene_path(file_name);
    std::cout<<scene_path.path()<<scene_path.base_name()<<scene_path.file_name()<<"\n";
    // set padding 
    int gpu_num = 0; 
    for(int dev_id = 0; dev_id < dev_num; dev_id++){
        auto target = target_list[dev_id];
        assert(target != Target(0));
        if (target == Target::NVVM_STREAMING   ||
            target == Target::NVVM_MEGAKERNEL  ||
            target == Target::AMDGPU_STREAMING ||
            target == Target::AMDGPU_MEGAKERNEL) {
            gpu_num++; 
        }
    }
    unsigned short padding_flag = 0;
    if(gpu_num > 0){
        padding_flag |= 2;
    }
    if(gpu_num != dev_num){
        padding_flag |= 1;
    }
    
    //get xml/obj file path 
    obj::ScenePath obj_file_paths(file_name);
    
    // Read all mtls 
    obj::MaterialLib mtl_lib;
    // .mtl must has the same name of obj
    for(int i = 0; i < obj_file_paths.mesh.size(); i ++) {
        auto mtl_name = obj_file_paths.mesh[i] + ".mtl";
        std::cout<<"read mtl "<<mtl_name <<" ";
        if (!obj::load_mtl(mtl_name, mtl_lib)) {
            error("Invalid MTL file '", mtl_name, "'");
            return false;
        }
    }
    for(int i = 0; i < obj_file_paths.simple.size(); i ++) {
        auto mtl_name = obj_file_paths.simple[i] + ".mtl";
        std::cout<<"read mtl "<<mtl_name <<" ";
        if (!obj::load_mtl(mtl_name, mtl_lib)) {
            error("Invalid MTL file '", mtl_name, "'");
            return false;
        }
    }
    std::cout<<"\n";
    
    size_t num_mats = cleanup_mtl(mtl_lib);
    
    // Generate images
    std::unordered_map<std::string, size_t> images;
    bool has_map_ke = false;
    for (auto& pair : mtl_lib.map) {
        auto & mat = pair.second;
        if (mat.map_kd != "") images.emplace(mat.map_kd, images.size());
        if (mat.map_ks != "") images.emplace(mat.map_kd, images.size());
        if (mat.map_ke != "") images.emplace(mat.map_ke, images.size()), has_map_ke = true;
    }

    std::vector<std::string> image_names(images.size());
    for (auto& pair : images)
        image_names[pair.second] = pair.first;

    // Lights global data
    obj::Light light;
     
    std::vector<int> num_lights(chunk_size);
    std::fill(num_lights.begin(), num_lights.end(), 0);

    std::vector<int> num_tris(chunk_size);
    std::fill(num_tris.begin(), num_tris.end(), 0);
    std::vector<int> tri_lights;

    // Record if triangles in a chunk is light, 

    //create data directory 
    //std::string data_path  = "data_grid" + std::to_string(chunk_size) + "/"; 
    std::string data_path  = "data/"; 
    std::cout<<data_path<<"\n";
    create_directory(data_path.c_str());
    create_directory("picture/");
    float elapsed_ms = 1;
    
    auto ticks = std::chrono::high_resolution_clock::now();
    
    int proc_rank = -1, proc_size = -1;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
   
    MeshChunk *chunks;
   // if(proc_size != chunk_size && proc_size != 1) 
   //     error("chunk size must equals to proc size\n");

    printf("rank %d size %d chunk size %d:\n", proc_rank, proc_size, chunk_size);
    for(int i = proc_rank; i < chunk_size; i+=proc_size) {
        int &chunk_num_lights = num_lights[i]; 
        int &chunk_num_tris = num_tris[i]; 
        
        obj::TriMesh tri_mesh;
        std::string chunk_path = data_path;
        if (i < 100) chunk_path += "0";
        if (i < 10) chunk_path += "0";
        chunk_path += std::to_string(i) + "/";

        create_directory(chunk_path.c_str());
        
        for(int j = 0; j < obj_file_paths.mesh.size(); j ++) {
            obj::File obj_file;
            std::string file_name =  obj_file_paths.mesh[j] + ".obj";
            statistics.start("load obj");
            if (!obj::load_obj(file_name, obj_file, mtl_lib)) {
                error("Invalid OBJ file '", file_name, "'");
                return false;
            }
            statistics.end("load obj");
            statistics.start("compute tri mesh");
            if(j == 0) {
                chunks = new MeshChunk(obj_file.bbox, chunk_size);
                tri_mesh = compute_tri_mesh(obj_file, mtl_lib, 0, chunks->list[i], true);
            } else {
                auto sub_mesh = compute_tri_mesh(obj_file, mtl_lib, 0, chunks->list[i], false);
                obj::mesh_add(tri_mesh, sub_mesh);
            }
            statistics.end("compute tri mesh");
            printf("rank %d chunk %d obj %d tri mesh num %ld \n", proc_rank, i, j, tri_mesh.indices.size() / 4);
            //if chunk is empty ??
        }
        
        printf("trimesh v %d tri %d nor %d face %d tex %d\n", tri_mesh.vertices.size()
                , tri_mesh.indices.size(), tri_mesh.normals.size(), tri_mesh.face_normals.size(), tri_mesh.texcoords.size());

        statistics.start("write tri mesh");
        write_tri_mesh(chunk_path, tri_mesh, padding_flag);
        statistics.end("write tri mesh");
        
        printf("after write tri mesh proc %d %d\n", proc_rank, i);
        
        //MPI_Barrier(MPI_COMM_WORLD);
        //build and write bvh
        printf("dev num %d \n", dev_num);
        unsigned short bvh_export = 0;
        for(int dev_id = 0; dev_id < dev_num; dev_id++) {
            printf("target %d \n", dev_id);
            Target target = target_list[dev_id];
            info("Generating BVH for '", file_name, "'");
            //obj::write_obj(tri_mesh, mtl_lib, i);
            write_bvh(tri_mesh, target, bvh_export, data_path, chunk_path);
        }
        printf("light\n"); 
        //MPI_Barrier(MPI_COMM_WORLD);
         
        //Local chunk Light data, light id
        
        printf("tri mesh size %ld\n", tri_mesh.indices.size() / 4);
        
        statistics.start("write light");
        for (size_t i = 0; i < tri_mesh.indices.size(); i += 4) {
            // Do not leave this array undefined, even if this triangle is not a light
            //if virtual portal continue
            //std::cout<<"tri idx "<<tri_mesh.indices[i]<<" "<<tri_mesh.indices[i + 1]<<" "<<tri_mesh.indices[i + 2]<< "mtl idx "<<tri_mesh.indices[i + 3]<<"\n";
            if(tri_mesh.indices[i + 3] >= num_mats)
               continue;
           
            auto& mtl_name = mtl_lib.list[tri_mesh.indices[i + 3]];
       //     std::cout<<"mtl "<<mtl_name<<"\n";
            if (mtl_name == "")
                continue;

            auto& mat = mtl_lib.map.find(mtl_name)->second;
            if (mat.ke == rgb(0.0f) && mat.map_ke == "")
                continue;

            auto& v0 = tri_mesh.vertices[tri_mesh.indices[i + 0]];
            auto& v1 = tri_mesh.vertices[tri_mesh.indices[i + 1]];
            auto& v2 = tri_mesh.vertices[tri_mesh.indices[i + 2]];
           
            printf("%d find %ld is light\n", proc_rank, i / 4); 
            tri_lights.emplace_back(i / 4);
            chunk_num_lights++;

            if (has_map_ke) {
            //    os << "    let light" << i << " = make_triangle_light(\n"
            //       << "        math,\n"
            //       << "        make_vec3(" << light.verts[i * 3].x << "f, " << light.verts[i * 3].y << "f, " << light.verts[i * 3].z << "f),\n"
            //       << "        make_vec3(" << light.verts[i * 3 + 1].x << "f, " << light.verts[i * 3 + 1].y << "f, " << light.verts[i * 3 + 1].z << "f),\n"
            //       << "        make_vec3(" << light.verts[i * 3 + 2].x << "f, " << light.verts[i * 3 + 2].y << "f, " << light.verts[i * 3 + 2].z << "f),\n";
            //    if (light.mats[i].map_ke != "") {
            //        os << "         make_texture(math, make_repeat_border(), make_bilinear_filter(), image_" << make_id(image_names[images[light.mats[i].map_ke]]) <<")\n";
            //    } else { 
            //        os << "         make_color(" << light.colors[i].x << "f, " << light.colors[i].y << "f, " << light.colors[i].z << "f)\n";
            //    } 
            //    os << "        );\n";
            } else {
                auto n = cross(v1 - v0, v2 - v0);
                auto inv_area = 1.0f / (0.5f * length(n));
                n = normalize(n);
                light.verts.emplace_back(v0);
                light.verts.emplace_back(v1);
                light.verts.emplace_back(v2);
                light.norms.emplace_back(n);
                light.areas.emplace_back(inv_area);
                light.colors.emplace_back(mat.ke);
            }
        }
        statistics.end("write light");
        //obj::write_obj(tri_mesh, mtl_lib, i);
        chunk_num_tris = tri_mesh.indices.size() / 4;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_gather_light(tri_lights, num_lights, num_tris, light.verts, light.norms, light.areas, light.colors, proc_rank, proc_size);
    printf("%d after mpi gather %ld\n", proc_rank, light.areas.size());
    
    int global_num_lights = light.areas.size();

    elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - ticks).count();
    // only proc 0 process the left part
    printf("proc rank %d process time %f\n", proc_rank, elapsed_ms);
    if(proc_rank != 0) 
        return true;

    float* verts_data = (float*)light.verts.data();
    float* norms_data = (float*)light.norms.data();
    float* areas_data = (float*)light.areas.data();
    float* colors_data = (float*)light.colors.data();
    
    printf("recv lisght verts %ld :\n", light.areas.size()); 
    for(int j = 0; j < light.areas.size(); j++) {
        for(int k = 0; k < 3; k ++) {
          //  printf("%d %f %f %f | ", j, verts_data[j * 9 + k * 3 + 0], verts_data[j * 9 + k * 3 + 1], verts_data[j * 9 + k * 3 + 2]); 
            printf("%f %f %f \n | normal %f | color %f| \n", 
                    verts_data[j * 9 + k * 3 + 0], verts_data[j * 9 + k * 3 + 1], verts_data[j * 9 + k * 3 + 2], 
                    norms_data[j * 3 + k], colors_data[j * 3 + k]);
        }
        printf("| area %f |\n", areas_data[j]); 
    }
    printf("\n");


    os << "//------------------------------------------------------------------------------------\n"
       << "// Generated from '" << file_name << "' with the scene conversion tool\n"
       << "//------------------------------------------------------------------------------------\n\n";

    os << "struct Settings {\n"
       << "    eye: Vec3,\n"
       << "    dir: Vec3,\n"
       << "    up: Vec3,\n"
       << "    right: Vec3,\n"
       << "    width: f32,\n"
       << "    height: f32,\n"
       << "    image_region: Vec4_i32,\n"
       << "    spp: i32,\n"
       << "    simple_trace: bool\n" 
       << "};\n";
    os << "\nextern fn get_spp() -> i32 { " << spp << " } \n"
       << "extern fn get_dev_num() -> i32{ " << dev_num << " }\n"
       << "extern fn get_chunk_num() -> i32{ "<< chunk_size << " }\n";
    os << "\nfn @make_scene(device: Device, settings: &Settings, file: File_path, chunk: i32, generateRays: bool)-> Scene {\n";
    os << "    let math     = device.intrinsics;\n"
       << "    // Camera\n"
       << "    let camera = make_perspective_camera(\n"
       << "        math,\n"
       << "        settings.eye,\n"
       << "        make_mat3x3(settings.right, settings.up, settings.dir),\n"
       << "        settings.width,\n"
       << "        settings.height,\n"
       << "        settings.image_region,\n"
       << "        settings.spp,\n"
       << "        generateRays\n"
       << "    );\n";
	
  //  printf("memory used %d\n", physical_memory_used_by_process());

    info("Generating triangle mesh for '", file_name, "'");
    const BBox &bbox = chunks->bbox; 
    os << "\n    // Triangle mesh\n"
       << "    let normals      = device.load_buffer(file.normals);\n"
       << "    let face_normals = device.load_buffer(file.face_normals);\n"
       << "    let texcoords    = device.load_buffer(file.texcoords);\n"
       << "    let indices      = device.load_buffer(file.indices);\n"
       << "    let tri_mesh     = TriMesh {\n"
       << "        normals:      @ |i| normals.load_vec3(i),\n"
       << "        face_normals: @ |i| face_normals.load_vec3(i),\n"
       << "        triangles:    @ |i| { let (i, j, k, _) = indices.load_int4(i); (i, j, k) },\n"
       << "        attrs:        @ |_| (false, @ |j| vec2_to_4(texcoords.load_vec2(j), 0.0f, 0.0f)),\n"
       << "        num_attrs:    1,\n"
       << "        num_tris:     27\n"  //wrong num tris
 //      << "        num_tris:     " << tri_mesh.indices.size() / 4 << "\n"
       << "    };\n"
       << "    let bvh = device.load_bvh(file.bvh);\n";
   


    os << "\n    // Lights\n";
    if (has_map_ke || global_num_lights == 0){
        if (global_num_lights != 0) {
            os << "       let lights = @ |i| match i {\n";
            for (size_t i = 0; i < global_num_lights; ++i) {
                if (i == global_num_lights - 1)
                    os << "                _ => light" << i << "\n";
                else
                    os << "                " << i << " => light" << i << ",\n";
            }
            os << "          };\n";
        } else {
            os << "        let lights = @ |_| make_point_light(math, make_vec3(0.0f, 0.0f, 0.0f), black);\n";
        }
    } else {
//        write_buffer("data/light_areas.bin",  light_areas);
        write_buffer_hetero(data_path, "light_areas.bin",  light.areas,  padding_flag, false);
        write_buffer_hetero(data_path, "light_verts.bin",  light.verts,  padding_flag, true);
        write_buffer_hetero(data_path, "light_normals.bin",  light.norms,  padding_flag, true);
        write_buffer_hetero(data_path, "light_colors.bin",  light.colors, padding_flag, true);
        os << "    let light_verts = device.load_buffer(file.light_verts);\n"
           << "    let light_areas = device.load_buffer(file.light_areas);\n"
           << "    let light_norms = device.load_buffer(file.light_norms);\n"
           << "    let light_colors = device.load_buffer(file.light_colors);\n"
           << "    let lights = @ |i| {\n"
           << "        make_precomputed_triangle_light(\n"
           << "            math,\n"
           << "            light_verts.load_vec3(i * 3 + 0),\n"
           << "            light_verts.load_vec3(i * 3 + 1),\n"
           << "            light_verts.load_vec3(i * 3 + 2),\n"
           << "            light_norms.load_vec3(i),\n"
           << "            light_areas.load_f32(i),\n"
           << "            vec3_to_color(light_colors.load_vec3(i))\n"
           << "        )\n"
           << "    };\n";
    }


    os << "\n    // Mapping from primitive to light source\n"
       << "    let light_ids = device.load_buffer(file.light_ids);\n";

    // Generate images
    info("Generating images for '", file_name, "'");
    os << "\n    // Images\n"
       << "    let dummy_image = make_image(@ |x, y| make_color(0.0f, 0.0f, 0.0f), 1, 1);\n";
    for (size_t i = 0; i < images.size(); i++) {
        auto name = fix_file(image_names[i]);
        copy_file(scene_path.base_name() + "/" + name, "data/" + name);
        os << "        let image_" << make_id(name) << " = ";
        if (ends_with(name, ".png")) {
            os << "    device.load_png(\"data/" << name << "\");\n";
        } else if (ends_with(name, ".tga")) {
            os << "    device.load_tga(\"data/" << name << "\");\n";
        } else if (ends_with(name, ".tiff")) {
            os << "    device.load_tga(\"data/" << name << "\");\n";
        } else if (ends_with(name, ".jpeg") || ends_with(name, ".jpg")) {
            os << "    device.load_jpg(\"data/" << name << "\");\n";
        } else {
            os << "    dummy_image; // Cannot determine image type for " << name << "\n";
        }
        printf("after copy\n");
    }

    // Generate shaders
    info("Generating materials for '", file_name, "'");
    os << "\n    // Shaders\n";
    for (auto& mtl_name : mtl_lib.list) {
        auto it = mtl_lib.map.find(mtl_name);
        assert(it != mtl_lib.map.end());

        auto& mat = it->second;
        // Stop at the first simple material (they have been moved to the end of the array)

        bool has_emission = mat.ke != rgb(0.0f) || mat.map_ke != "";
        os << "    let shader_" << make_id(mtl_name) << " : Shader = @ |ray, hit, surf| {\n";
        if (mat.illum == 5) {
            os << "        let bsdf = make_mirror_bsdf(math, surf, make_color(" << mat.ks.x << "f, " << mat.ks.y << "f, " << mat.ks.z << "f));\n";
        } else if (mat.illum == 7) {
            os << "        let bsdf = make_glass_bsdf(math, surf, 1.0f, " << mat.ni << "f, " << "make_color(" << mat.ks.x << "f, " << mat.ks.y << "f, " << mat.ks.z << "f), make_color(" << mat.tf.x << "f, " << mat.tf.y << "f, " << mat.tf.z << "f));\n";
        } else {
            bool has_diffuse  = mat.kd != rgb(0.0f) || mat.map_kd != "";
            bool has_specular = mat.ks != rgb(0.0f) || mat.map_ks != "";

            if (has_diffuse) {
                if (mat.map_kd != "") {
                    os << "        let diffuse_texture = make_texture(math, make_repeat_border(), make_bilinear_filter(), image_" << make_id(image_names[images[mat.map_kd]]) << ");\n";
                    os << "        let kd = diffuse_texture(vec4_to_2(surf.attr(0)));\n";
                } else {
                    os << "        let kd = make_color(" << mat.kd.x << "f, " << mat.kd.y << "f, " << mat.kd.z << "f);\n";
                }
                os << "        let diffuse = make_diffuse_bsdf(math, surf, kd);\n";
            }
            if (has_specular) {
                if (mat.map_ks != "") {
                    os << "       let specular_texture = make_texture(math, make_repeat_border(), make_bilinear_filter(), image_" << make_id(image_names[images[mat.map_ks]]) << ");\n";
                    os << "       let ks = specular_texture(vec4_to_2(surf.attr(0)));\n";
                } else {
                    os << "       let ks = make_color(" << mat.ks.x << "f, " << mat.ks.y << "f, " << mat.ks.z << "f);\n";
                }
                os << "        let ns = " << mat.ns << "f;\n";
                os << "        let specular = make_phong_bsdf(math, surf, ks, ns);\n";
            }
            os << "        let bsdf = ";
            if (has_diffuse && has_specular) {
                os << "{\n"
                   << "            let lum_ks = color_luminance(ks);\n"
                   << "            let lum_kd = color_luminance(kd);\n"
                   << "            let k = select(lum_ks + lum_kd == 0.0f, 0.0f, lum_ks / (lum_ks + lum_kd));\n"
                   << "            make_mix_bsdf(diffuse, specular, k)\n"
                   << "        };\n";
            } else if (has_diffuse || has_specular) {
                if (has_specular) os << "specular;\n";
                else              os << "diffuse;\n";
            } else {
                os << "make_black_bsdf();\n";
            }
        }
        if (has_emission) {
            os << "        make_emissive_material(surf, bsdf, lights(light_ids.load_i32(hit.prim_id)))\n";
        } else {
            os << "        make_material(bsdf)\n";
        }
        os << "    };\n";
    }

    // Generate geometries
    os << "\n    // Geometries\n"
       << "    let geometries = @ |i| match i {\n";
    for (uint32_t mat = 0; mat < num_mats; ++mat) {
        os << "    ";
        if (mat != num_mats - 1)
            os << mat;
        else
            os << "_";
        os << "        => make_tri_mesh_geometry(math, tri_mesh, shader_" << make_id(mtl_lib.list[mat]) << "),\n";
    }
    os << "    };\n";
    
    // Scene
    os << "\n       // Scene\n"
       << "    Scene {\n"
       << "        num_geometries: " << num_mats << ",\n"
       << "        num_lights:     " << global_num_lights << ",\n"
       << "        geometries:     @ |i| geometries(i),\n"
       << "        lights:         @ |i| lights(i),\n"
       << "        camera:         camera,\n"
       << "        bvh:            bvh, \n"
       << "        bbox:           make_bbox(make_vec3(" << bbox.min.x << "f, " << bbox.min.y << "f, " << bbox.min.z << "f), (make_vec3(" 
                                                         << bbox.max.x << "f, " << bbox.max.y << "f, " << bbox.max.z << "f))),\n"
       << "        chunk:          make_vec3("<< chunks->scale.x <<"f, "<<chunks->scale.y << "f, "<< chunks->scale.z <<"f),\n"     
       << "        cur_chk:        chunk, \n"
       << "        simple_trace:    settings.simple_trace\n"
       << "    }\n"
       << "}\n\n";
     
    os << "extern fn load_bvh_test() -> () { \n"   
       << "    let device = make_avx2_device(false); \n"
       << "    device.load_bvh(make_file_path(1, get_chunk_num()).bvh); \n"
       << "} \n";

    os << "extern fn render(settings: &Settings, iter: i32, dev: i32, chunk: i32, next_chunk: i32, generate_rays: bool) -> () { \n"   
       << "    let renderer = make_path_tracing_renderer( "<<max_path_len<<" /*path length*/, " << spp << " /*spp*/); \n";
    for(int dev_id = 0; dev_id < dev_num; dev_id++){                            
                                                                                
        Target target = target_list[dev_id];
        auto dev = dev_list[dev_id];
        if(dev_id == 0){
            os << "    if dev == 0 { \n"; 
        }else if (dev_id == dev_num){
            os << "    else { \n"; 
        } else {
            os << "    else if dev == "<< dev_id  <<" { \n"; 
        }
        switch (target) {
            case Target::GENERIC:           os << "        let device = make_cpu_default_device();\n";               break;
            case Target::AVX2:              os << "        let device = make_avx2_device(false);\n";     break;
            case Target::AVX2_EMBREE:       os << "        let device = make_avx2_device(true);\n";      break;
            case Target::AVX:               os << "        let device = make_avx_device();\n";             break;
            case Target::SSE42:             os << "        let device = make_sse42_device();\n";                     break;
            case Target::ASIMD:             os << "        let device = make_asimd_device();\n";                     break;
            case Target::NVVM_STREAMING:    os << "        let device = make_nvvm_device(" << dev <<", true);\n";    break;
            case Target::NVVM_MEGAKERNEL:   os << "        let device = make_nvvm_device(" << dev <<", false);\n";   break;
            case Target::AMDGPU_STREAMING:  os << "        let device = make_amdgpu_device(" << dev <<", true);\n";  break;
            case Target::AMDGPU_MEGAKERNEL: os << "        let device = make_amdgpu_device(" << dev <<", false);\n"; break;
            default:
                assert(false);
                break;
        }
        if (target == Target::NVVM_STREAMING   || target == Target::NVVM_MEGAKERNEL ||
            target == Target::AMDGPU_STREAMING || target == Target::AMDGPU_MEGAKERNEL) {
            os << "        let scene  = make_scene(device, settings, make_file_path(0, chunk), chunk, generate_rays);\n";    
        } else if (target == Target::GENERIC || target == Target::ASIMD || target == Target::SSE42) {
            os << "        let scene  = make_scene(device, settings, make_file_path(2, chunk), chunk, generate_rays);\n";    
        } else {
            os << "        let scene  = make_scene(device, settings, make_file_path(1, chunk), chunk, generate_rays);\n";    
        }
        os << "        for i in parallel(2, 0, 2) {\n"
           << "            if i == 0 && next_chunk >= 0 { device.load_bvh(make_file_path(1, next_chunk).bvh); } \n"
           << "            if i == 1 { renderer(scene, device, iter);} \n"
           << "        }\n"
           << "        //device.present();\n"
           << "    }\n";
    }
        
    os << "}\n\n";
     
    os << "extern fn get_bbox() -> &[f32] { &["<< bbox.min.x << "f, " << bbox.min.y << "f, " << bbox.min.z << "f, "
       << bbox.max.x << "f, " << bbox.max.y << "f, " << bbox.max.z << "f] }\n";
    os << "extern fn get_chunk() -> &[f32] { &[" << chunks->scale.x <<"f, "<<chunks->scale.y << "f, "<< chunks->scale.z <<"f] }\n";

    bool export_simplify_mesh = SIMPLE_TRACE;
    // if use simple mesh
    if (export_simplify_mesh)
        convert_simple_mesh(obj_file_paths, mtl_lib, light, dev_num, target_list, dev_list, chunk_size, padding_flag); 

    info("Scene was converted successfully");
    statistics.print(0, 0);
    return true;
}

static void usage() {
    std::cout << "converter [options] file\n"
              << "Available options:\n"
              << "    -h     --help                Shows this message\n"
              << "    -t     --target              Sets the target platform (default: autodetect CPU)\n"
              << "    -d     --device              Sets the device to use on the selected platform (default: 0)\n"
              << "           --max-path-len        Sets the maximum path length (default: 64)\n"
              << "    -spp   --samples-per-pixel   Sets the number of samples per pixel (default: 4)\n"
#ifdef ENABLE_EMBREE_BVH
              << "           --embree-bvh          Use Embree to build the BVH (default: disabled)\n"
#endif
              << "Available target:\n"
              << "    generic, sse42, avx, avx2, avx2-embree, asimd,\n"
              << "    nvvm = nvvm-streaming, nvvm-megakernel,\n"
              << "    amdgpu = amdgpu-streaming, amdgpu-megakernel\n"
              << std::flush;
}

static bool check_option(int i, int argc, char** argv) {
    if (i + 1 >= argc) {
        std::cerr << "Missing argument for '" << argv[i] << "'. Aborting." << std::endl;
        return false;
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc <= 1) {
        std::cerr << "Not enough arguments. Run with --help to get a list of options." << std::endl;
        return 1;
    }
    std::string obj_file;
    
    size_t dev_num = 1;
    size_t dev[5];   // 5 device max
    Target target[5];

    size_t chunk = 4;
    size_t spp = 4;
    size_t max_path_len = 64;
    bool embree_bvh = false;
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
                usage();
                return 0;
            } else if (!strcmp(argv[i], "--dev-num")){
                if (!check_option(i++, argc, argv)) return 1;
                dev_num = strtoul(argv[i++], NULL, 10); 
                for(int j = 0; j < dev_num; j++){
                    if (!strcmp(argv[i], "-t") || !strcmp(argv[i], "--target")) {
                         if (!check_option(i++, argc, argv)) return 1;
                         if (!strcmp(argv[i], "sse42"))
                             target[j] = Target::SSE42;
                         else if (!strcmp(argv[i], "avx"))
                             target[j] = Target::AVX;
                         else if (!strcmp(argv[i], "avx2"))
                             target[j] = Target::AVX2;
                         else if (!strcmp(argv[i], "avx2-embree"))
                             target[j] = Target::AVX2_EMBREE;
                         else if (!strcmp(argv[i], "asimd"))
                             target[j] = Target::ASIMD;
                         else if (!strcmp(argv[i], "nvvm") || !strcmp(argv[i], "nvvm-streaming"))
                             target[j] = Target::NVVM_STREAMING;
                         else if (!strcmp(argv[i], "nvvm-megakernel"))
                             target[j] = Target::NVVM_MEGAKERNEL;
                         else if (!strcmp(argv[i], "amdgpu") || !strcmp(argv[i], "amdgpu-streaming"))
                             target[j] = Target::AMDGPU_STREAMING;
                         else if (!strcmp(argv[i], "amdgpu-megakernel"))
                             target[j] = Target::AMDGPU_MEGAKERNEL;
                         else if (!strcmp(argv[i], "generic"))
                             target[j] = Target::GENERIC;
                         else {
                             std::cerr << "Unknown target '" << argv[i] << "'. Aborting." << std::endl;
                             return 1;
                         }
                         ++i; 
                    } else {
                        target[j] = j==0?cpuid():target[j - 1];
                        if (target[j] == Target::GENERIC)
                            warn("No vector instruction set detected. Select the target platform manually to improve performance.");
                    }
                    if(!strcmp(argv[i], "-d") || !strcmp(argv[i], "--device")) { 
                        if (!check_option(i++, argc, argv)) return 1;
                        dev[j] = strtoul(argv[i++], NULL, 10);
                    } 
                    else {
                        dev[j] = j==0?0:dev[j - 1];
                    }
                }
                --i; 
            } else if (!strcmp(argv[i], "-spp") || !strcmp(argv[i], "--samples-per-pixel")) {
                if (!check_option(i++, argc, argv)) return 1;
                spp = strtol(argv[i], NULL, 10);
            } else if(!strcmp(argv[i], "-c") || !strcmp(argv[i], "--chunk")) { 
                if (!check_option(i++, argc, argv)) return 1;
                chunk = strtol(argv[i], NULL, 10);
            } else if (!strcmp(argv[i], "--max-path-len")) {
                if (!check_option(i++, argc, argv)) return 1;
                max_path_len = strtol(argv[i], NULL, 10);
#ifdef ENABLE_EMBREE_BVH
            } else if (!strcmp(argv[i], "--embree-bvh")) {
                embree_bvh = true;
#endif
            } else {
                std::cerr << "Unknown option '" << argv[i] << "'. Aborting." << std::endl;
                return 1;
            }
        } else {
            if (obj_file != "") {
                std::cerr << "Only one OBJ file can be converted. Aborting." << std::endl;
                return 1;
            }
            obj_file = argv[i];
        }
    }


    if (obj_file == "") {
        std::cerr << "Please specify an OBJ file to convert. Aborting." << std::endl;
        return 1;
    }

    std::ofstream of("main.impala");
    convert_obj(obj_file, dev_num, target, dev, max_path_len, spp, embree_bvh, of, chunk);
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("%d finish\n", rank);
    return 0;
}
