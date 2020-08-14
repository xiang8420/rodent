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
#include "obj.h"
#include "buffer.h"
#include "../distributed/decomposition.h" //domain split
#include "simplify.h"
#ifdef WIN32
#include <direct.h>
#define create_directory(d) _mkdir(d)
#else
#include <sys/stat.h>
#define create_directory(d) { umask(0); mkdir(d, 0777); }
#endif

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
    write_buffer_hetero(path,  "tri_vert.bin", tri_mesh.vertices,     target, true);
    write_buffer_hetero(path,  "tri_norm.bin", tri_mesh.normals,      target, true);
    write_buffer_hetero(path,  "face_nor.bin", tri_mesh.face_normals, target, true);
    write_buffer_hetero(path,  "text_cod.bin", tri_mesh.texcoords,    target, true);
    write_buffer(path + "tri_indx.bin",      tri_mesh.indices);
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
    write_buffer(of, nodes);
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

static size_t cleanup_obj(obj::File& obj_file, obj::MaterialLib& mtl_lib) {
    // Create a dummy material
    auto& dummy_mat = mtl_lib[""];
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

    // Check that all materials exist
    for (auto& mtl_name : obj_file.materials) {
        if (mtl_name != "" && !mtl_lib.count(mtl_name)) {
            warn("Missing material definition for '", mtl_name, "'. Replaced by dummy material.");
            mtl_name = "";
        }
    }

    // Remap identical materials (avoid duplicates)
    std::unordered_map<std::string, std::string> mtl_remap;
    for (size_t i = 0; i < obj_file.materials.size(); ++i) {
        auto& mtl1_name = obj_file.materials[i];
        auto& mtl1 = mtl_lib[mtl1_name];
        if (mtl_remap.count(mtl1_name) != 0)
            continue;
        for (size_t j = i + 1; j < obj_file.materials.size(); ++j) {
            auto& mtl2_name = obj_file.materials[j];
            auto& mtl2 = mtl_lib[mtl2_name];
            if (mtl1 == mtl2)
                mtl_remap.emplace(mtl2_name, mtl1_name);
        }
    }
    // Record unused materials
    std::unordered_set<std::string> used_mtls;
    for (auto& obj : obj_file.objects) {
        for (auto& group : obj.groups) {
            for (auto& face : group.faces) {
                auto mtl_name = obj_file.materials[face.material];
                auto it = mtl_remap.find(mtl_name);
                if (it != mtl_remap.end())
                    mtl_name = it->second;
                used_mtls.emplace(mtl_name);
            }
        }
    }

    // Remap indices/materials
    size_t num_complex = obj_file.materials.size();
    if (used_mtls.size() != obj_file.materials.size()) {
        std::vector<std::string> new_materials = obj_file.materials;
        new_materials.erase(std::remove_if(new_materials.begin(), new_materials.end(), [&] (auto& mtl_name) {
            return used_mtls.count(mtl_name) == 0;
        }), new_materials.end());
        // Put simple materials at the end
        num_complex = std::partition(new_materials.begin(), new_materials.end(), [&] (auto& mtl_name) {
            return !is_simple(mtl_lib[mtl_name]);
        }) - new_materials.begin();
        std::vector<uint32_t> mtl_id_remap;
        for (auto mtl_name : obj_file.materials) {
            auto it = mtl_remap.find(mtl_name);
            if (it != mtl_remap.end())
                mtl_name = it->second;

            auto new_mtl_index = std::find(new_materials.begin(), new_materials.end(), mtl_name) - new_materials.begin();
            mtl_id_remap.emplace_back(new_mtl_index);
        }
        for (auto& obj : obj_file.objects) {
            for (auto& group : obj.groups) {
                for (auto& face : group.faces) {
                    assert(face.material < mtl_id_remap.size());
                    face.material = mtl_id_remap[face.material];
                    assert(face.material < new_materials.size());
                }
            }
        }
        std::swap(obj_file.materials, new_materials);
        info("Removed ", new_materials.size() - obj_file.materials.size(), " unused/duplicate material(s)");
        info("The scene has ", num_complex, " complex material(s), and ", new_materials.size() - num_complex, " simple material(s)");
    }
    return num_complex;
}

void write_bvh(obj::TriMesh &tri_mesh, Target target, unsigned short &bvh_export, std::string &data_path, std::string &chunk_path) 
{
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
#ifdef ENABLE_EMBREE_BVH
            if (embree_bvh) build_embree_bvh<4>(tri_mesh, nodes, tris);
            else
#endif
            build_bvh<4, 4>(tri_mesh, nodes, tris);
            write_bvh_buffer(nodes, tris, chunk_path + "bvh_sse_.bin");
            bvh_export += 2;
        }
    } else {
        printf("avx \n");
        if(!(bvh_export&4)){
            create_directory((data_path + "cpu/").c_str());
            std::vector<typename BvhNTriM<8, 4>::Node> nodes;
            std::vector<typename BvhNTriM<8, 4>::Tri> tris;
#ifdef ENABLE_EMBREE_BVH
            if (embree_bvh) build_embree_bvh<8>(tri_mesh, nodes, tris);
            else
#endif
            build_bvh<8, 4>(tri_mesh, nodes, tris);
            write_bvh_buffer(nodes, tris, chunk_path + "bvh_avx_.bin");
            bvh_export += 4;
        }
    }
}

static bool convert_obj(const std::string& file_name, size_t dev_num, Target* target_list, size_t* dev_list, 
                        bool fusion, size_t max_path_len, size_t spp, bool embree_bvh, std::ostream& os, size_t grid_num) 
{
    info("Converting OBJ file '", file_name, "'");
    obj::File obj_file;
    obj::MaterialLib mtl_lib;
    FilePath scene_path(file_name);

    if (!obj::load_obj(scene_path, obj_file)) {
        error("Invalid OBJ file '", file_name, "'");
        return false;
    }

    for (auto lib_name : obj_file.mtl_libs) {
        auto mtl_name = scene_path.base_name() + "/" + lib_name;
        if (!obj::load_mtl(mtl_name, mtl_lib)) {
            error("Invalid MTL file '", mtl_name, "'");
            return false;
        }
    }

    size_t num_complex = cleanup_obj(obj_file, mtl_lib);
    size_t num_mats = obj_file.materials.size();

    std::unordered_map<std::string, size_t> images;
    bool has_map_ke = false;
    for (auto& pair : mtl_lib) {
        auto & mat = pair.second;
        if (mat.map_kd != "") images.emplace(mat.map_kd, images.size());
        if (mat.map_ks != "") images.emplace(mat.map_kd, images.size());
        if (mat.map_ke != "") images.emplace(mat.map_ke, images.size()), has_map_ke = true;
    }
    
    // Generate images
    std::vector<std::string> image_names(images.size());
    for (auto& pair : images)
        image_names[pair.second] = pair.first;

    os << "//------------------------------------------------------------------------------------\n"
       << "// Generated from '" << scene_path.file_name() << "' with the scene conversion tool\n"
       << "//------------------------------------------------------------------------------------\n\n";

    os << "struct Settings {\n"
       << "    eye: Vec3,\n"
       << "    dir: Vec3,\n"
       << "    up: Vec3,\n"
       << "    right: Vec3,\n"
       << "    width: f32,\n"
       << "    height: f32,\n"
       << "    image_region: Vec4_i32,\n"
       << "    spp: i32\n"
       << "};\n";

    os << "\nextern fn get_spp() -> i32 { " << spp << " } \n"
       << "extern fn get_dev_num() -> i32{ " << dev_num << " }\n"
       << "extern fn get_chunk_num() -> i32{ "<< grid_num << " }\n";


    os << "\nfn @make_scene(device: Device, settings: &Settings, file: File_path, chunk: i32, generateRays: bool)-> Scene {\n";

    os << "    let math     = device.intrinsics;\n"
       << "    let spp      = settings.spp;\n";

    os << "    // Camera\n"
       << "    let camera = make_perspective_camera(\n"
       << "        math,\n"
       << "        settings.eye,\n"
       << "        make_mat3x3(settings.right, settings.up, settings.dir),\n"
       << "        settings.width,\n"
       << "        settings.height,\n"
       << "        settings.image_region,\n"
       << "        spp,\n"
       << "        generateRays\n"
       << "    );\n";
    
    int gpu_num = 0; 
    for(int dev_id = 0; dev_id < dev_num; dev_id++){
        auto target = target_list[dev_id];
        assert(target != Target(0));
        if (target == Target::NVVM_STREAMING   ||
            target == Target::NVVM_MEGAKERNEL  ||
            target == Target::AMDGPU_STREAMING ||
            target == Target::AMDGPU_MEGAKERNEL){
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
	
  //  printf("memory used %d\n", physical_memory_used_by_process());
    // Setup camera

    MeshChunk chunk = MeshChunk(obj_file.bbox, grid_num);
    // Setup triangle mesh
    info("Generating triangle mesh for '", file_name, "'");
    const BBox &bbox = chunk.bbox; 
    os << "\n    // Triangle mesh\n"
       << "    let vertices     = device.load_buffer(file.vertices);\n"
       << "    let normals      = device.load_buffer(file.normals);\n"
       << "    let face_normals = device.load_buffer(file.face_normals);\n"
       << "    let texcoords    = device.load_buffer(file.texcoords);\n"
       << "    let indices      = device.load_buffer(file.indices);\n"
       << "    let tri_mesh     = TriMesh {\n"
       << "        vertices:     @ |i| vertices.load_vec3(i),\n"
       << "        normals:      @ |i| normals.load_vec3(i),\n"
       << "        face_normals: @ |i| face_normals.load_vec3(i),\n"
       << "        triangles:    @ |i| { let (i, j, k, _) = indices.load_int4(i); (i, j, k) },\n"
       << "        attrs:        @ |_| (false, @ |j| vec2_to_4(texcoords.load_vec2(j), 0.0f, 0.0f)),\n"
       << "        num_attrs:    1,\n"
       << "        num_tris:     27\n"  //wrong num tris
 //      << "        num_tris:     " << tri_mesh.indices.size() / 4 << "\n"
       << "    };\n"
       << "    let bvh = device.load_bvh(file.bvh);\n"
       << "    let node = bvh.node(0); \n"
       << "    let bbox = node.bbox(0);   \n";
    
    create_directory("data/");
    std::string data_path  = "data/"; 
    bool has_simple = false;
    
    // Lights global data
    size_t num_lights = 0;
    std::vector<rgb>    light_colors;
    std::vector<float3> light_verts;
    std::vector<float3> light_norms;
    std::vector<float>  light_areas;
    std::vector<obj::Material> light_mats;
   
    bool PreRendering = true; 
    obj::TriMesh simple_mesh; 
   
    int size = -1, rank = -1; 
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     
    for(int c = 0, n = chunk.size; c < n; c++) {
        auto tri_mesh = compute_tri_mesh(obj_file, mtl_lib, 0, chunk.list.at(c), false);
        printf("tri mesh num %ld \n", tri_mesh.indices.size() / 4);
        //if chunk is empty ??
		
		if(PreRendering) {
            int target_count = round((tri_mesh.indices.size() / 4) * 0.2);
            auto sub_mesh = Simplify::simplify_TriMesh(&tri_mesh, target_count, 7, true);
            obj::mesh_add(simple_mesh, sub_mesh);
      //      obj::write_obj(&simple_mesh, c);
        }
		std::string chunk_path = (c > 9 ? "data/0" : "data/00") + std::to_string(c) + "/";
        create_directory(chunk_path.c_str());
        num_complex = fusion ? num_complex : num_mats;
        has_simple |= num_complex < num_mats;
        // Simplify materials if necessary
        if (has_simple) {
            info("Simple materials will be fused");
            size_t num_tris = tri_mesh.indices.size() / 4;
            std::vector<rgb> simple_kd(num_tris);
            std::vector<rgb> simple_ks(num_tris);
            std::vector<float> simple_ns(num_tris);
            for (size_t i = 0, j = 0; i < tri_mesh.indices.size(); i += 4, j++) {
                auto& geom_id = tri_mesh.indices[i + 3];
                if (geom_id >= num_complex) {
                    auto& mat = mtl_lib[obj_file.materials[geom_id]];
                    assert(is_simple(mat));
                    simple_kd[j] = mat.kd;
                    simple_ks[j] = mat.ks;
                    simple_ns[j] = mat.ns;
                    geom_id = num_complex;
                } else if(geom_id == num_mats) {
                    geom_id = std::min(num_complex + 1, num_mats);
                } else {
                    simple_kd[j] = rgb(0.1f, 0.05f, 0.01f);
                    simple_ks[j] = rgb(0.1f, 0.05f, 0.01f);
                    simple_ns[j] = 1.0f;
                }
            }
            write_buffer_hetero(chunk_path, "simpl_kd.bin", simple_kd, padding_flag, true);
            write_buffer_hetero(chunk_path, "simpl_ks.bin", simple_ks, padding_flag, true);
            write_buffer(chunk_path + "simpl_ns.bin", simple_ns);
        }
        write_tri_mesh(chunk_path, tri_mesh, padding_flag);
        
        printf("after write tri mesh\n");
        unsigned short bvh_export = 0;
        for(int dev_id = 0; dev_id < dev_num; dev_id++) {
            Target target = target_list[dev_id];
            info("Generating BVH for '", file_name, "'");
            write_bvh(tri_mesh, target, bvh_export, data_path, chunk_path);
        }
        printf("light\n"); 
        //Local Light data light id
        std::vector<int> light_ids(tri_mesh.indices.size() / 4, 0);
        printf("tri mesh size %ld\n", tri_mesh.indices.size() / 4);
        for (size_t i = 0; i < tri_mesh.indices.size(); i += 4) {
            // Do not leave this array undefined, even if this triangle is not a light
            light_ids[i / 4] = 0;
            if(tri_mesh.indices[i + 3] >= obj_file.materials.size())
               continue; 

            auto& mtl_name = obj_file.materials[tri_mesh.indices[i + 3]];
            if (mtl_name == "")
                continue;
//            printf("mesh indices %d %d %d %d %d search mtl %s\n", i / 4, tri_mesh.indices[i], tri_mesh.indices[i+1], tri_mesh.indices[i+2], tri_mesh.indices[i+3], mtl_name.c_str());

            auto& mat = mtl_lib.find(mtl_name)->second;
            if (mat.ke == rgb(0.0f) && mat.map_ke == "")
                continue;
            printf("after check emission\n");

            auto& v0 = tri_mesh.vertices[tri_mesh.indices[i + 0]];
            auto& v1 = tri_mesh.vertices[tri_mesh.indices[i + 1]];
            auto& v2 = tri_mesh.vertices[tri_mesh.indices[i + 2]];

            int find_light = 0;
            for(find_light; find_light < num_lights; find_light++){
                int vert_id = find_light * 3;
            //    printf("vert id %d ", vert_id); 
                if(v0 == light_verts[vert_id] && v1 == light_verts[vert_id + 1] && v2 == light_verts[vert_id + 2])
                    break;
            } 
            printf("after find\n");
            light_ids[i / 4] = find_light;
            if(find_light != num_lights)
                continue;
            num_lights++;

            if (has_map_ke){
                os << "    let light" << i << " = make_triangle_light(\n"
                   << "        math,\n"
                   << "        make_vec3(" << light_verts[i * 3].x << "f, " << light_verts[i * 3].y << "f, " << light_verts[i * 3].z << "f),\n"
                   << "        make_vec3(" << light_verts[i * 3 + 1].x << "f, " << light_verts[i * 3 + 1].y << "f, " << light_verts[i * 3 + 1].z << "f),\n"
                   << "        make_vec3(" << light_verts[i * 3 + 2].x << "f, " << light_verts[i * 3 + 2].y << "f, " << light_verts[i * 3 + 2].z << "f),\n";
                if (light_mats[i].map_ke != "") {
                    os << "         make_texture(math, make_repeat_border(), make_bilinear_filter(), image_" << make_id(image_names[images[light_mats[i].map_ke]]) <<")\n";
                } else { 
                    os << "         make_color(" << light_colors[i].x << "f, " << light_colors[i].y << "f, " << light_colors[i].z << "f)\n";
                } 
                os << "        );\n";
            } else {
                auto n = cross(v1 - v0, v2 - v0);
                auto inv_area = 1.0f / (0.5f * length(n));
                n = normalize(n);
                light_verts.emplace_back(v0);
                light_verts.emplace_back(v1);
                light_verts.emplace_back(v2);
                light_norms.emplace_back(n);
                light_areas.emplace_back(inv_area);
                light_colors.emplace_back(mat.ke);
            }
        }
        printf("light over\n"); 
        write_buffer(chunk_path + "ligt_ids.bin", light_ids);
        printf("chunk over\n"); 
    }
    MPI_Finalize();
    
    ///write simple mesh
	if(PreRendering) {
        std::string simple_mesh_path = "data/999/";
        create_directory(simple_mesh_path.c_str());
        unsigned short t = 0; 
        write_tri_mesh(simple_mesh_path, simple_mesh, padding_flag);
        write_bvh(simple_mesh, Target::AVX2, t, data_path, simple_mesh_path);
        std::vector<int> light_ids(simple_mesh.indices.size() / 4, 0);
        for (size_t i = 0; i < simple_mesh.indices.size(); i += 4) {
            // Do not leave this array undefined, even if this triangle is not a light
            light_ids[i / 4] = 0;
            if(simple_mesh.indices[i + 3] >= obj_file.materials.size())
               continue; 
            auto& mtl_name = obj_file.materials[simple_mesh.indices[i + 3]];
            if (mtl_name == "")
                continue;
            auto& mat = mtl_lib.find(mtl_name)->second;
            if (mat.ke == rgb(0.0f) && mat.map_ke == "")
                continue;

            auto& v0 = simple_mesh.vertices[simple_mesh.indices[i + 0]];
            auto& v1 = simple_mesh.vertices[simple_mesh.indices[i + 1]];
            auto& v2 = simple_mesh.vertices[simple_mesh.indices[i + 2]];

            int find_light = 0;
            for(find_light; find_light < num_lights; find_light++){
                int vert_id = find_light * 3;
                if(v0 == light_verts[vert_id] && v1 == light_verts[vert_id + 1] && v2 == light_verts[vert_id + 2])
                    break;
            } 
            light_ids[i / 4] = find_light;
        }
        write_buffer(simple_mesh_path + "ligt_ids.bin", light_ids);
    }

    os << "\n    // Lights\n";
    if (has_map_ke || num_lights == 0){
        if (num_lights != 0) {
            os << "       let lights = @ |i| match i {\n";
            for (size_t i = 0; i < num_lights; ++i) {
                if (i == num_lights - 1)
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
        write_buffer_hetero(data_path, "ligt_are.bin",  light_areas,  padding_flag, true);
        write_buffer_hetero(data_path, "ligt_ver.bin",  light_verts,  padding_flag, true);
        write_buffer_hetero(data_path, "ligt_nor.bin",  light_norms,  padding_flag, true);
        write_buffer_hetero(data_path, "ligt_col.bin",  light_colors, padding_flag, true);
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
    }

    // Generate shaders
    info("Generating materials for '", file_name, "'");
    os << "\n    // Shaders\n";
    for (auto& mtl_name : obj_file.materials) {
        auto it = mtl_lib.find(mtl_name);
        assert(it != mtl_lib.end());

        auto& mat = it->second;
        // Stop at the first simple material (they have been moved to the end of the array)
        if (has_simple && is_simple(mat))
            break;

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

    if (has_simple) {
        os << "\n    // Simple materials data\n"
           << "    let simple_kd = device.load_buffer(\"data/simple_kd.bin\");\n"
           << "    let simple_ks = device.load_buffer(\"data/simple_ks.bin\");\n"
           << "    let simple_ns = device.load_buffer(\"data/simple_ns.bin\");\n";
    }

    // Generate geometries
    os << "\n    // Geometries\n"
       << "    let geometries = @ |i| match i {\n";
    for (uint32_t mat = 0; mat < num_complex; ++mat) {
        os << "    ";
        if (mat != num_complex - 1 || has_simple)
            os << mat;
        else
            os << "_";
        os << "        => make_tri_mesh_geometry(math, tri_mesh, shader_" << make_id(obj_file.materials[mat]) << "),\n";
    }
    if (has_simple)
        os << "    _ => make_tri_mesh_geometry(math, tri_mesh, @ |ray, hit, surf| {\n"
           << "        let kd = vec3_to_color(simple_kd.load_vec3(hit.prim_id));\n"
           << "        let ks = vec3_to_color(simple_ks.load_vec3(hit.prim_id));\n"
           << "        let ns = simple_ns.load_f32(hit.prim_id);\n"
           << "        let diffuse = make_diffuse_bsdf(math, surf, kd);\n"
           << "        let specular = make_phong_bsdf(math, surf, ks, ns);\n"
           << "        let lum_ks = color_luminance(ks);\n"
           << "        let lum_kd = color_luminance(kd);\n"
           << "        make_material(make_mix_bsdf(diffuse, specular, lum_ks / (lum_ks + lum_kd)))\n"
           << "    })\n";
    os << "    };\n";
    
    // Scene
    os << "\n       // Scene\n"
       << "    Scene {\n"
       << "        num_geometries: " << std::min(num_complex + 1, num_mats) << ",\n"
       << "        num_lights:     " << num_lights << ",\n"
       << "        geometries:     @ |i| geometries(i),\n"
       << "        lights:         @ |i| lights(i),\n"
       << "        camera:         camera,\n"
       << "        bvh:            bvh, \n"
       << "        bbox:           make_bbox(make_vec3(" << bbox.min.x << "f, " << bbox.min.y << "f, " << bbox.min.z << "f), (make_vec3(" 
                                                         << bbox.max.x << "f, " << bbox.max.y << "f, " << bbox.max.z << "f))),\n"
       << "        grid:           make_vec3("<< chunk.scale.x <<"f, "<<chunk.scale.y << "f, "<< chunk.scale.z <<"f),\n"     
       << "        chunk:          chunk \n"
       << "    }\n"
       << "}\n\n";
    
    os << "extern fn render(settings: &Settings, iter: i32, dev: i32, chunk: i32, generate_rays: bool) -> () { \n"   
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
        os << "        renderer(scene, device, iter);\n"
           << "        device.present();\n"
           << "    }\n";
    }
        
    os << "}\n\n";
     
    bool preprocess  = true; 

    if(preprocess) {
        os << "extern fn prerender(settings: &Settings) -> () {\n"   
           << "    let renderer = make_path_tracing_renderer(4 /*max_path_len*/, 1 /*spp*/); \n"
           << "    let device = make_prerender_avx2_device(); \n"
           << "    let scene    = make_scene(device, settings, make_file_path(1, 999), 0, true);\n"
           << "    renderer(scene, device, 0);\n"
           << "}\n";
    } else {
        os << "extern fn preprocess(settings: &Settings) -> () {}\n";
    }
    os << "extern fn get_bbox() -> &[f32] { &["<< bbox.min.x << "f, " << bbox.min.y << "f, " << bbox.min.z << "f, "
       << bbox.max.x << "f, " << bbox.max.y << "f, " << bbox.max.z << "f] }\n";
    os << "extern fn get_grid() -> &[f32] { &[" << chunk.scale.x <<"f, "<<chunk.scale.y << "f, "<< chunk.scale.z <<"f] }\n";

    info("Scene was converted successfully");
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
              << "           --fusion              Enables megakernel shader fusion (default: disabled)\n"
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
    bool fusion;

    size_t grid = 4;
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
                    if (!strcmp(argv[i], "--fusion")) {
                        fusion = true;
                        if (target[j] != Target::NVVM_MEGAKERNEL && target[j] != Target::AMDGPU_MEGAKERNEL) {
                            std::cerr << "Fusion is only available for megakernel targets. Aborting." << std::endl;
                            return 1;
                        }
                        ++i;
                    } else {
                        fusion = false;
                    }
                }
                --i; 
            } else if (!strcmp(argv[i], "-spp") || !strcmp(argv[i], "--samples-per-pixel")) {
                if (!check_option(i++, argc, argv)) return 1;
                spp = strtol(argv[i], NULL, 10);
            } else if(!strcmp(argv[i], "-g") || !strcmp(argv[i], "--grid")) { 
                if (!check_option(i++, argc, argv)) return 1;
                grid = strtol(argv[i], NULL, 10);
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
    if (!convert_obj(obj_file, dev_num, target, dev, fusion, max_path_len, spp, embree_bvh, of, grid))
        return 1;
    return 0;
    MPI_Finalize();
}
