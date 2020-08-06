#ifndef LOAD_OBJ_H
#define LOAD_OBJ_H

#include <vector>
#include <string>
#include <unordered_map>

#include "float3.h"
#include "color.h"
#include "FilePath.h"
#include "bbox.h"

namespace obj {

struct Index {
    int v, n, t;
};

struct Face {
    std::vector<Index> indices;
    int material;
};

struct Group {
    std::vector<Face> faces;
};

struct Object {
    std::vector<Group> groups;
};

struct Material {
    rgb ka;
    rgb kd;
    rgb ks;
    rgb ke;
    float ns;
    float ni;
    rgb tf;
    float tr;
    float d;
    int illum;
    std::string map_ka;
    std::string map_kd;
    std::string map_ks;
    std::string map_ke;
    std::string map_bump;
    std::string map_d;
};

struct File {
    std::vector<Object>      objects;
    std::vector<float3>      vertices;
    std::vector<float3>      normals;
    std::vector<float2>      texcoords;
    std::vector<std::string> materials;
    std::vector<std::string> mtl_libs;
    
    BBox bbox;
};

typedef std::unordered_map<std::string, Material> MaterialLib;

struct TriMesh {
    std::vector<float3>   vertices;
    std::vector<uint32_t> indices;
    std::vector<float3>   normals;
    std::vector<float3>   face_normals;
    std::vector<float2>   texcoords;
};

bool load_obj(const FilePath&, File&);
bool load_mtl(const FilePath&, MaterialLib&);
void write_obj(TriMesh *, int);
void mesh_add(TriMesh&, TriMesh&);
TriMesh compute_tri_mesh(const File&, const MaterialLib&, size_t, BBox&, bool);
bool chunk_division(File& file);

void compute_vertex_normals(const std::vector<uint32_t>& ,const std::vector<float3>&, std::vector<float3>&, size_t);
void compute_face_normals(const std::vector<uint32_t>&, const std::vector<float3>&, std::vector<float3>&, size_t);
} // namespace obj

#endif // LOAD_OBJ_H
