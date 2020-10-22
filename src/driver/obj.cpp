#include <fstream>
#include <iostream> 
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <iomanip>

#include "common.h"
#include "obj.h"
#include "../distributed/decomposition.h"
namespace obj {

struct TriIdx {
    int32_t v0, v1, v2, m;
    TriIdx(int32_t v0, int32_t v1, int32_t v2, int32_t m)
        : v0(v0), v1(v1), v2(v2), m(m)
    {}
};

struct HashIndex {
    size_t operator () (const obj::Index& i) const {
        unsigned h = 0, g;

        h = (h << 4) + i.v;
        g = h & 0xF0000000;
        h = g ? (h ^ (g >> 24)) : h;
        h &= ~g;

        h = (h << 4) + i.t;
        g = h & 0xF0000000;
        h = g ? (h ^ (g >> 24)) : h;
        h &= ~g;

        h = (h << 4) + i.n;
        g = h & 0xF0000000;
        h = g ? (h ^ (g >> 24)) : h;
        h &= ~g;

        return h;
    }
};

struct CompareIndex {
    bool operator () (const obj::Index& a, const obj::Index& b) const {
        return a.v == b.v && a.t == b.t && a.n == b.n;
    }
};

inline void remove_eol(char* ptr) {
    int i = 0;
    while (ptr[i]) i++;
    i--;
    while (i > 0 && std::isspace(ptr[i])) {
        ptr[i] = '\0';
        i--;
    }
}

inline char* strip_text(char* ptr) {
    while (*ptr && !std::isspace(*ptr)) { ptr++; }
    return ptr;
}

inline char* strip_spaces(char* ptr) {
    while (std::isspace(*ptr)) { ptr++; }
    return ptr;
}

inline int find(char * ptr, char t) {
    int i = 0;
    while(ptr[i] != t && ptr[i] != ' ' && ptr[i] != '\0') {
        i++;
    }
    return i;
}

inline bool read_index(char** ptr, obj::Index& idx) {
    char* base = *ptr;

    // Detect end of line (negative indices are supported) 
    base = strip_spaces(base);
    if (!std::isdigit(*base) && *base != '-') return false;

    idx.v = 0;
    idx.t = 0;
    idx.n = 0;

    idx.v = std::strtol(base, &base, 10);

    base = strip_spaces(base);

    if (*base == '/') {
        base++;

        // Handle the case when there is no texture coordinate
        if (*base != '/') {
            idx.t = std::strtol(base, &base, 10);
        }

        base = strip_spaces(base);

        if (*base == '/') {
            base++;
            idx.n = std::strtol(base, &base, 10);
        }
    }

    *ptr = base;

    return true;
}

static bool parse_obj(std::istream& stream, obj::File& file, obj::MaterialLib& mtl_lib) {
    // Add an empty object to the scene
    int cur_object = 0;
    file.objects.emplace_back();

    // Add an empty group to this object
    int cur_group = 0;
    file.objects[0].groups.emplace_back();

    // 
    int cur_mtl = 0;
    
    // Add dummy vertex, normal, and texcoord
    file.vertices.emplace_back();
    file.normals.emplace_back();
    file.texcoords.emplace_back();

    int err_count = 0, cur_line = 0;
    const int max_line = 1024;
    char line[max_line];

    // global bounding box    
    file.bbox = BBox::empty(); 

    while (stream.getline(line, max_line)) {
        cur_line++;

        // Strip spaces
        char* ptr = strip_spaces(line);

        // Skip comments and empty lines
        if (*ptr == '\0' || *ptr == '#')
            continue;

        remove_eol(ptr);

        // Test each command in turn, the most frequent first
        if (*ptr == 'v') {
            switch (ptr[1]) {
                case ' ':
                case '\t':
                    {
                        float3 v;
                        v.x = std::strtof(ptr + 1, &ptr);
                        v.y = std::strtof(ptr, &ptr);
                        v.z = std::strtof(ptr, &ptr);
                        file.vertices.push_back(v);
                        file.bbox.extend(v);
                    }
                    break;
                case 'n':
                    {
                        float3 n;
                        n.x = std::strtof(ptr + 2, &ptr);
                        n.y = std::strtof(ptr, &ptr);
                        n.z = std::strtof(ptr, &ptr);
                        file.normals.push_back(n);
                    }
                    break;
                case 't':
                    {
                        float2 t;
                        t.x = std::strtof(ptr + 2, &ptr);
                        t.y = std::strtof(ptr, &ptr);
                        file.texcoords.push_back(t);
                    }
                    break;
                default:
                    error("Invalid vertex (line ", cur_line, ").");
                    err_count++;
                    break;
            }
        } else if (*ptr == 'f' && std::isspace(ptr[1])) {
            obj::Face f;

            f.material = cur_mtl;
            bool valid = true;
            ptr += 2;
            while (valid) {
                obj::Index index;
                valid = read_index(&ptr, index);
                if (valid)
                    f.indices.push_back(index);
            }

            if (f.indices.size() < 3) {
                error("Invalid face (line ", cur_line, ").");
                err_count++;
            } else {
                // Convert relative indices to absolute
                for (size_t i = 0; i < f.indices.size(); i++) {
                    f.indices[i].v = (f.indices[i].v < 0) ? file.vertices.size()  + f.indices[i].v : f.indices[i].v;
                    f.indices[i].t = (f.indices[i].t < 0) ? file.texcoords.size() + f.indices[i].t : f.indices[i].t;
                    f.indices[i].n = (f.indices[i].n < 0) ? file.normals.size()   + f.indices[i].n : f.indices[i].n;
                }

                // Check if the indices are valid or not
                valid = true;
                for (size_t i = 0; i < f.indices.size(); i++) {
                    if (f.indices[i].v <= 0 || f.indices[i].t < 0 || f.indices[i].n < 0) {
                        valid = false;
                        break;
                    }
                }

                if (valid) {
                    file.objects[cur_object].groups[cur_group].faces.push_back(f);
                } else {
                    error("Invalid indices in face definition (line ", cur_line, ").");
                    err_count++;
                }
            }
        } else if (*ptr == 'g' && std::isspace(ptr[1])) {
            file.objects[cur_object].groups.emplace_back();
            cur_group++;
        } else if (*ptr == 'o' && std::isspace(ptr[1])) {
            file.objects.emplace_back();
            cur_object++;

            file.objects[cur_object].groups.emplace_back();
            cur_group = 0;
        } else if (!std::strncmp(ptr, "usemtl", 6) && std::isspace(ptr[6])) {
            ptr += 6;

            ptr = strip_spaces(ptr);
            char* base = ptr;
            ptr = strip_text(ptr);

            const std::string mtl_name(base, ptr);
            cur_mtl = mtl_lib.ids.find(mtl_name)->second;// std::find(mtl_lib.list.begin(), mtl_lib.list.end(), mtl_name) - mtl_lib.list.begin();
           // std::cout<<"cur mtl "<<mtl_name<<" "<<cur_mtl<<"\n";
            if (cur_mtl == (int)mtl_lib.list.size()) {
                warn("Load Unknown material'", mtl_name, " set to dummy ' (line ", cur_line, ").");
                cur_mtl = std::find(mtl_lib.list.begin(), mtl_lib.list.end(), "") - mtl_lib.list.begin();
            }
        } else if (!std::strncmp(ptr, "mtllib", 6) && std::isspace(ptr[6])) {
            ptr += 6;

            ptr = strip_spaces(ptr);
            char* base = ptr;
            ptr = strip_text(ptr);

            const std::string lib_name(base, ptr);

            file.mtl_libs.push_back(lib_name);
        } else if (*ptr == 's' && std::isspace(ptr[1])) {
            // Ignore smooth commands
        } else {
            error("Unknown command '", ptr, "' (line ", cur_line, ").");
            err_count++;
        }
    }
    file.bbox.extend(0.1); 
    
    return (err_count == 0);
}

static bool parse_mtl(std::istream& stream, obj::MaterialLib& mtl_lib) {
    const int max_line = 1024;
    int err_count = 0, cur_line = 0;
    char line[max_line];

    std::string mtl_name;
    auto current_material = [&] () -> obj::Material& {
        return mtl_lib.map[mtl_name];
    };

    while (stream.getline(line, max_line)) {
        cur_line++;

        // Strip spaces
        char* ptr = strip_spaces(line);

        // Skip comments and empty lines
        if (*ptr == '\0' || *ptr == '#')
            continue;

        remove_eol(ptr);

        if (!std::strncmp(ptr, "newmtl", 6) && std::isspace(ptr[6])) {
            ptr = strip_spaces(ptr + 7);
            char* base = ptr;
            ptr = strip_text(ptr);

            mtl_name = std::string(base, ptr);
            if (mtl_lib.map.find(mtl_name) != mtl_lib.map.end()) {
                warn("Material redefinition for '", mtl_name, "' (line ", cur_line, ").");
                continue;
            }
        } else if (ptr[0] == 'K') {
            if (ptr[1] == 'a' && std::isspace(ptr[2])) {
                auto& mat = current_material();
                mat.ka[0] = std::strtof(ptr + 3, &ptr);
                mat.ka[1] = std::strtof(ptr, &ptr);
                mat.ka[2] = std::strtof(ptr, &ptr);
            } else if (ptr[1] == 'd' && std::isspace(ptr[2])) {
                auto& mat = current_material();
                mat.kd[0] = std::strtof(ptr + 3, &ptr);
                mat.kd[1] = std::strtof(ptr, &ptr);
                mat.kd[2] = std::strtof(ptr, &ptr);
            } else if (ptr[1] == 's' && std::isspace(ptr[2])) {
                auto& mat = current_material();
                mat.ks[0] = std::strtof(ptr + 3, &ptr);
                mat.ks[1] = std::strtof(ptr, &ptr);
                mat.ks[2] = std::strtof(ptr, &ptr);
            } else if (ptr[1] == 'e' && std::isspace(ptr[2])) {
                auto& mat = current_material();
                mat.ke[0] = std::strtof(ptr + 3, &ptr);
                mat.ke[1] = std::strtof(ptr, &ptr);
                mat.ke[2] = std::strtof(ptr, &ptr);
            } else {
                error("Invalid command '", ptr, "' (line ", cur_line , ").");
                err_count++;
            }
        } else if (ptr[0] == 'N') {
            if (ptr[1] == 's' && std::isspace(ptr[2])) {
                auto& mat = current_material();
                mat.ns = std::strtof(ptr + 3, &ptr);
            } else if (ptr[1] == 'i' && std::isspace(ptr[2])) {
                auto& mat = current_material();
                mat.ni = std::strtof(ptr + 3, &ptr);
            } else {
                error("Invalid command '", ptr, "' (line ", cur_line , ").");
                err_count++;
            }
        } else if (ptr[0] == 'T') {
            if (ptr[1] == 'f' && std::isspace(ptr[2])) {
                auto& mat = current_material();
                mat.tf.x = std::strtof(ptr + 3, &ptr);
                mat.tf.y = std::strtof(ptr, &ptr);
                mat.tf.z = std::strtof(ptr, &ptr);
            } else if (ptr[1] == 'r' && std::isspace(ptr[2])) {
                auto& mat = current_material();
                mat.tr = std::strtof(ptr + 3, &ptr);
            } else {
                error("Invalid command '", ptr, "' (line ", cur_line , ").");
                err_count++;
            }
        } else if (ptr[0] == 'd' && std::isspace(ptr[1])) {
            auto& mat = current_material();
            mat.d = std::strtof(ptr + 2, &ptr);
        } else if (!std::strncmp(ptr, "illum", 5) && std::isspace(ptr[5])) {
            auto& mat = current_material();
            mat.illum = std::strtof(ptr + 6, &ptr);
        } else if (!std::strncmp(ptr, "map_Ka", 6) && std::isspace(ptr[6])) {
            auto& mat = current_material();
            mat.map_ka = std::string(strip_spaces(ptr + 7));
        } else if (!std::strncmp(ptr, "map_Kd", 6) && std::isspace(ptr[6])) {
            auto& mat = current_material();
            mat.map_kd = std::string(strip_spaces(ptr + 7));
        } else if (!std::strncmp(ptr, "map_Ks", 6) && std::isspace(ptr[6])) {
            auto& mat = current_material();
            mat.map_ks = std::string(strip_spaces(ptr + 7));
        } else if (!std::strncmp(ptr, "map_Ke", 6) && std::isspace(ptr[6])) {
            auto& mat = current_material();
            mat.map_ke = std::string(strip_spaces(ptr + 7));
        } else if (!std::strncmp(ptr, "map_bump", 8) && std::isspace(ptr[8])) {
            auto& mat = current_material();
            mat.map_bump = std::string(strip_spaces(ptr + 9));
        } else if (!std::strncmp(ptr, "bump", 4) && std::isspace(ptr[4])) {
            auto& mat = current_material();
            mat.map_bump = std::string(strip_spaces(ptr + 5));
        } else if (!std::strncmp(ptr, "map_d", 5) && std::isspace(ptr[5])) {
            auto& mat = current_material();
            mat.map_d = std::string(strip_spaces(ptr + 6));
        } else {
            warn("Unknown command '", ptr, "' (line ", cur_line , ").");
        }
    }


    return (err_count == 0);
}

bool load_obj(const std::string path, obj::File& obj_file, obj::MaterialLib& mtl_lib) {
    // Parse the OBJ file
    std::ifstream stream(path);
    bool res = stream && parse_obj(stream, obj_file, mtl_lib);
    printf("load obj  vertices %d face %d \n", obj_file.vertices.size());
    return res;
}

bool load_mtl(const FilePath& path, obj::MaterialLib& mtl_lib) {
    // Parse the MTL file
    std::ifstream stream(path);
    return stream && parse_mtl(stream, mtl_lib);
}

void write_obj(const TriMesh &tri_mesh, const MaterialLib& mtl_lib , int c) {
	std::ofstream outfile;
    printf("write obj %d \n", c);
	outfile.open("obj_chunk_" + std::to_string(c) + ".obj");	

	int trx_size = tri_mesh.indices.size() / 4;
	int vtx_size = tri_mesh.vertices.size();
	printf("vtx_size %d\n", vtx_size);
	for(int i = 0; i < vtx_size; i++) {
		outfile << "v " <<std::fixed<<std::setprecision(4) 
						<< tri_mesh.vertices[i].x << " " 
						<< tri_mesh.vertices[i].y << " " 
						<< tri_mesh.vertices[i].z<<std::endl; 
	}

//	for(int i = 0; i < vtx_size; i++) {
//		outfile << "vn " <<std::fixed<<std::setprecision(4) 
//						<< tri_mesh->normals[i].x << " " 
//						<< tri_mesh->normals[i].y << " " 
//						<< tri_mesh->normals[i].z<<std::endl; 
//	}
//	for(int i = 0; i < size * 3; i++) {
//		outfile << "vt " <<std::fixed<<std::setprecision(4) 
//						<< tri_mesh->texcoords[i].x << " " 
//						<< tri_mesh->texcoords[i].y << " "
//						<< 0.0 << std::endl; 
//	}
	
	for(int i = 0; i < trx_size; i++){
//		if (has_uv)
//		{
//			fprintf(file, "f %d/%d %d/%d %d/%d\n", triangles[i].v[0]+1, uv, triangles[i].v[1]+1, uv+1, triangles[i].v[2]+1, uv+2);
//			uv += 3;
//		}
//		else
//		{
			outfile << "f " << tri_mesh.indices[i * 4] + 1 <<" "
					        << tri_mesh.indices[i * 4 + 1] + 1 <<" "
							<< tri_mesh.indices[i * 4 + 2] + 1<<std::endl;
            
			outfile << "mtl " <<tri_mesh.indices[i * 4 + 3]<<" "<< mtl_lib.list[tri_mesh.indices[i * 4 + 3]] <<"\n";
//		}
		//fprintf(file, "f %d// %d// %d//\n", triangles[i].v[0]+1, triangles[i].v[1]+1, triangles[i].v[2]+1); //more compact: remove trailing zeros
	}
	outfile.close();
}

void write_light_obj(std::vector<float3> verts) {
	std::ofstream outfile;
	outfile.open("light.obj");	

	int vtx_size = verts.size();
    int tri_size = verts.size() / 3;
	printf("vtx_size %d\n", vtx_size);
	for(int i = 0; i < vtx_size; i++) {
		outfile << "v " <<std::fixed<<std::setprecision(4) 
						<< verts[i].x << " " 
						<< verts[i].y << " " 
						<< verts[i].z<<std::endl; 
	}

//	for(int i = 0; i < vtx_size; i++) {
//		outfile << "vn " <<std::fixed<<std::setprecision(4) 
//						<< tri_mesh->normals[i].x << " " 
//						<< tri_mesh->normals[i].y << " " 
//						<< tri_mesh->normals[i].z<<std::endl; 
//	}
//	for(int i = 0; i < size * 3; i++) {
//		outfile << "vt " <<std::fixed<<std::setprecision(4) 
//						<< tri_mesh->texcoords[i].x << " " 
//						<< tri_mesh->texcoords[i].y << " "
//						<< 0.0 << std::endl; 
//	}
	
	for(int i = 0; i < tri_size; i++){
        outfile << "f " << i * 3 + 1 <<" "
                        << i * 3 + 2 <<" "
                        << i * 3 + 3<<std::endl;
            
		//fprintf(file, "f %d// %d// %d//\n", triangles[i].v[0]+1, triangles[i].v[1]+1, triangles[i].v[2]+1); //more compact: remove trailing zeros
	}
	outfile.close();
}

void compute_face_normals(const std::vector<uint32_t>& indices,
                                 const std::vector<float3>& vertices,
                                 std::vector<float3>& face_normals,
                                 size_t first_index) 
{
    for (auto i = first_index, k = indices.size(); i < k; i += 4) {
        const float3& v0 = vertices[indices[i + 0]];
        const float3& v1 = vertices[indices[i + 1]];
        const float3& v2 = vertices[indices[i + 2]];
        face_normals[i / 4] = normalize(cross(v1 - v0, v2 - v0));
    }
}

void compute_vertex_normals(const std::vector<uint32_t>& indices,
                                   const std::vector<float3>& face_normals,
                                   std::vector<float3>& normals,
                                   size_t first_index) 
{
    for (auto i = first_index, k = indices.size(); i < k; i += 4) {
        float3& n0 = normals[indices[i + 0]];
        float3& n1 = normals[indices[i + 1]];
        float3& n2 = normals[indices[i + 2]];
        const float3& n = face_normals[i / 4];
        n0 += n;
        n1 += n;
        n2 += n;
    }
}

void mesh_add(TriMesh &tri_mesh, TriMesh &sub_mesh) {
    auto vtx_offset = tri_mesh.vertices.size();
    auto idx_offset = tri_mesh.indices.size() / 4;
    auto vtx_increment = sub_mesh.vertices.size();
    auto idx_increment = sub_mesh.indices.size() / 4;
    
    tri_mesh.indices.resize(4 * (idx_offset + idx_increment));
    tri_mesh.vertices.resize(vtx_offset + vtx_increment);
    tri_mesh.texcoords.resize(vtx_offset + vtx_increment);
    tri_mesh.normals.resize(vtx_offset + vtx_increment);
    tri_mesh.face_normals.resize(idx_offset + idx_increment);

    memcpy(&tri_mesh.indices[idx_offset * 4], sub_mesh.indices.data(), 4 * idx_increment * sizeof(uint32_t)); 
    int idx_size = 4 * (idx_offset + idx_increment);
    printf("vtx offset %ld idx_offset%ld idx increment%ld vtx increment %ld %ld\n", vtx_offset, idx_offset, idx_increment, vtx_increment, tri_mesh.indices[idx_offset * 4] + vtx_offset);
    for(int i = idx_offset * 4; i < idx_size; i += 4) {
       // std::cout<<"sub mesh tri "<<i / 4<<" "<<tri_mesh.indices[i]<<" "<<tri_mesh.indices[i + 1]<<" "<<tri_mesh.indices[i + 2]<<" "<<tri_mesh.indices[i + 3] <<"\n";
        tri_mesh.indices[i]     += vtx_offset; 
        tri_mesh.indices[i + 1] += vtx_offset; 
        tri_mesh.indices[i + 2] += vtx_offset;
       // std::cout<<"sub mesh tri "<<i / 4<<" "<<tri_mesh.indices[i]<<" "<<tri_mesh.indices[i + 1]<<" "<<tri_mesh.indices[i + 2]<<" "<<tri_mesh.indices[i + 3] <<"\n";
    }
//    printf("normals %ld, texcoords %ld, face normals %ld \n", sub_mesh.normals.size(), sub_mesh.texcoords.size(), sub_mesh.face_normals.size());
    memcpy(&tri_mesh.vertices[vtx_offset], sub_mesh.vertices.data(), vtx_increment * sizeof(float3)); 
    memcpy(&tri_mesh.normals[vtx_offset], sub_mesh.normals.data(), vtx_increment * sizeof(float3)); 
    memcpy(&tri_mesh.texcoords[vtx_offset], sub_mesh.texcoords.data(), vtx_increment * sizeof(float2)); 
    memcpy(&tri_mesh.face_normals[idx_offset], sub_mesh.face_normals.data(), idx_increment * sizeof(float2)); 
}

void virtual_face(TriMesh &tri_mesh, size_t mtl_size, BBox& bbox, int axis) {
    std::vector<TriIdx> triangles;
    auto vtx_offset = tri_mesh.vertices.size();                                                       
    auto idx_offset = tri_mesh.indices.size();                                                        
    tri_mesh.indices.resize(idx_offset + 4 * 2);                                                      
    tri_mesh.vertices.resize(vtx_offset + 4);                                                         
    tri_mesh.texcoords.resize(vtx_offset + 4);                                                        
    tri_mesh.normals.resize(vtx_offset + 4);                                                          
   
    printf("face %f %f %f %f %f %f\n", bbox.min.x, bbox.min.y, bbox.min.z, bbox.max.x, bbox.max.y, bbox.max.z);

//    for(int i = 0; i < 3; i++) {
//        bbox.min[i] += 0.0001f;
//        bbox.max[i] -= 0.0001f;
//    } 

    //Add 2 triangles
    size_t idx[] = {0 + vtx_offset, 1 + vtx_offset, 2 + vtx_offset, mtl_size, 0 + vtx_offset, 2 + vtx_offset, 3 + vtx_offset, mtl_size};
    for(int i = 0; i < 8; i++) {
        tri_mesh.indices[idx_offset + i] = idx[i];                                                
    }

    if (axis == 0) {
        tri_mesh.vertices[vtx_offset + 0] = float3(bbox.max.x, bbox.max.y, bbox.max.z);                   
        tri_mesh.vertices[vtx_offset + 1] = float3(bbox.max.x, bbox.min.y, bbox.max.z);                   
        tri_mesh.vertices[vtx_offset + 2] = float3(bbox.max.x, bbox.min.y, bbox.min.z);                   
        tri_mesh.vertices[vtx_offset + 3] = float3(bbox.max.x, bbox.max.y, bbox.min.z);                   
    } else if(axis == 1) {
        tri_mesh.vertices[vtx_offset + 0] = float3(bbox.max.x, bbox.max.y, bbox.max.z);                   
        tri_mesh.vertices[vtx_offset + 1] = float3(bbox.min.x, bbox.max.y, bbox.max.z);                   
        tri_mesh.vertices[vtx_offset + 2] = float3(bbox.min.x, bbox.max.y, bbox.min.z);                   
        tri_mesh.vertices[vtx_offset + 3] = float3(bbox.max.x, bbox.max.y, bbox.min.z);                   
    } else {
        tri_mesh.vertices[vtx_offset + 0] = float3(bbox.max.x, bbox.max.y, bbox.max.z);                   
        tri_mesh.vertices[vtx_offset + 1] = float3(bbox.min.x, bbox.max.y, bbox.max.z);                   
        tri_mesh.vertices[vtx_offset + 2] = float3(bbox.min.x, bbox.min.y, bbox.max.z);                   
        tri_mesh.vertices[vtx_offset + 3] = float3(bbox.max.x, bbox.min.y, bbox.max.z);                   
    }
    std::fill(tri_mesh.texcoords.begin() + vtx_offset, tri_mesh.texcoords.end(), float2(0.0f));       
//    printf("before face normal\n");    
    auto face_offset = tri_mesh.face_normals.size();
    tri_mesh.face_normals.resize(face_offset + 2);                      
    float3 face_normal(0.0f);
    face_normal[axis] = bbox.max[axis] == bbox.min[axis]? 1: -1;    

    std::fill(tri_mesh.face_normals.begin() + face_offset, tri_mesh.face_normals.end(), face_normal);           
//    printf("face normal  %f %f %f \n %f %f %f\n ",
//            tri_mesh.face_normals[tri_mesh.face_normals.size() - 1].x, tri_mesh.face_normals[tri_mesh.face_normals.size() - 1].y, tri_mesh.face_normals[tri_mesh.face_normals.size() - 1].z, 
//            tri_mesh.face_normals[tri_mesh.face_normals.size() - 2].x, tri_mesh.face_normals[tri_mesh.face_normals.size() - 2].y, tri_mesh.face_normals[tri_mesh.face_normals.size() - 2].z 
//            );   
    std::fill(tri_mesh.normals.begin() + vtx_offset, tri_mesh.normals.end(), float3(0.0f));           
    compute_vertex_normals(tri_mesh.indices, tri_mesh.face_normals, tri_mesh.normals, idx_offset); 
//    printf("normal  %f %f %f \n %f %f %f\n %f %f %f\n %f %f %f \n",
//            tri_mesh.normals[vtx_offset + 0].x, tri_mesh.normals[vtx_offset + 0].y, tri_mesh.normals[vtx_offset + 0].z, 
//            tri_mesh.normals[vtx_offset + 1].x, tri_mesh.normals[vtx_offset + 1].y, tri_mesh.normals[vtx_offset + 1].z, 
//            tri_mesh.normals[vtx_offset + 2].x, tri_mesh.normals[vtx_offset + 2].y, tri_mesh.normals[vtx_offset + 2].z, 
//            tri_mesh.normals[vtx_offset + 3].x, tri_mesh.normals[vtx_offset + 3].y, tri_mesh.normals[vtx_offset + 3].z 
//            );   
}

void compute_virtual_portal(TriMesh &tri_mesh, size_t mtl_size, BBox& bbox_local, BBox& bbox_global) {
    printf("virtural portal %f %f %f | %f %f %f\n", bbox_local.min[0], bbox_local.min[1], bbox_local.min[2], bbox_local.max[0], bbox_local.max[1], bbox_local.max[2]);
    printf("                %f %f %f | %f %f %f\n", bbox_global.min[0], bbox_global.min[1], bbox_global.min[2], bbox_global.max[0], bbox_global.max[1], bbox_global.max[2]);
    for(int i = 0; i < 3; i++) {
        if(bbox_local.min[i] > bbox_global.min[i]) {
            printf("min %d \n", i);
            BBox tmp = bbox_local;
            tmp.max[i] = tmp.min[i];
            virtual_face(tri_mesh, mtl_size, tmp, i); 
        }
        if(bbox_local.max[i] < bbox_global.max[i]) {
            printf("max %d \n", i);
            BBox tmp = bbox_local;
            virtual_face(tri_mesh, mtl_size, tmp, i); 
        }
    }
}

TriMesh compute_tri_mesh(const File& obj_file, const MaterialLib& mtl_lib, size_t mtl_offset, BBox& bbox, bool virtual_portal) {
    TriMesh tri_mesh;
    for (auto& obj: obj_file.objects) {
        // Convert the faces to triangles & build the new list of indices
        std::vector<TriIdx> triangles;
        std::unordered_map<obj::Index, size_t, HashIndex, CompareIndex> mapping;

        bool has_normals = false;
        bool has_texcoords = false;
        for (auto& group : obj.groups) {
            for (auto& face : group.faces) {
                for (size_t i = 0; i < face.indices.size(); i++) {
                    auto map = mapping.find(face.indices[i]);
                    if (map == mapping.end()) {
                        has_normals |= (face.indices[i].n != 0);
                        has_texcoords |= (face.indices[i].t != 0);
                        mapping.insert(std::make_pair(face.indices[i], mapping.size()));
                    }
                }

                auto v0        = mapping[face.indices[0]];
                auto prev      = mapping[face.indices[1]];
                bool v0_inside = bbox.is_inside(obj_file.vertices[face.indices[0].v]);
                bool vp_inside = bbox.is_inside(obj_file.vertices[face.indices[1].v]);

                for (size_t i = 1; i < face.indices.size() - 1; i++) {
                    auto next = mapping[face.indices[i + 1]];
                    bool vn_inside = bbox.is_inside(obj_file.vertices[face.indices[i + 1].v]);
                    if(v0_inside || vp_inside || vn_inside) 
                        triangles.emplace_back(v0, prev, next, face.material + mtl_offset);
                    prev = next;
                    vp_inside = vn_inside;
                }
            }
        }
        if (triangles.size() == 0) continue;

        // Add this object to the mesh
        auto vtx_offset = tri_mesh.vertices.size();
        auto idx_offset = tri_mesh.indices.size();
        tri_mesh.indices.resize(idx_offset + 4 * triangles.size());
        tri_mesh.vertices.resize(vtx_offset + mapping.size());
        tri_mesh.texcoords.resize(vtx_offset + mapping.size());
        tri_mesh.normals.resize(vtx_offset + mapping.size());

//        printf("triangles size %ld idx offset %ld\n", triangles.size(), idx_offset);
        for (size_t i = 0, n = triangles.size(); i < n; i++) {
            auto& t = triangles[i];
            tri_mesh.indices[idx_offset + i * 4 + 0] = t.v0 + vtx_offset;
            tri_mesh.indices[idx_offset + i * 4 + 1] = t.v1 + vtx_offset;
            tri_mesh.indices[idx_offset + i * 4 + 2] = t.v2 + vtx_offset;
            tri_mesh.indices[idx_offset + i * 4 + 3] = t.m;
//            printf("tri mesh %d %d %d %d \n", t.v0 + vtx_offset, t.v1 + vtx_offset, t.v2 + vtx_offset, t.m);
        }
//        printf("tri mesh size %ld over \n", tri_mesh.indices.size());

        for (auto& p : mapping) {
            tri_mesh.vertices[vtx_offset + p.second] = obj_file.vertices[p.first.v];
        }

        if (has_texcoords) {
            for (auto& p : mapping) {
                tri_mesh.texcoords[vtx_offset + p.second] = obj_file.texcoords[p.first.t];
            }
        } else {
            warn("No texture coordinates are present, using default value.");
            std::fill(tri_mesh.texcoords.begin() + vtx_offset, tri_mesh.texcoords.end(), float2(0.0f));
        }

        // Compute the geometric normals for this mesh
        tri_mesh.face_normals.resize(tri_mesh.face_normals.size() + triangles.size());
        compute_face_normals(tri_mesh.indices, tri_mesh.vertices, tri_mesh.face_normals, idx_offset);

        if (has_normals) {
            // Set up mesh normals
            for (auto& p : mapping) {
                const auto& n = obj_file.normals[p.first.n];
                tri_mesh.normals[vtx_offset + p.second] = n;
            }
        } else {
            // Recompute normals
            warn("No normals are present, recomputing smooth normals from geometry.");
            std::fill(tri_mesh.normals.begin() + vtx_offset, tri_mesh.normals.end(), float3(0.0f));
            compute_vertex_normals(tri_mesh.indices, tri_mesh.face_normals, tri_mesh.normals, idx_offset);
        }
    }
    if(virtual_portal) {
        BBox global = obj_file.bbox;
        compute_virtual_portal(tri_mesh, mtl_lib.list.size(), bbox, global);
    }

    // Re-normalize all the values in the OBJ file to handle invalid meshes
    bool fixed_normals = false;
    for (auto& n : tri_mesh.normals) {
        auto len2 = lensqr(n);
        if (len2 <= std::numeric_limits<float>::epsilon() || std::isnan(len2)) {
            fixed_normals = true;
            n = float3(0.0f, 1.0f, 0.0f);
        } else
            n = n * (1.0f / std::sqrt(len2));
    }

    if (fixed_normals)
        warn("Some normals were incorrect and thus had to be replaced with arbitrary values.");

    return tri_mesh;
}

void read_obj_paths(std::string path, std::vector<std::string> &obj_files) {
    std::string suffix = path.substr(path.size() - 3, path.size());
    if(suffix == "obj") {
        obj_files.emplace_back(path.substr(0, path.size() - 4));
        std::cout<<" "<<path.substr(0, path.size() - 4)<<"\n";
    } else {
        std::ifstream stream;
        stream.open(path, std::ios::in);
     
        std::vector<std::string> tmp;
        char line[4096];
        while (stream.getline(line, 4096)) {
            char* ptr = strip_spaces(line);
            printf("read obj path get line %s\n", ptr);
            if(ptr[0] == '<') {
                ptr = strip_spaces(ptr);
                int ed = find(ptr, '>');
                if(ptr[1] == '/') {
                    std::string s(ptr + 2, ptr + ed);
                    if(s == "scene") {
                        break;
                    } else if (s == "mesh") {
                        obj_files.push_back(tmp.back());
                        std::cout<<s<<" "<<tmp.back()<<"\n";
                    } else {
                        std::cout<<"error tag "<<s<<"\n";
                    }
                    tmp.pop_back();
                    tmp.pop_back();
                } else {
                    std::string s(ptr + 1, ptr + ed);
                    std::cout<<"s "<<s<<"\n";
                    tmp.emplace_back(s);
                }
            } else {
                ptr = strip_spaces(ptr);
                int ed = find(ptr, '<');
                std::string s(ptr, ptr + ed);
                std::cout<<"s "<<s<<"\n";
                tmp.emplace_back(s);
            }
        }
        std::cout<<" "<<obj_files[0]<<"\n";
    }
}

} // namespace obj
