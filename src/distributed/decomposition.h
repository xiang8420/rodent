#pragma once
#include <vector>
#include <memory>
#include <mutex>
#include "../driver/float3.h"
#include "../driver/bbox.h"
#include "../driver/obj.h"
#include "../driver/image.h"
#include "../driver/interface.h"

#ifndef PI 
#define PI 3.14159265359f
#endif


inline void splat(size_t n, float* grid, int d) {
    for(int i = 0; i < d; i++) 
        grid[i] = 1; 
    
    int axit = d - 1;
    int cur_n = n;
    //choose longest axit splat
    while(cur_n != 1){
        grid[axit] *= 2;
        axit = ++axit % d;
        cur_n /= 2;
    }
    printf("grid n %d :", n);
    for(int i = 0; i < d; i++)
        printf("%f ", grid[i]);
    printf("\n");
}

//Splitting bounding box for distributed computing

struct MeshChunk{
    std::vector<BBox>      list;
    BBox                   bbox;
    float3                 scale;
    float3                 step; 
    int                    size;

    MeshChunk(BBox bbox, int grid_num):bbox(bbox) {
        splat(grid_num, &scale[0], 3);
        printf("mesh chunk scale %f %f %f\n", scale[0], scale[1], scale[2]);
        chunk_division();
        size = scale[0] * scale[1] * scale[2];
    }
    
    MeshChunk() {
        bbox  = BBox(get_bbox());
        scale = float3(get_grid());
        chunk_division();
        size = scale[0] * scale[1] * scale[2];
    }
    
    bool chunk_division() {
      //  printf("min %f %f %f max %f %f %f\n", bbox.min.x, bbox.min.y, bbox.min.z, bbox.max.x, bbox.max.y, bbox.max.z);
        
        // shortest axis and cut it
        float3 length = {bbox.max.x - bbox.min.x, bbox.max.y - bbox.min.y, bbox.max.z - bbox.min.z};
        int axis = 0;
        float minlength = length[0];
        for(int i = 1; i< 3; i++){
           if(minlength > length[i]){
               minlength = length[i];
               axis = i;
           }
        }
        for(int i = 0; i < 3; i++){
            step[i] = length[i] / scale[i];
        }

        for(int i = 0; i < scale[0]; i++){
            float x = bbox.min.x + i * step.x;
            for(int j = 0; j < scale[1]; j++){
                float y = bbox.min.y + j * step.y;
                for(int k = 0; k < scale[2]; k++){
                    float z = bbox.min.z + k * step.z;
                    struct BBox bb;
                    bb.min.x = x - 0.1f;
                    bb.min.y = y - 0.1f;
                    bb.min.z = z - 0.1f;
                    bb.max.x = ((i == scale[0] - 1)? bbox.max.x:x + step.x) + 0.1f;
                    bb.max.y = ((j == scale[1] - 1)? bbox.max.y:y + step.y) + 0.1f;
                    bb.max.z = ((k == scale[2] - 1)? bbox.max.z:z + step.z) + 0.1f;
                  //  printf("min %f %f %f max %f %f %f\n", bb.min.x, bb.min.y, bb.min.z, bb.max.x, bb.max.y, bb.max.z);
                    list.push_back(bb);  //emplace_back
                }
            }
        }
        return true;
    }
};

// Grid hierarchy for scheduling
struct ImageDecomposition {
    float2 scale;
    
    int width, height, spp, proc_spp;
    int region[4]; 
    float3 eye, dir, right, up; 
    float w, h;
    std::vector<int>   domain;
    std::vector<float> depth;
    
    ImageDecomposition(){}

    float3 project_point_to_image(float3 p) {
        float3 d = normalize(p - eye); 
        float3 res(dot(d, right) / w, dot(d, up) / h, length(p - eye)/*-dot(d, dir)*/); 
       // printf("d %f %f %f p %f %f %f p2 %f %f %f\n", d.x, d.y, d.z, p.x, p.y, p.z, res.x, res.y, res.z); 
        return res ;
    }
  
    void write_project_result() {
        ImageRgba32 img;
        img.width = width;
        img.height = height;
        img.pixels.reset(new uint8_t[width * height * 4]);
        float3 color[] = { float3(60, 0, 0), float3(120, 0, 0), float3(180, 0, 0),
                           float3(0, 60, 0), float3(0, 120, 0), float3(0, 180, 0),
                           float3(0, 0, 60), float3(0, 0, 120), float3(0, 0, 180)};

        for (size_t y = 0; y < height; ++y) {
            for (size_t x = 0; x < width; ++x) {
                if(domain[y * width + x] != -1) {
                    int c = domain[y * width + x];
                    img.pixels[4 * (y * width + x) + 0] = color[c].x;
                    img.pixels[4 * (y * width + x) + 1] = color[c].y;
                    img.pixels[4 * (y * width + x) + 2] = color[c].z;
                    img.pixels[4 * (y * width + x) + 3] = 255;
                }
            }
        }
        std::string out = "chunk.png";
        if (!save_png(out, img))
            error("Failed to save PNG file\n");
    }

    void project_cube_to_image(BBox box, int chunk) {
        float3 p_min(1), p_max(-1);
        p_min.z = 0xff;
        float3 p;
        for(int x = 0; x < 2; x++) {
            p.x = x == 0 ? box.min.x : box.max.x;
            for(int y = 0; y < 2; y++) {
                p.y = y == 0 ? box.min.y : box.max.y;
                for(int z = 0; z < 2; z++) {
                    
                    p.z = z == 0 ? box.min.z : box.max.z;
                    float3 p2 = project_point_to_image(p); 
                    float depth = p2.z;
                    p_min = max(float3(-1), min(p_min, p2));
                    p_max = min(float3(1),  max(p_max, p2));
                    //z 
                    p_min.z = depth < p_min.z ? depth : p_min.z;
                }
            }
        }
       // printf("project min %f %f %f max %f %f %f\n", p_min.x, p_min.y, p_min.z, p_max.x, p_max.y, p_max.z); 
        //左上 0 0 
        int pid[4]; 
        pid[0] = (p_min.x + 1) * width / 2;
        pid[1] = (p_max.x + 1) * width / 2;
        pid[2] = (1 - p_max.y) * height / 2;
        pid[3] = (1 - p_min.y) * height / 2;

       
        for(int y = pid[2]; y < pid[3]; y++) {
            for(int x = pid[0]; x < pid[1]; x++) {
                int id = y * width + x;
                if(domain[id] == -1 || p_min.z < depth[id]) {
                    domain[id] = chunk;
                    depth[id] = p_min.z;
                }
            }
        } 
    }

    int get_most_unloaded_chunk(int *region, int *chunk_map, int chunk_size) {
        int* chunks = new int [chunk_size];
        for(int i = 0; i< chunk_size; i++)
            chunks[i] = 0;
        for(int y = region[1]; y < region[3]; y++) {
            for(int x = region[0]; x < region[2]; x ++) {
                if(domain[y * width + x] >= 0)
                    chunks[domain[y * width + x]] ++;
            }
        }
        int max = 0, c = -1;
        int unloaded = -1;
        for(int i = 0; i < chunk_size; i++) {
            if(chunks[i] > max && chunk_map[i] == -1) {
                max = chunks[i];
                c = i;
            } else if(chunk_map[i] == -1){
                unloaded = i;
            }
        }
      //  printf("get most unloaded %d %d %d %d chunk %d\n", region[0], region[1], region[2], region[3], c);
        c = c==-1 ? unloaded : c;
        return c;
    }
    
    void image_domain_decomposition(int* render_region, int* chunk_map, int proc_rank, int proc_size) {

        MeshChunk chunks;
        for(int i = 0; i < width * height; i++)
            domain[i] = -1;
        int chunk_size = chunks.size; 
        
        splat(proc_size, &scale[0], 2);
        printf("image scale %f %f\n", scale[0], scale[1]);
        
        for(int i = 0; i < chunk_size; i++) 
            chunk_map[i] = -1;

        for(int i = 0; i < chunk_size; i++) {
            project_cube_to_image(chunks.list[i], i);
        }
        
        int step_width = width / scale[0];
        int step_height = height / scale[1];   
        printf("step width %d height %d\n", step_width, step_height); 

        for(int i = 0; i < proc_size; i++) {
            int region[4];
            int x = i % int(scale[0]);
            int y = i / int(scale[0]);

            region[0] = x * step_width; 
            region[1] = y * step_height;  
            region[2] = std::min(width, region[0] + step_width); 
            region[3] = std::min(width, region[1] + step_height); 
            
            int chunk = chunk_size == proc_size ? i : get_most_unloaded_chunk(&region[0], chunk_map, chunks.size); 
           
            if(chunk != -1); 
                chunk_map[chunk] = i;

            if(i == proc_rank) {
                for(int k = 0; k < 4; k++)
                    render_region[k] = region[k]; 
            }
        }
        printf("renderregion %d %d %d %d\n", render_region[0], render_region[1], render_region[2], render_region[3]);
        write_project_result(); 
    }
    
    void decomposition(int* chunk_map, bool imageDecompose, int rank, int size) {
        region[0] = 0;     region[1] = 0;
        region[2] = width; region[3] = height;
        if(imageDecompose) {
            image_domain_decomposition(region, chunk_map, rank, size); 
        } else {
            spp = spp / size;
        }
    }
   


    ImageDecomposition(float3 e, float3 d, float3 u, float fov, int width, int height, int spp)
            : width(width), height(height), spp(spp) 
    {
        float ratio = (float) width / (float)height;
        eye = e;
        dir = normalize(d);
        right = normalize(cross(dir, u));
        up = normalize(cross(right, dir));

        w = std::tan(fov * PI / 360.0f);
        h = w / ratio;
        
        domain.resize(width * height);
        std::fill(domain.begin(), domain.end(), -1);
        depth.resize(width * height);    
    }
    
    int* get_region() { return &region[0]; }

    int get_spp() { return spp; }

    void rotate(float yaw, float pitch) {
        dir = ::rotate(dir, right,  -pitch);
        dir = ::rotate(dir, up,     -yaw);
        dir = normalize(dir);
        right = normalize(cross(dir, up));
        up = normalize(cross(right, dir));
    }

    void move(float x, float y, float z) {
        eye += right * x + up * y + dir * z;
    }
    
};



