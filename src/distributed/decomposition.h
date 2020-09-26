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
    
    int axit = 1;
    int cur_n = n;
    //choose longest axit splat
    while(cur_n != 1){
        grid[axit] *= 2;
        axit = ++axit % d;
        cur_n /= 2;
    }
    printf("grid n %ld :", n);
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
        printf("mesh div min %f %f %f max %f %f %f\n", bbox.min.x, bbox.min.y, bbox.min.z, bbox.max.x, bbox.max.y, bbox.max.z);
        
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
                    bb.min.x = x - 0.001f;
                    bb.min.y = y - 0.001f;
                    bb.min.z = z - 0.001f;
                    bb.max.x = std::min(bbox.max.x, x + step.x + 0.001f);
                    bb.max.y = std::min(bbox.max.y, y + step.y + 0.001f);
                    bb.max.z = std::min(bbox.max.z, z + step.z + 0.001f);
                    printf("min %f %f %f max %f %f %f\n", bb.min.x, bb.min.y, bb.min.z, bb.max.x, bb.max.y, bb.max.z);
                    list.push_back(bb);  //emplace_back
                }
            }
        }
        return true;
    }
};

struct ImageBlock {
    union {
        struct { int xmin, ymin, xmax, ymax; };
        int values[4];
    };

    ImageBlock() {xmin = 0; ymin = 0; xmax = 0; ymax = 0; } 
    ImageBlock(int xmax, int ymax) : xmax(xmax), ymax(ymax), xmin(0), ymin(0) 
    {}
    
    ImageBlock(int xmin, int ymin, int xmax, int ymax) : xmin(xmin), ymin(ymin), xmax(xmax), ymax(ymax)  
    {}
    
    int operator [] (size_t i) const { return values[i]; }
    int& operator [] (size_t i) { return values[i]; }
};

// Grid hierarchy for scheduling
struct ImageDecomposition {
    float2 scale;
    
    int width, height, spp, proc_spp;
    
    //multi process load same chunk. 
    ImageBlock render_block, chunk_project_block; 

    float3 eye, dir, right, up; 
    float w, h;
    float w_sin, h_sin;
    std::vector<int>   domain;
    std::vector<float> depth;
    
    ImageDecomposition(){}

    float3 project_point_to_image(float3 p) {
          //  let d = vec3_normalize(math, vec3_sub(p, eye));
          //  make_vec3(vec3_dot(d, right) / w, vec3_dot(d, up) / h, -vec3_dot(d, dir))
        float3 d = normalize(p - eye);
        float3 res(dot(d, right) / w_sin, dot(d, up) / h_sin, dot(dir, (p - eye))); 
        printf("d %f %f %f p %f %f %f p2 %f %f %f\n", d.x, d.y, d.z, p.x, p.y, p.z, res.x, res.y, res.z); 
        return res ;
    }
  
    void write_project_result() {
        ImageRgba32 img;
        img.width = width;
        img.height = height;
        img.pixels.reset(new uint8_t[width * height * 4]);
        float3 color[] = { float3(0, 60, 0), float3(120, 0, 0), float3(180, 0, 0),
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

    ImageBlock project_cube_to_image(BBox box, int chunk, bool Record) {
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
                    printf("project p2 %f %f %f\n", p2.x, p2.y, p2.z);
                    float depth = p2.z;
                    p_min = min(p_min, p2);
                    p_max = max(p_max, p2);
                    //z 
                    p_min.z = depth < p_min.z ? depth : p_min.z;
                }
            }
        }
        p_min = max(float3(-1), p_min);
        p_max = min(float3(1), p_max);
        printf("project min %f %f %f max %f %f %f\n", p_min.x, p_min.y, p_min.z, p_max.x, p_max.y, p_max.z); 
        //左上 0 0 
        ImageBlock block( (p_min.x + 1) * width  / 2
                        , (1 - p_max.y) * height / 2
                        , (p_max.x + 1) * width  / 2
                        , (1 - p_min.y) * height / 2);
       
        if(Record) { 
            for(int y = block.ymin; y < block.ymax; y++) {
                for(int x = block.xmin; x < block.xmax; x++) {
                    int id = y * width + x;
                    if(domain[id] == -1 || p_min.z < depth[id]) {
                        domain[id] = chunk;
                        depth[id] = p_min.z;
                    }
                }
            }
        }

        return block;
    }

    int get_most_unloaded_chunk(int *block, int *chunk_map, int chunk_size) {
        int* chunks = new int [chunk_size];
        for(int i = 0; i< chunk_size; i++)
            chunks[i] = 0;
        for(int y = block[1]; y < block[3]; y++) {
            for(int x = block[0]; x < block[2]; x ++) {
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
      //  printf("get most unloaded %d %d %d %d chunk %d\n", block[0], block[1], block[2], block[3], c);
        c = c==-1 ? unloaded : c;
        return c;
    }
    
    void image_domain_decomposition(ImageBlock image, int* chunk_map, int proc_rank, int proc_size) {

        MeshChunk chunks; //MeshChunk get bbox automaticly
        for(int i = 0; i < width * height; i++)
            domain[i] = -1;
        int chunk_size = chunks.size; 
        
        splat(proc_size, &scale[0], 2);
        printf("image scale %f %f\n", scale[0], scale[1]);
        
        for(int i = 0; i < chunk_size; i++) 
            chunk_map[i] = -1;

        for(int i = 0; i < chunk_size; i++) {
            project_cube_to_image(chunks.list[i], i, true);
        }
       
        //global available block 
        ImageBlock gblock = image;
        //ImageBlock gblock = project_cube_to_image(chunks.bbox, 0, false);
        printf("global avail block %d %d %d %d\n", gblock[0], gblock[1], gblock[2], gblock[3]);
        
       // int step_width = width / scale[0];
       // int step_height = height / scale[1];   
        int step_width = (gblock.xmax - gblock.xmin) / scale[0];
        int step_height = (gblock.ymax - gblock.ymin) / scale[1];   
        printf("step width %d height %d\n", step_width, step_height); 

        for(int i = 0; i < proc_size; i++) {
            int x = i % int(scale[0]);
            int y = i / int(scale[0]);

            int xmin = gblock.xmin + x * step_width;
            int ymin = gblock.ymin + y * step_height;
            ImageBlock block(xmin, ymin, std::min(gblock.xmax, xmin + step_width), std::min(gblock.ymax, ymin + step_height));  
            
            int chunk = chunk_size == proc_size ? i : get_most_unloaded_chunk(&block[0], chunk_map, chunks.size); 
           
            printf("rank %d x %d y %d block %d %d %d %d\n", i, x, y, block[0], block[1], block[2], block[3]);
            if(chunk != -1); 
                chunk_map[chunk] = i;

            if(i == proc_rank) {
                render_block = block; 
            }
        }
        printf("renderblock %d %d %d %d\n", render_block[0], render_block[1], render_block[2], render_block[3]);
        write_project_result(); 
    }
    
    void decomposition(int* chunk_map, bool imageDecompose, int rank, int size) {
        if(imageDecompose) {
            image_domain_decomposition(ImageBlock(width, height), chunk_map, 0 /*rank*/, 1/*size*/); 
        } else {
            MeshChunk chunks;
            spp = spp / size;
            for(int i = 0; i < chunks.size; i++) 
                chunk_map[i] = i;
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
        
        w_sin = std::sin(fov * PI / 360.0f);
        h_sin = w_sin / ratio;
        printf("fov w %f h %f w sin %f h sin %f\n", w, h, w_sin, h_sin);
        
        domain.resize(width * height);
        std::fill(domain.begin(), domain.end(), -1);
        depth.resize(width * height);    
    }
    
    int* get_render_block() { return &render_block[0]; }

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



