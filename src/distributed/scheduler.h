#pragma once
#include <vector>
#include <memory>
#include <mutex>
#include <utility>
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

    int axit = 0;
    int cur_n = n;
    //choose longest axit splat
    while(cur_n != 1){
        grid[axit] *= 2;
        axit = ++axit % d;
        cur_n /= 2;
    }
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
        chunk_division();
        size = scale[0] * scale[1] * scale[2];
    }
    
    MeshChunk() {
        bbox  = BBox(get_bbox());
        scale = float3(get_chunk());
        chunk_division();
        size = get_chunk_num();
    }
    
    bool chunk_division() {
        
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
                    list.push_back(bb);  //emplace_back
                }
            }
        }
        return true;
    }
};
  
static void write_project_result(int width, int height, int* domain) {
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
//    if (!save_png(out, img))
//        error("Failed to save PNG file\n");
}

static void save_image_its(int* reduce_buffer, int chunk_size, float spp) {
    int res = LIGHT_FIELD_RES;
    ImageRgba32 img;
    img.width = res * 6 + 20;
    img.height = res * chunk_size;
    img.pixels.reset(new uint8_t[img.width * img.height * 4]);

    int width = img.width; 
    int height = img.height;
    
    std::vector<rgb> color;
    for(int i = 0; i < chunk_size; i++) {
        rgb col(rand() % 255, rand() % 255, rand() % 255 ); 
        color.emplace_back(col);
        for(int h = 0; h < 128; h++) {
            int ph = i * 128 + h; 
            for(int w = 0; w < 19; w ++) {

                img.pixels[4 * (ph * width + w) + 0] = col.x; 
                img.pixels[4 * (ph * width + w) + 1] = col.y;
                img.pixels[4 * (ph * width + w) + 2] = col.z;
                img.pixels[4 * (ph * width + w) + 3] = 255;
            }
        }
    }

    rgb black(0.0, 0.0, 0.0);

    int all_face_size = res * res * 6;
    int face_size = res * res;
    int inv_spp = 1;
    for(int i = 0; i < chunk_size; i++) {
        int h_st = res * i; 
        for(int j = 0; j < 6; j++) {
            int w_st = 20 + res * j;
            for(int u = 0; u < res; u++) {
                int pu = w_st + u;
                for(int v = 0; v < res; v++) {
                    int pv = h_st + v;
                    if(v == res - 1 || u == res - 1) {
                        img.pixels[4 * (pv * width + pu) + 0] = 0;
                        img.pixels[4 * (pv * width + pu) + 1] = 255;
                        img.pixels[4 * (pv * width + pu) + 2] = 0;
                        img.pixels[4 * (pv * width + pu) + 3] = 255;
                    } else { 
                        int id = all_face_size * i + face_size * j + (res - 1 - v) * res + u;
                        int recv_its = reduce_buffer[id];
                    //    int its = (recv_its & 0xFF) - 1;
                        int its = ((recv_its >> 8) & 0xFF ) - 1;
                    //    int its = ((recv_its >> 16) & 0xFF ) - 1;
                    //    int its = ((recv_its >> 24) & 0xFF ) - 1;
                
                        rgb *col;
                        if(its >= 254 || its == -1) col = &black;
                        else { 
                            if(its < 0|| its > chunk_size - 1) {
                             //   printf("error its %d\n", its);
                                col = &black;
                            } else
                                col = &color[its];
                        }
                        img.pixels[4 * (pv * width + pu) + 0] = col->x;
                        img.pixels[4 * (pv * width + pu) + 1] = col->y;
                        img.pixels[4 * (pv * width + pu) + 2] = col->z;
                        img.pixels[4 * (pv * width + pu) + 3] = 255;
                    }
                }
            }
        }
    }
    if (!save_png(std::string("picture/light_field_its.png"), img))
        error("Failed to save PNG file light_field.png");
}

static void save_image_ctrb(int* reduce_buffer, int chunk_size, float spp) {
    int res = LIGHT_FIELD_RES;
    ImageRgba32 img;
    img.width = res * 6;
    img.height = res * chunk_size;
    img.pixels.reset(new uint8_t[img.width * img.height * 4]);

    int width = img.width; 
    int height = img.height;
    int all_face_size = res * res * 6;
    int face_size = res * res;
    int inv_spp = 1;//5 / spp;

    for(int i = 0; i < chunk_size; i++) {
        int h_st = res * i; 
        for(int j = 0; j < 6; j++) {
            int w_st = res * j;
            for(int u = 0; u < res; u++) {
                int pu = w_st + u;
                for(int v = 0; v < res; v++) {
                    int pv = h_st + v;
                    
                    int lu = reduce_buffer[all_face_size * i + face_size * j + (res - 1 - v) * res + u];
                    int lu_q0 = lu & 0xFF;
                    int lu_q1 = lu >> 8  & 0xFF;
                    int lu_q2 = lu >> 16 & 0xFF;
                    int lu_q3 = lu >> 24 & 0xFF;
                    
                    int lu_q = lu_q0;// + lu_q1 + lu_q2 + lu_q3; 
                    img.pixels[4 * (pv * width + pu) + 0] = lu_q; // (lu_q >> 4) * inv_spp; 
                    img.pixels[4 * (pv * width + pu) + 1] = lu_q; // ((lu_q >> 2) & 0x3) * inv_spp;
                    img.pixels[4 * (pv * width + pu) + 2] = lu_q; // (lu_q & 0x3) * inv_spp;
                    img.pixels[4 * (pv * width + pu) + 3] = 255;
                    if(v == res - 1 || u == res - 1) {
                        img.pixels[4 * (pv * width + pu) + 0] = 255;
                        img.pixels[4 * (pv * width + pu) + 1] = 0;
                        img.pixels[4 * (pv * width + pu) + 2] = 0;
                        img.pixels[4 * (pv * width + pu) + 3] = 255;
                    } 
                } 
            }
        } 
    }
    if (!save_png(std::string("picture/light_field_ctrb.png"), img))
        error("Failed to save PNG file light_field.png");
}
//
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

struct Camera {
    float3 eye, dir, right, up; 
    float w, h;
    float w_sin, h_sin;
    int width, height, spp;
    int iter;

    Camera(float3 e, float3 d, float3 u, float fov, int width, int height) 
        : width(width), height(height)
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
        iter = 0;
    }

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

    float3 project_point_to_image(float3 p) {
        float3 d = normalize(p - eye);
        float3 res(dot(d, right) / w_sin, dot(d, up) / h_sin, dot(dir, (p - eye))); 
        return res ;
    }
};

// Grid hierarchy for scheduling
struct Scheduler {
    float2 scale;
    
    int width, height, spp, proc_spp;
    int block_count; 
    
    //multi process load same chunk. 
    ImageBlock render_block, chunk_project_block; 

    std::vector<int>   domain;
    std::vector<float> depth;
    std::vector<ImageBlock> blockList;

    Camera* camera;
    std::vector<std::pair<int, int>> chunk_map;

    Scheduler(){}

    ImageBlock project_cube_to_image(Camera *camera, BBox box, int chunk, bool Record, ImageBlock image);

    void image_domain_decomposition(Camera *cam, int block, int proc_rank, int proc_size, bool simple_mesh, bool sync);
    
    void split_image_block(Camera *cam, int block, int proc_rank, int proc_size);
    
    int chunk_loaded(int chunk); 

    int get_most_unloaded_chunk(int *block, int chunk_size);
   
    void write_chunk_proc(int* chunk_proc); 

    int* get_render_block() { return &render_block[0]; }
    
    void set_render_block(int i) { render_block = blockList[i]; }

    size_t get_block_count() { return block_count; }
    
    int get_spp() { return spp; }
    
    Scheduler( int width, int height, int spp);
};

int Scheduler::chunk_loaded(int chunk) {
    int i = 0;
    for(auto& iter: chunk_map) {
        if(chunk == iter.first)
            return iter.second;
    }
    return -1;
}

int Scheduler::get_most_unloaded_chunk(int *block, int chunk_size) {
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
        if(chunk_loaded(i) == -1) {
            if(chunks[i] > max ) {
                max = chunks[i];
                c = i;
            } else {
                //这是个随意的unloaded
                unloaded = i;
            }
        }
    }
    c = c==-1 ? unloaded : c;
    return c;
}

void Scheduler::write_chunk_proc(int* chunk_proc) {
    int i = 0;
    for(auto& iter: chunk_map) {
        chunk_proc[i * 2] = iter.first;
        chunk_proc[i * 2 + 1] = iter.second;
        i++;
    }
    
    for(int j = 0; j < i; j++) 
        printf("chunk map %d %d\n", chunk_proc[j * 2], chunk_proc[j * 2 + 1]);

    chunk_proc[i * 2] = -1;
}

void Scheduler::image_domain_decomposition(Camera *cam, int block, int proc_rank, int proc_size, bool simple_mesh, bool sync) {
    camera = cam;
    ImageBlock image(width, height);
    block_count = block;
    MeshChunk chunks;//MeshChunk get bbox automaticly
    int chunk_size = chunks.size; 
    if (block_count == 1) {
        spp = spp / proc_size;
        for(int i = 0; i < chunks.size; i++) {
            chunk_map.push_back(std::make_pair(i, i));
        }
        return;
    } 

    for(int i = 0; i < width * height; i++)
        domain[i] = -1;
    
    splat(proc_size, &scale[0], 2);

    for(int i = 0; i < chunk_size; i++) {
        printf("chunks list %f %f %f max %f %f %f\n", chunks.list[i].min.x, chunks.list[i].min.y, chunks.list[i].min.z
                                                    , chunks.list[i].max.x, chunks.list[i].max.y, chunks.list[i].max.z);
        project_cube_to_image(camera, chunks.list[i], i, true, image);
    }
   
    //global available block 
    ImageBlock gblock = project_cube_to_image(camera, chunks.bbox, 0, false, image);

    // int step_width = width / scale[0];
    // int step_height = height / scale[1];   
    int step_width = (gblock.xmax - gblock.xmin) / scale[0];
    int step_height = (gblock.ymax - gblock.ymin) / scale[1];   

    for(int i = 0; i < proc_size; i++) {
        int x = i % int(scale[0]);
        int y = i / int(scale[0]);

        int xmin = gblock.xmin + x * step_width;
        int ymin = gblock.ymin + y * step_height;
        ImageBlock block(xmin, ymin, std::min(gblock.xmax, xmin + step_width), std::min(gblock.ymax, ymin + step_height));  
        
        if(chunk_size < proc_size)
            error("chunk size ", chunk_size, " proc size ", proc_size, "invalid");

        int chunk;
        if(chunk_size == proc_size) {
            chunk = i;
        } else if(chunk_size > proc_size) {
            chunk = get_most_unloaded_chunk(&block[0], chunks.size); 
        } else if(chunk_size == 1) {
            chunk = 0;
        } else {
            error("chunk size ", chunk_size, " proc size ", proc_size, "invalid");
        }
       
        if(chunk != -1){ 
            chunk_map.push_back(std::make_pair(chunk, i));
            printf("make pair %d %d\n", chunk, i);
            blockList.emplace_back(block);
        }
        if(i == proc_rank) {
            render_block = block; 
        }
    }
    if(proc_size < chunk_size) {
        std::vector<int> unloaded_chunk;
        for(int i = 0; i < chunk_size; i++) {
            bool find = false;
            for(auto& iter: chunk_map) { 
                if(iter.first == i) {
                    find = true;
                    break;
                }
            }
            if(!find) 
                unloaded_chunk.emplace_back(i);
        }
        int proc = 0;
        for(auto& iter :unloaded_chunk) {
            if(sync)
                chunk_map.push_back(std::make_pair(iter, -1));
            else 
                chunk_map.push_back(std::make_pair(iter, proc % proc_size));
            proc++;
        } 
    }
//        printf("renderblock %d %d %d %d\n", render_block[0], render_block[1], render_block[2], render_block[3]);
    //write_project_result(width, height, domain.data()); 
}

ImageBlock Scheduler::project_cube_to_image(Camera *camera, BBox box, int chunk, bool Record, ImageBlock image) {
    if(box.is_inside(camera->eye)) { 
        int num_pixels = width * height;
        for(int i = 0; i < num_pixels; i++) {
            domain[i] = chunk;
            depth[i] = 0;
        }
        return image;
    }

    float3 p_min(1), p_max(-1);
    p_min.z = 0xff;
    float3 p;
    for(int x = 0; x < 2; x++) {
        p.x = x == 0 ? box.min.x : box.max.x;
        for(int y = 0; y < 2; y++) {
            p.y = y == 0 ? box.min.y : box.max.y;
            for(int z = 0; z < 2; z++) {
                
                p.z = z == 0 ? box.min.z : box.max.z;
                float3 p2 = camera->project_point_to_image(p); 
//                    printf("project p2 %f %f %f\n", p2.x, p2.y, p2.z);
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
//        printf("project min %f %f %f max %f %f %f\n", p_min.x, p_min.y, p_min.z, p_max.x, p_max.y, p_max.z); 
    //top-left 0 0 
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

void Scheduler::split_image_block(Camera *cam, int block, int proc_rank, int proc_size) {
    camera = cam;
    block_count = block;
    ImageBlock image(width, height);

    splat(block_count, &scale[0], 2);
   
    //global available block 
    for(int i = 0; i < proc_size; i++) 
        chunk_map.push_back(std::make_pair(0, proc_rank));
    ImageBlock gblock = project_cube_to_image(camera, BBox(get_bbox()), 0, false, image);
    if(block == 1) { 
        render_block = gblock; 
        return;
    }
    
    // int step_width = width / scale[0];
    // int step_height = height / scale[1];   
    int step_width = (gblock.xmax - gblock.xmin) / scale[0];
    int step_height = (gblock.ymax - gblock.ymin) / scale[1];   
//        printf("step width %d height %d\n", step_width, step_height); 

    for(int i = 0; i < block_count; i++) {
        int x = i % int(scale[0]);
        int y = i / int(scale[0]);

        int xmin = gblock.xmin + x * step_width;
        int ymin = gblock.ymin + y * step_height;
        ImageBlock block(xmin, ymin, std::min(gblock.xmax, xmin + step_width), std::min(gblock.ymax, ymin + step_height));  
//            printf("rank %d x %d y %d block %d %d %d %d\n", i, x, y, block[0], block[1], block[2], block[3]);

        blockList.emplace_back(block);

    }
    render_block = blockList[proc_rank]; 
//        printf("renderblock %d %d %d %d\n", render_block[0], render_block[1], render_block[2], render_block[3]);
    //write_project_result(width, height, domain.data()); 
}
    
Scheduler::Scheduler( int width, int height, int spp)
       : width(width), height(height), spp(spp) 
{
    domain.resize(width * height);
    std::fill(domain.begin(), domain.end(), -1);
    depth.resize(width * height);    
}
