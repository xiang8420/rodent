#pragma once
#include <vector>
#include <memory>
#include <cstring>
#include <mutex>
#include <utility>
#include "../driver/float3.h"
#include "../driver/bbox.h"
#include "../driver/obj.h"
#include "../driver/image.h"
#include "../driver/interface.h"
#include "MeshChunk.h"

#ifndef PI 
#define PI 3.14159265359f
#endif

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

static void save_image_its(int* reduce_buffer, int chunk_size) {
    int res = CHUNK_HIT_RES;
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
    if (!save_png(std::string("picture/chunk_hit_its.png"), img))
        error("Failed to save PNG file chunk_hit.png");
}

static void save_image_ctrb(int* reduce_buffer, int chunk_size) {
    int res = CHUNK_HIT_RES;
    ImageRgba32 img;
    img.width = res * 6;
    img.height = res * chunk_size;
    img.pixels.reset(new uint8_t[img.width * img.height * 4]);

    int width = img.width;
    int height = img.height;
    int all_face_size = res * res * 6;
    int face_size = res * res;
    int inv_spp = 1;//5 / spp;

    rgb col1(0, 125, 125); 
    rgb col2(255, 0, 0); 
    rgb col3(0, 0, 0); 

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
                    rgb *col; 
                    if(lu_q == 1) col = &col1;
                    else if(lu_q > 1) col = &col2;
                    else col = &col3;

                    img.pixels[4 * (pv * width + pu) + 0] = col->x;
                    img.pixels[4 * (pv * width + pu) + 1] = col->y;
                    img.pixels[4 * (pv * width + pu) + 2] = col->z;
                    img.pixels[4 * (pv * width + pu) + 3] = 255;
                //    img.pixels[4 * (pv * width + pu) + 0] = lu_q; // (lu_q >> 4) * inv_spp; 
                //    img.pixels[4 * (pv * width + pu) + 1] = lu_q; // ((lu_q >> 2) & 0x3) * inv_spp;
                //    img.pixels[4 * (pv * width + pu) + 2] = lu_q; // (lu_q & 0x3) * inv_spp;
                //    img.pixels[4 * (pv * width + pu) + 3] = 255;
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
    if (!save_png(std::string("picture/chunk_hit_ctrb.png"), img))
        error("Failed to save PNG file chunk_hit.png");
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

struct PassRecord {
    int * data; //(chunk_size + 2) * chunk_size [cur_chk] miss
    int * reduce_data;
    int * cur_chk_send;
    int * cur_chk_recv; //[0] recv 
    int chunk_size;
    
    PassRecord(int chunk_size) : chunk_size(chunk_size) {
        data = new int[(chunk_size + 1) * chunk_size];
        reduce_data = new int[(chunk_size + 1) * chunk_size]; 
    }

    ~PassRecord(){
        delete[] data;
        delete[] reduce_data;
    }

    void set_cur_chk(int chk) {
        cur_chk_send = data + (chunk_size + 1) * chk; 
        cur_chk_recv = cur_chk_send + chunk_size; 
    }

    void write_send(float * rays, int size,  bool primary) {
        int width = primary ? PRIMARY_WIDTH : SECONDARY_WIDTH;
        int *iptr = (int*) rays;
        for(int i = 0; i < size; i++) {
            int org_chk = iptr[i * width + 10] & 0xFF;;
            cur_chk_send[org_chk] ++;
        }
    }

    void write_recv(int n) {
        cur_chk_recv[0] += n;
    }

    void print() {
        printf("pass record : \n");
        int col = chunk_size + 1;
        for(int i = 0; i < chunk_size; i++) {
            for(int j = 0; j < chunk_size + 1; j++) {
                printf(" %d |", data[i * col + j]);
                data[i * col + j] = 0;
            }
            printf("\n");
        }
        printf("\n");
    }
    
    void gather() {
        MPI_Reduce(data, reduce_data, (chunk_size + 1) * chunk_size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
        int* tmp = data;
        data = reduce_data;
        reduce_data = tmp;
    }
};

//(outlist[next].rays_size() 
//+ next.other_send_left / local_chunk_size 
//+ next.recv_map[current] * outlist[current].rays_size() 
////- for(outlist[i] * recv_map[i]) / local_chunk_size
//) * |1 - (outlist[current] / compute_speed - load_time) / max_size| 
// make sure the formula > 0
struct Chunk {
    int id;
    int load_time; //loaded time
    int compute_speed; // run_time = rays_size / compute_speed 
    int other_send_left; //recv from other proc 
    bool simple_chunk;   //min priority
    std::vector<int> recv_map // recv from other local chunk, possiable recv size = chunk_n_size * recv_map[n]

    Chunk(int id, int priority)
        :id(id), priority(priority) 
    {}
};

struct LocalChunks {
    std::vector<Chunk> chunks;
    std::queue<int> history;

    int current, next;
    bool new_loaded;
    
    LocalChunks() {
        new_loaded= true;
        next = -1;
    }

    void insert(int chk, int pri) { chunks.push_back(Chunk(chk, pri)); }
    int size() { return chunks.size(); }
    void set_new_loaded(bool a) { new_loaded = a; }
    bool get_new_loaded() { return new_loaded; }
    Chunk& get_chunk(int i) { return chunks[i]; }
    
    bool find(int chk) {
        for(int i = 0; i < chunks.size(); i++) 
            if(chk == chunks[i].id) return true;
        return false;
    }
    
    void record_history() {
        if(history.size() >= chunks.size()) 
            history.pop();
        history.push(current); 
    }

    void init() {
        for(int i = 0; i < chunks.size(); i ++) {
            Chunk& chk = chunks[i];
            printf("chunk %d priority %d \n ", chk.id, chk.priority);
            if(chk.priority == 0) current = chk.id;
            if(chk.priority == 1) next = -1;//iter.first;
        } 
        record_history(); 
        new_loaded = true;
    }
    
    int select_new_chunk(RayStreamList * outlist) {
        int max = 0, new_chk = -1;
        float weight = 1;
        printf("select new chunk :");
        int local_size = SIMPLE_TRACE ? chunks.size() - 1 : chunks.size();
        for(int i = 0; i < local_size; i ++) {
            Chunk& chk = chunks[i]; 
            int cur_size = outlist[chk.id].size();
            printf("%d %d %d |", chk.id, cur_size, weight * cur_size);
            cur_size = cur_size > 100 ? weight * cur_size : cur_size;
            weight *= 0.8; // 
            if(chk.id != current && cur_size > max) {
                new_chk = chk.id;
                max = outlist[chk.id].size();
            }
        }
        if(SIMPLE_TRACE && new_chk == -1 && outlist[chunks.back().id].size() > 0) 
            new_chk = chunks.back().id;        
        if(next == -1) {
            current = new_chk;
        } else {
            current = next;
            next = new_chk;
        }
        new_loaded = true;
    }
    bool new_chunk() { return next != current && next != -1; }
};

struct ChunkProcs {
    std::vector<int> proc_list;
    int iter;
    
    ChunkProcs () { iter = 0; }
    int get_proc() { return proc_list[iter++ % proc_list.size()]; }
    int set_proc(int proc) { proc_list.emplace_back(proc); }
    int size() {return proc_list.size(); }
    int back() {return proc_list.back(); }
};

struct ChunkManager {
    ChunkProcs * chunk_list;
    LocalChunks local_chunks;  // chunk id, chunk priority
    int chunk_size;

    ChunkManager(int size) : chunk_size(size) { chunk_list = new ChunkProcs[chunk_size]; }
    ~ChunkManager() { delete[] chunk_list; }

    int switch_current_chunk(RayStreamList * outlist) { local_chunks.select_new_chunk(outlist); }
    
    int get_proc(int c) {
        if(c >= chunk_size) { 
            warn(c, " chunk not loaded");
            return -1;
        }
        printf("get chunk %d  proc %d \n", c, chunk_list[c].get_proc());
        return chunk_list[c].get_proc();
    }

    bool is_local_chunk(int c) { return local_chunks.find(c); }
    
    bool update_chunk(int* schedule, int proc_rank) { error("function update chunk not implement"); }

    bool update_chunk(int candidate_proc, int candidate_chunk, int proc_rank) {
        error("check function update chunk first.");
        chunk_list[candidate_chunk].set_proc(candidate_proc);
        if(candidate_proc == proc_rank) {
            local_chunks.insert(candidate_chunk, 1/*priority*/);
            local_chunks.set_new_loaded(true);
        }
        return local_chunks.get_new_loaded();
    }

    void get_start_chunk() { local_chunks.init(); }
    int get_next_chunk() {return local_chunks.next;}
    int get_current_chunk() {return local_chunks.current;}
    bool new_chunk(){return local_chunks.new_chunk();}
};

// Grid hierarchy for scheduling
class Scheduler {
public:
    ImageBlock project_cube_to_image(Camera *camera, BBox box, int chunk, bool Record, ImageBlock image);

    void get_projected_chunk(int proc_rank, int proc_size, ImageBlock& image); 
    void process_left_chunk(int, bool, bool);
    void image_domain_decomposition(bool simple_mesh, bool sync);
    void split_image_block();
    void chunk_reallocation(int* reduce_buffer); 
    int  chunk_loaded(int chunk); 
    int  get_most_unloaded_chunk(int *block, int chunk_size);
    void generate_chunk_manager(Camera *, int, bool, bool);    
    void get_neighbor(int chk_id); 

    int* get_render_block() { return &render_block[0]; }
    void set_render_block(int i) { render_block = blockList[i]; }
    size_t get_block_count() { return block_count; }
    int get_spp() { return spp; }
    Scheduler( int, int, int, int, int);
    
    ~Scheduler();

    Camera* camera;
    ChunkManager * chunk_manager;
private:    
    int width, height, spp, proc_spp;
    int block_count; 
    float2 scale;
    bool load_chunk_hit;
    int proc_size, proc_rank;

    //multi process load same chunk. 
    ImageBlock render_block, chunk_project_block;
    MeshChunk chunks;

    std::vector<int>   domain;
    std::vector<float> depth;
    std::vector<ImageBlock> blockList;

    std::vector<std::pair<int, int>> chunk_map;
    std::vector<int> send_recv_map;
};

Scheduler::Scheduler( int width, int height, int spp, int prank, int psize)
       : width(width), height(height), spp(spp), proc_rank(prank), proc_size(psize)
{
    domain.resize(width * height);
    std::fill(domain.begin(), domain.end(), -1);
    depth.resize(width * height);
    int chunk_size = chunks.size;
    send_recv_map.resize(chunk_size * chunk_size);
    load_chunk_hit = false;
    chunk_manager = new ChunkManager(SIMPLE_TRACE ? chunk_size + 1 : chunk_size);
}

Scheduler::~Scheduler() {
    //delete chunk_manager;
}

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
            }
        }
    }
    return c;
}

void Scheduler::generate_chunk_manager(Camera *cam, int block, bool simple_trace, bool sync) {    
    camera = cam;
    block_count = block;
    if(block_count == proc_size) 
        image_domain_decomposition(simple_trace, sync);
    else
        split_image_block(); 
    
    int priority = 0;
    printf("generate chunk manager: ");
    for(auto& iter: chunk_map) {
        chunk_manager->chunk_list[iter.first].set_proc(iter.second);
        printf("%d %d |", iter.first, chunk_manager->chunk_list[iter.first].back());
        if(iter.second == proc_rank) {
            chunk_manager->local_chunks.insert(iter.first, priority++/*priority*/);
        }
    }
    chunk_manager->get_start_chunk(); 
}

void Scheduler::get_projected_chunk(int proc_rank, int proc_size, ImageBlock& image) {
    int chunk_size = chunks.size; 
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
    int step_width  = (gblock.xmax - gblock.xmin) / scale[0];
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
      
        //if no chunk project to this image block, push (-1, proc_rank)  
        chunk_map.push_back(std::make_pair(chunk, i));
        blockList.emplace_back(block);
        
        if(i == proc_rank) {
            render_block = block; 
        }
    }
}

struct ChunkHitInfo {
    int id;
    int send_count, send_rank;
    int recv_count, recv_rank;
    int center;
    //int divergence_rank;
};

bool send_greater_mark(const ChunkHitInfo& a, const ChunkHitInfo& b) {
    return a.send_count > b.send_count;
}
 
bool recv_greater_mark(const ChunkHitInfo& a, const ChunkHitInfo& b) {
    return a.recv_count > b.recv_count;
}

bool id_greater_mark(const ChunkHitInfo& a, const ChunkHitInfo& b) {
    return a.id < b.id;
}
//void Scheduler::process_chunk_hit() {
//    
//}

void Scheduler::process_left_chunk(int proc_size, bool sync, bool simple_trace) {
    
    int chunk_size = chunks.size; 
    if(chunk_map.size() >= chunk_size) return;
   
    //load_chunk_hit = false; 
    // put left chunk to proc
    if(!sync && load_chunk_hit) {
 //   if(false) {
        printf(" process left chunk second frame\n");
        //[2 * i] count, [2 * i + 1] divergence degree 
        std::vector<ChunkHitInfo> chunk_hit_info(chunk_size);
        memset(chunk_hit_info.data(), 0, sizeof(ChunkHitInfo) * chunk_size);
        //get recv / send, 
        for(int i = 0; i < chunk_size; i++) {
            int st = i * chunk_size; 
            chunk_hit_info[i].id = i;
            chunk_hit_info[i].center = -1;
            for(int j = 0; j < chunk_size; j ++) {
                if(send_recv_map[st + j] > 0 ) {
                    chunk_hit_info[i].send_count += send_recv_map[st + j];
                    if(j != i) {
                        chunk_hit_info[j].recv_count += send_recv_map[st + j];
                    }
                }
            }
        }
        
        sort(chunk_hit_info.begin(), chunk_hit_info.end(), recv_greater_mark);
        for(int i = 0; i < chunk_size; i++)
            chunk_hit_info[i].recv_rank = i;

        sort(chunk_hit_info.begin(), chunk_hit_info.end(), send_greater_mark);
        for(int i = 0; i < chunk_size; i++)
            chunk_hit_info[i].send_rank = i;
        
        sort(chunk_hit_info.begin(), chunk_hit_info.end(), id_greater_mark);
        for(int i = 0; i < chunk_size; i++) {
            ChunkHitInfo& chk = chunk_hit_info[i];
            printf("chunk_hit_info chunk %d send %d rank %d  recv %d rank %d\n", 
                    chk.id, chk.send_count, chk.send_rank, chk.recv_count, chk.recv_rank);
        }
        //Get cluster center
        std::vector<int> center;
        for(auto& chk_proc: chunk_map) {
            if(chk_proc.second != -1) {
                for(int i = 0; i < chunk_size; i ++) {
                    if(chunk_hit_info[i].id == chk_proc.first) {
                        chunk_hit_info[i].center = chk_proc.first;
                        center.emplace_back(chk_proc.first);
                        break;
                    }
                }
            }
        }
        
        printf("cluster center size %d list: ", center.size());
        for(int i = 0; i < center.size(); i++)
            printf(" %d |", center[i]);
        printf("\n");
        
        // Avoid center empty, push most send chunk to center
        if(center.empty()) {
            center.emplace_back(chunk_hit_info[0].id);
            chunk_hit_info[0].center = 0;
        }

        printf("new chunk allocation :");
        for(auto& chk: chunk_hit_info) {
            printf("%d  %d | ", chk.id, chk.center); 
        }
        while(center.size() < proc_size) {
            int max_dis = 0, cand_chk = -1;
            for(int i = 0; i < chunk_size; i++) {
                // Find the most frequent data transfer with center, to fill centers
                std::vector<int>::iterator it = find(center.begin(), center.end(), i);
                if (it != center.end()) continue; // it is already in centers
                
                int dis = 0;
                for(auto& cent_chk : center) {
                     dis += send_recv_map[i * chunk_size + cent_chk];
                     dis += send_recv_map[cent_chk * chunk_size + i];
                }
                if(dis > max_dis) {
                    max_dis = dis;
                    cand_chk = i;
                }
            }
            center.emplace_back(cand_chk);
            chunk_hit_info[cand_chk].center = cand_chk;
        }
        printf("cluster center : ");
        for(int i = 0; i < proc_size; i++)
            printf(" %d |", center[i]);
        printf("\n");
        printf("new chunk allocation :");
        for(auto& chk: chunk_hit_info) {
            printf("%d  %d | ", chk.id, chk.center); 
        }

        printf("\n");
        std::vector<int> cluster_size(chunk_size);
        int max_cluster_size = chunk_size / proc_size - 1;
        for(int chk = 0; chk < chunk_size; chk ++) {
            //计算到各个聚类的距离 
            if(chunk_hit_info[chk].center == -1) {
                int min_dist = (1<<31)-1; 
                int cand_center = -1;
                for(auto& cent_chk : center) {
                    int dist = 0; 
                    for(int clst_chk = 0; clst_chk < chunk_size; clst_chk ++) {
                        if(chunk_hit_info[chk].center == cent_chk)
                            dist += send_recv_map[chk * chunk_size + clst_chk] + send_recv_map[cent_chk * chunk_size + clst_chk];
                    }
                    
                    if(dist < min_dist && cluster_size[cent_chk] < max_cluster_size) {
                        min_dist = dist;
                        cand_center = cent_chk;
                    }
                }
                printf("chunk hit %d center %d \n", chk, cand_center);
                cluster_size[cand_center] ++;
                chunk_hit_info[chk].center = cand_center;
            }
        }

        printf("\nnew chunk allocation :");
        for(auto& chk: chunk_hit_info) {
            printf("%d  %d | ", chk.id, chk.center); 
        }
        printf("\n");
   
        chunk_map.clear(); 
        for(int i = 0; i < center.size(); i++) {
            int cent_chk = center[i];
            printf("center %d | ", cent_chk);
            for(auto& chk: chunk_hit_info) {
                if(chk.center == cent_chk) {
                    printf("%d ", chk.id);
                    chunk_map.push_back(std::make_pair(chk.id, i));
                }
            }
            printf("\n");
        }
        for(auto& chk :chunk_map) 
            printf("chunk %d loaded by %d \n", chk.first, chk.second);
        printf("\n");

    } else {
        // find unloaded chunk 
        std::vector<int> unloaded_chunk;
        for(int i = 0; i < chunk_size; i++) {
            bool find = false;
            for(auto& iter: chunk_map) 
                if(iter.first == i) {
                    find = true; break; 
                }

            if(!find)
                unloaded_chunk.emplace_back(i);
        }
        if(simple_trace) {
            // Every proc will load simple mesh 
            for(int i = 0; i < proc_size; i++)
                unloaded_chunk.emplace_back(chunk_size);
        }

        printf(" process left chunk first frame\n");
        int proc = 0;
        bool no_proc_empty = false;
        for(auto& chk :unloaded_chunk) {
            if(!no_proc_empty) {
                no_proc_empty = true;
                for(auto& chk_proc: chunk_map) {
                    if(chk_proc.first == -1) {
                        chk_proc.first = chk;
                        no_proc_empty = false;
                        break; 
                    }
                }
            }
            if(no_proc_empty) {
                if(sync)
                    chunk_map.push_back(std::make_pair(chk, -1));
                else
                    chunk_map.push_back(std::make_pair(chk, proc % proc_size));
                proc++;
            }
        }
        for(auto& chk :chunk_map) 
            printf("chunk %d loaded by %d \n", chk.first, chk.second);
        printf("\n");
    }
}

void Scheduler::image_domain_decomposition(bool simple_mesh, bool sync) {
    ImageBlock image(width, height);
    int chunk_size = chunks.size;
    if (block_count == 1) {
        //Allcopy mode
        spp = spp / proc_size;
        for(int i = 0; i < chunks.size; i++) {
            chunk_map.push_back(std::make_pair(i, i));
        }
        return;
    }
    chunk_map.clear();
    get_projected_chunk(proc_rank, proc_size, image);
    //得到聚类中心
    process_left_chunk(proc_size, sync, simple_mesh);
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
//                printf("project p2 %f %f %f\n", p2.x, p2.y, p2.z);
                float depth = p2.z;
                p_min = min(p_min, p2);
                p_max = max(p_max, p2);
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

void Scheduler::split_image_block() {
    ImageBlock image(width, height);

    splat(block_count, &scale[0], 2);
   
    //global available block 
    for(int i = 0; i < proc_size; i++) 
        chunk_map.push_back(std::make_pair(0, proc_rank));
    ImageBlock gblock = project_cube_to_image(camera, BBox(get_bbox()), 0, false, image);
    if(block_count == 1) { 
        render_block = gblock; 
        return;
    }
    
    // int step_width = width / scale[0];
    // int step_height = height / scale[1];   
    int step_width = (gblock.xmax - gblock.xmin) / scale[0];
    int step_height = (gblock.ymax - gblock.ymin) / scale[1];   

    for(int i = 0; i < block_count; i++) {
        int x = i % int(scale[0]);
        int y = i / int(scale[0]);

        int xmin = gblock.xmin + x * step_width;
        int ymin = gblock.ymin + y * step_height;
        ImageBlock block(xmin, ymin, std::min(gblock.xmax, xmin + step_width), std::min(gblock.ymax, ymin + step_height));  
        blockList.emplace_back(block);

    }
    render_block = blockList[proc_rank]; 
    //write_project_result(width, height, domain.data()); 
}

inline void record_send_its(int * send, int self, int its_compact) {
    
    int its = (its_compact & 0xFF) - 1;
    if(its == 254)  send[self] ++; 
    else if(its >= 0)  send[its] ++; 

    its = ((its_compact >> 8) & 0xFF ) - 1;
    if(its == 254)  send[self] ++; else if(its >= 0) send[its] ++;
    
    its = ((its_compact >> 16) & 0xFF ) - 1;
    if(its == 254)  send[self] ++; else if(its >= 0) send[its] ++; 
    
    its = ((its_compact >> 24) & 0xFF ) - 1;
    if(its == 254)  send[self] ++; else if(its >= 0) send[its] ++;
}

inline void record_send_num(int * send, int chk, int ctrb) {
    if((ctrb & 0xFF) > 1)        send[chk] += ((ctrb & 0xFF) >> 1); 
    if((ctrb >> 8) & 0xFF > 1)   send[chk] += (((ctrb >> 8) & 0xFF) >> 1); 
    if((ctrb >> 16) & 0xFF > 1)  send[chk] += (((ctrb >> 16) & 0xFF) >> 1); 
    if((ctrb >> 24) & 0xFF > 1)  send[chk] += (((ctrb >> 24) & 0xFF) >> 1); 
}

inline int get_chunk_id(int x, int y, int z) {
    
}

void Scheduler::chunk_reallocation(int* reduce_buffer) {
    int chunk_size = chunks.size; 
    int res = CHUNK_HIT_RES;
    int size = res * res * 6 * chunk_size;
    int face_size = CHUNK_HIT_RES * CHUNK_HIT_RES;
    int* its = &reduce_buffer[size];
    int* ctrb = reduce_buffer;
    int scale[] = {(int)chunks.scale.x, (int)chunks.scale.y, (int)chunks.scale.z};
    int scale_yz = (int)chunks.scale.y * (int)chunks.scale.z;
    int scale_z = (int)chunks.scale.z;
    
    auto func_get_chunk = [=](int x, int y, int z) -> int { return x*scale_yz + y*scale_z + z;};  
    printf("func_get_chunk xyz %d %d %d chunk %d \n", 1, 2, 3, func_get_chunk(1, 2, 3));

    for(int chk = 0; chk < chunk_size; chk ++) {
        int * send = new int[chunk_size];
        memset(send, 0, sizeof(int) * chunk_size);

        int x = chk / scale_yz ;
        int y = (chk % scale_yz) / scale_z;
        int z = chk % scale_z;
        
        int faces[] = {x-1, x+1, y-1, y+1, z-1, z+1};
        int neighbors[] = {func_get_chunk(x-1, y, z), func_get_chunk(x+1, y, z)
                         , func_get_chunk(x, y-1, z), func_get_chunk(x, y+1, z)
                         , func_get_chunk(x, y, z-1), func_get_chunk(x, y, z+1)
                         };
        printf("chunk %d xyz %d %d %d neighbor: ", chk, x, y, z);

        for(int i = 0; i < 6; i++) {
            if(faces[i] < 0 || faces[i] >= scale[i / 2]) neighbors[i] = -1;
            printf("%d ", neighbors[i]);
        }
        printf("\n");
        int cst = chk * face_size * 6;
        for(int f = 0; f < 6; f++) {
            int neighbor = neighbors[f];
            if(neighbor < 0) continue;
            int fst = cst + f * face_size;
            for(int u = 0; u < res; u++) {  //may be not u, but whatever
                int ust = fst + u * res;
                for(int v = 0; v < res; v++) {
                    int id = ust + v;
                    record_send_num(send, neighbor, ctrb[ust + v]);
                }
            }
        }

        for(int i = 0; i < chunk_size; i++) {
            send_recv_map[chk * chunk_size + i] = send[i];
        }
        delete[] send;
    }
    printf("send recv map\n");
    for(int i = 0; i < chunk_size; i ++) {
        for(int j = 0; j < chunk_size; j ++) {
            printf("%d ", send_recv_map[i * chunk_size + j]);
        }
        printf("\n");
    }
    load_chunk_hit = true;
    printf("chunk_reallocation\n");
}


