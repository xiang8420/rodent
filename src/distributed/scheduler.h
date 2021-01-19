#pragma once
#include <vector>
#include <memory>
#include <cstring>
#include <mutex>
#include <utility>
#include <sys/stat.h>

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
    float * time;
    int chunk_size;
    int cur_chk;
    
    void reset() {
        memset(data, 0, sizeof(int) * (chunk_size + 1) * chunk_size);
        memset(reduce_data, 0, sizeof(int) * (chunk_size + 1) * chunk_size);
        memset(time, 0, sizeof(float) * chunk_size);
    }

    PassRecord(int chunk_size) : chunk_size(chunk_size) {
        data = new int[(chunk_size + 1) * chunk_size];
        reduce_data = new int[(chunk_size + 1) * chunk_size]; 
        time = new float[chunk_size]; 
        reset();
    }

    ~PassRecord(){
        delete[] data;
        delete[] reduce_data;
        delete[] time;
    }

    void set_cur_chk(int chk) {
        cur_chk = chk;
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

    void write_time(float t) { time[cur_chk] += t; }

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
    
    // gather pass record, get ray process speed rays/ms 
    float gather() {
       // MPI_Allreduce(data, reduce_data, (chunk_size + 1) * chunk_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
       // int* tmp = data;
       // data = reduce_data;
       // reduce_data = tmp;
        float *reduce_time = new float[chunk_size];
        printf("pass record time before gather: ");
        for(int i = 0; i < chunk_size; i++)
            std::cout<<time[i]<<" | ";
    //    MPI_Allreduce(time, reduce_time, chunk_size, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); 
    //    delete[] time;
    //    time = reduce_time;
        
        float total_time = 0;
        float total_rays = 0;
        
        printf("\nchunk speed : ");
        for(int i = 0; i < chunk_size; i++) {
            printf(" %f %f |", data[i * (chunk_size + 1)+ chunk_size], time[i]);
            total_time += time[i];
            total_rays += data[i * (chunk_size + 1)+ chunk_size];
        }
        float speed = total_rays / total_time;
        printf("speed %f\n", speed);
        return speed;
    }
};

int64_t get_chunk_storage_size(int chunk) {
    struct stat statbuf;
    std::string path = "data/"; 
                path += ((chunk > 9) ? "0" : "00") 
                     + std::to_string(chunk) + "/"
                     + "bvh_avx.bin";
                std::cout<<path<<"\n";
    long size = 0;  
    if(stat(path.data(),&statbuf)==0)
        return statbuf.st_size; 
    else
        error("no simple mesh provided\n");
} 

//(outlist[next].rays_size() 
//+ next.other_send_left / local_chunk_size 
//+ next.recv_map[current] * outlist[current].rays_size() 
//- for(outlist[i] * recv_map[i]) / local_chunk_size
//) / compute_speed(计算时间) - load_time(加载时间) 希望数据加载被隐藏,所以这个值越大越会被换入 
// make sure the formula > 0
// simple_mesh always stay in memory
struct Chunk {
    int id;
    int priority;
    int load_time;          //loaded time = bvh_size / bandwidth
    int other_send_left;    //recv from other proc (ray msg) 
    std::vector<float> recv_rate; //recv from other local chunk, possiable recv size = outlist[n].size() * recv_rate[n]

    Chunk(int id, int priority) :id(id), priority(priority) { }

    Chunk(int id, int priority, int64_t ltime) :id(id), priority(priority), load_time(ltime) {
        printf("new chunk %d pri %d load time %ld\n", id, priority, ltime);
    }
    
//    int64_t get_runtime_priotity(int ) {
//        return 
//    }


};

struct LocalChunks {
    std::vector<Chunk> chunks;
    long load_speed;  //b / ms
    float render_speed; //ray msgs / ms
    std::queue<int> history;

    bool first_frame;
    int current, next;
    bool new_loaded;
    
    LocalChunks() {
        //get load speed;
        new_loaded= true;
        next = -1;
        load_speed = 0;
        render_speed = 0;
        first_frame = true;
    }

    void reset(int simple_chunk) {
        chunks.clear(); 
        new_loaded= true;
        next = -1;

        if(load_speed == 0) {
            // first frame 
            statistics.start("preload simple mesh");
            load_bvh_test();        
            float t = statistics.end("preload simple mesh");
            load_speed = get_chunk_storage_size(simple_chunk) / t;
            printf("test load bvh speed %ld\n", load_speed);
        } else {
            // second frame
            first_frame = false;
        }
    }

    void insert(int chk, int pri, int64_t chunk_size) { 
        //get chunk storage size
        if(load_speed != 0) 
            chunks.push_back(Chunk(chk, pri, chunk_size / load_speed)); 
        else
            chunks.push_back(Chunk(chk, pri)); 
    }

    int size() { return chunks.size(); }
    void set_new_loaded(bool a) { new_loaded = a; }
    bool get_new_loaded() { return new_loaded; }
    Chunk& get_chunk(int i) { return chunks[i]; }
    
    bool find(int chk) {
        for(int i = 0; i < chunks.size(); i++) 
            if(chk == chunks[i].id) return true;
        return false;
    }
   
    void recv_rays(int size) {
        int i = 0;
        for(i; i < chunks.size(); i++) 
            if(current == chunks[i].id)
                break;
        printf("recv_rays other send left %d %d %d\n", current, chunks[i].other_send_left, size);
        chunks[i].other_send_left -= size;
    }

    void record_history() {
        if(history.size() >= chunks.size()) 
            history.pop();
        history.push(current); 
    }

    //get priority chunk, new chunk can't be previous chunk
    int get_max_priority_chunk(int cur_chk, RayStreamList * outlist) {
        float weight = 1;
        float max = -10000000; 
        int next_chk = -1;
        int local_size = chunks.size(); 
        for(auto& chk : chunks) {
            float priority;
            if(chk.id == cur_chk || outlist[chk.id].ray_size() == 0) continue;
            printf("chk %d :\n", chk.id);
            if(SIMPLE_TRACE && !first_frame) {
                int predict_ray_size = outlist[chk.id].ray_size() 
                                     + chk.other_send_left / local_size;
                                    // - other local chunk ray size * chk.recv_rate
                if(cur_chk >= 0)
                    predict_ray_size += outlist[cur_chk].ray_size() * chk.recv_rate[cur_chk];
                printf("render predict size %d speed %f  load time %ld ", predict_ray_size, render_speed, chk.load_time);
                priority = predict_ray_size / render_speed - chk.load_time; 

             //   printf("chunk ray size %d, next size %d recv_rate %f, other left %d  "
             //           , outlist[chk.id].ray_size(), outlist[cur_chk].ray_size(), chk.recv_rate[cur_chk], chk.recv_rate[cur_chk], chunks[i].other_send_left);
                printf(" chunk %d get priority %f\n", chk.id, priority);
            } else {
                int cur_size = outlist[chk.id].ray_size();
                priority = cur_size;
            }
            priority *= weight;
            weight *= 0.8; // 

            if(priority > max) {
                next_chk = chk.id;
                max = priority;
            }
        }
        if(next_chk == -1)
            printf("cur %d next %d %d %d\n", cur_chk, next_chk, chunks.back().id, outlist[chunks.back().id].ray_size());
        if(cur_chk == next_chk) 
            return -1;
        return next_chk;
    }
    int select_new_chunk(RayStreamList * outlist) {
        if(next == -1) { 
            current = get_max_priority_chunk(current, outlist); 
            next = get_max_priority_chunk(current, outlist); 
            printf("schedule1  c %d n %d\n", current, next);
        } else {
            if(outlist[next].size() > 0) 
                current = next;
            else 
                current = get_max_priority_chunk(current, outlist);
            next = get_max_priority_chunk(current, outlist); 
            printf("schedule2 c %d n %d\n", current, next);
        }
        
        if(next == -1 && current == -1) {
            printf("check new chunk %d ", current);
            for(int i = 0; i < chunks.size(); i ++) {
                printf("%d chk %d size %d", i, chunks[i].id, outlist[i].ray_size());
            }
        }
        new_loaded = true;
        if(!first_frame)
            printf("select current %d rays %d next %d rays %d \n", current, outlist[current].size(), next, outlist[next].size());
    }

    void init() {
        current = chunks[0].id;
        next = chunks[1].id;
        printf("init current %d priority %f  next %d priority %f\n "
                , current, chunks[0].priority, next, chunks[1].priority);
        record_history(); 
        new_loaded = true;
    }
    bool new_chunk() { return next != current && next != -1; }
};

struct ChunkProcs {
    std::vector<int> proc_list;
    int64_t storage;
    int recv, send;
    int iter;
    
    ChunkProcs () { iter = 0; }
    int get_proc() { return proc_list[iter++ % proc_list.size()]; }
    int set_proc(int proc) { proc_list.push_back(proc); }
    int size() {return proc_list.size(); }
    int back() {return proc_list.back(); }
    bool loaded() {return !proc_list.empty(); }
};

struct ChunkManager {
    ChunkProcs * chunk_list;
    LocalChunks local_chunks;  // chunk id, chunk priority
    int64_t chunk_size;
    std::vector<int> scheduler_history; //current, next

    ChunkManager(int size) : chunk_size(size) {
        chunk_list = new ChunkProcs[chunk_size + 1]; //+1 for simple mesh
        for(int i = 0; i < chunk_size; i++)
            chunk_list[i].storage = get_chunk_storage_size(i);
    }
    ~ChunkManager() { }

    int switch_current_chunk(RayStreamList * outlist) { 
        local_chunks.select_new_chunk(outlist); 
        scheduler_history.emplace_back(local_chunks.current);
        scheduler_history.emplace_back(local_chunks.next);
    }
    
    bool is_local_chunk(int c) { return local_chunks.find(c); }

    int get_proc(int c) {
        if(c >= chunk_size) { 
            warn(c, " chunk not loaded");
            return -1;
        }
        printf("get chunk %d  proc %d \n", c, chunk_list[c].get_proc());
        return chunk_list[c].get_proc();
    }

    void set_render_speed(float speed) { 
        printf("set render speed %f \n", speed );
        local_chunks.render_speed = speed; }

    bool update_chunk(int* schedule, int proc_rank) { error("function update chunk not implement"); }

    bool update_chunk(int candidate_proc, int candidate_chunk, int proc_rank) {
        error("check function update chunk first.");
//        chunk_list[candidate_chunk].set_proc(candidate_proc);
//        if(candidate_proc == proc_rank) {
//            local_chunks.insert(candidate_chunk, 1/*priority*/);
//            local_chunks.set_new_loaded(true);
//        }
//        return local_chunks.get_new_loaded();
    }

    void update_local(int proc_rank) {
        int priority = 0;
        // set proc load speed
        local_chunks.reset(SIMPLE_TRACE ? chunk_size-1 : chunk_size);

        for(int i = 0; i < chunk_size + 1; i++) {
            ChunkProcs &chk = chunk_list[i];
            for(auto& proc : chk.proc_list)
                if(proc == proc_rank)
                    local_chunks.insert(i, priority++/*priority*/, chunk_list[i].storage);
        }
    }
    
    void insert(int chk, int proc) { chunk_list[chk].proc_list.emplace_back(proc); }
    
    bool all_assigned() {
        for(int i = 0; i < chunk_size; i++) 
            if(chunk_list[i].size() == 0)
                return false;
        return true;
    }
    void reset_chunk_list() {
        delete[] chunk_list;
        chunk_list = new ChunkProcs[chunk_size + 1]; //+1 for simple mesh
    }

    void print_schedule(int frame, int rank) {
        std::ofstream os = std::ofstream("out/schedule_"+ std::to_string(frame) + "_" + std::to_string(rank));
        printf("local chunk : ");
        for(int i = 0; i < local_chunks.chunks.size(); i++) {
            os << " " << local_chunks.chunks[i].id;
        }
        os<<"\n";
        for(int i = 0; i < scheduler_history.size(); i += 2) {
            os<<" c "<< scheduler_history[i] <<" n "<< scheduler_history[i + 1] <<" \n";
        } 
        os<<"\n";
        
        scheduler_history.clear();
    }

    void get_start_chunk() {local_chunks.init(); printf("init %d %d\n", local_chunks.current, local_chunks.next);}
    int get_next_chunk() {return local_chunks.next;}
    int get_current_chunk() {return local_chunks.current;}
    bool new_chunk(){return local_chunks.new_chunk();}
};

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

// Grid hierarchy for scheduling
class Scheduler {
public:
    ImageBlock project_cube_to_image(Camera *camera, BBox box, int chunk, bool Record, ImageBlock image);

    void get_projected_chunk(int proc_rank, int proc_size, ImageBlock& image); 
    void fill_local_chunks(int, bool, bool);
    void image_domain_decomposition(bool simple_mesh, bool sync);
    void split_image_block();
    void chunk_reallocation(int* reduce_buffer); 
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
    ImageBlock render_block;
    MeshChunk chunks;

    // project domain to image space
    std::vector<int>   domain;
    std::vector<float> depth;
    std::vector<ImageBlock> blockList;

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
        if(!chunk_manager->chunk_list[i].loaded()) {
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
    int chunk_size = chunks.size;
    if(block_count == proc_size) 
        image_domain_decomposition(simple_trace, sync);
    else
        split_image_block(); 
    
    printf("generate chunk manager: ");
    chunk_manager->update_local(proc_rank); 
    //get local chunk scheduler info: other left, proc_recv..
    if(load_chunk_hit) {
        printf("\nlocal chunks: \n");
        for(int i = 0; i < chunk_size; i++) {
            chunk_manager->chunk_list[i].recv = 0;
            for(int j = 0; j < chunk_size; j++)
                chunk_manager->chunk_list[i].recv += send_recv_map[j * chunk_size + i]; //recv size
        }  

        for(auto& chunk : chunk_manager->local_chunks.chunks) {
            int cur_chk = chunk.id;
            printf(" %d : ", cur_chk);
            chunk.other_send_left = 0;
            for(int i = 0; i < chunk_size; i++) {
                if(i == cur_chk) continue;
                if(!chunk_manager->local_chunks.find(i)) {
                    chunk.other_send_left += send_recv_map[chunk_size * i + cur_chk];  //chunk i send to cur_chk
                } 
                chunk.recv_rate.emplace_back(float(send_recv_map[chunk_size * i + cur_chk]) / float(chunk_manager->chunk_list[i].recv));
            }
            printf("other send left %d \n", chunk.other_send_left);
        }
        for(int i = 0; i < chunk_size + 1; i++) {
            printf("chunk %d loaded by: ", i);
            for(auto& proc: chunk_manager->chunk_list[i].proc_list)
                printf("%d ", proc);
            printf("\n");
        } 
        printf("local chunk  : ");
        for(auto& chunk : chunk_manager->local_chunks.chunks) {
            printf(" %d recv rate :", chunk.id);
            for(int i = 0; i < chunk_size; i++) {
                printf(" %f |", chunk.recv_rate[i]);
            }
            printf("\n");
        }
        printf("\n");
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
        if(chunk >= -1)
            chunk_manager->insert(chunk, i);

        blockList.emplace_back(block);
        
        if(i == proc_rank) {
            render_block = block; 
        }
    }
}

void Scheduler::fill_local_chunks(int proc_size, bool sync, bool simple_trace) {
    
    int chunk_size = chunks.size; 
    if(chunk_manager->all_assigned()) return;
    
    ChunkProcs * chunk_list = chunk_manager->chunk_list;

    //load_chunk_hit = false; 
    // put left chunk to proc
    if(!sync && load_chunk_hit) {
//    if(false) {
        printf(" process left chunk second frame\n");
        std::vector<ChunkHitInfo> chunk_hit_info;
        chunk_hit_info.resize(chunk_size);
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
        //Get cluster center
        std::vector<int> center;
        for(int i = 0; i < chunk_size; i++) {
            if(chunk_list[i].size() > 0) {
                chunk_hit_info[i].center = i;
                center.emplace_back(i);
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
   
        chunk_manager->reset_chunk_list(); 
        chunk_list = chunk_manager->chunk_list;
        for(int i = 0; i < center.size(); i++) {
            int cent_chk = center[i];
            printf("center %d | ", cent_chk);
            for(auto& chk: chunk_hit_info) {
                if(chk.center == cent_chk) {
                    printf("%d ", chk.id);
                    chunk_manager->insert(chk.id, i);
                }
            }
            printf("\n");
        }
    } else {
        // find unloaded chunk 
        printf(" process left chunk first frame\n");
        if(!sync) {
            int proc = 0;
            for(int i = 0; i < chunk_size; i++) 
                if(chunk_list[i].size() == 0)
                    chunk_manager->insert(i, proc++ % proc_size);
        }
    }
    if(simple_trace) {
        // Every proc will load simple mesh 
        for(int i = 0; i < proc_size; i++)
            chunk_manager->insert(chunk_size, i);
    }
    for(int i = 0; i < chunk_size + 1; i++) {
        printf("chunk %d loaded by: ", i);
        for(auto& proc: chunk_list[i].proc_list)
            printf("%d ", proc);
        printf("\n");
    }
}

void Scheduler::image_domain_decomposition(bool simple_mesh, bool sync) {
    ImageBlock image(width, height);
    int chunk_size = chunks.size;
    chunk_manager->reset_chunk_list();
    if (block_count == 1) {
        //Allcopy mode
        spp = spp / proc_size;
        for(int i = 0; i < chunks.size; i++)
            chunk_manager->insert(i, i);
        return;
    }
    get_projected_chunk(proc_rank, proc_size, image);
    //put other not projected chunks to procs 
    fill_local_chunks(proc_size, sync, simple_mesh);
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
        chunk_manager->insert(0, i);
    ImageBlock gblock = project_cube_to_image(camera, BBox(get_bbox()), 0, false, image);
    if(block_count == 1) { 
        render_block = gblock; 
        return;
    }
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
}

inline void record_send_num(int * send, int chk, int ctrb) {
    if((ctrb & 0xFF) > 1)        send[chk] += ((ctrb & 0xFF) >> 1); 
    if((ctrb >> 8) & 0xFF > 1)   send[chk] += (((ctrb >> 8) & 0xFF) >> 1); 
    if((ctrb >> 16) & 0xFF > 1)  send[chk] += (((ctrb >> 16) & 0xFF) >> 1); 
    if((ctrb >> 24) & 0xFF > 1)  send[chk] += (((ctrb >> 24) & 0xFF) >> 1); 
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

    memset(send_recv_map.data(), 0, sizeof(int) * chunk_size * chunk_size);
    for(int chk = 0; chk < chunk_size; chk ++) {
        int * send = &send_recv_map[chk * chunk_size];

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


