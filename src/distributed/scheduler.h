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

static void save_image_ctrb(int* reduce_buffer, int chunk_size, int iter) {
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

    rgb col1_1(0, 100, 25); 
    rgb col1_2(20, 80, 25); 
    rgb col1_3(40, 60, 25); 
    rgb col1_4(60, 40, 25); 
    rgb col1_5(80, 20, 25); 
    rgb col1_6(100, 0, 25); 

    rgb col1(0, 125, 125); 
    rgb col2(255, 0, 0); 
    rgb col3(0, 0, 0); 
    rgb col4(8, 8, 8); 
    //rgb col4(15, 15, 15); 
    int test = 0x1 << 4;
    printf("test %d \n", test );
    for(int i = 0; i < chunk_size; i++) {
        int h_st = res * i; 
        for(int j = 0; j < 6; j++) {
            int w_st = res * j;
            for(int u = 0; u < res; u++) {
                int pu = w_st + u;
                for(int v = 0; v < res; v++) {
                    int pv = h_st + v;
                    
                    int lu = reduce_buffer[all_face_size * i + face_size * j + (res - 1 - v) * res + u];
                    int max_lu = 0;
                    int gradient = 0;
                    int last = 0;
                    for(int k = 0; k < 8; k++) {
                        //if(k < 2) continue;
                        int bit_cur = k * 4;
                        int lu_cur = (lu >> bit_cur) & 0xF;
                        int lu_dir = lu_cur;//median(lu_pre, lu_cur, lu_next); 

                        max_lu = std::max(max_lu, lu_dir);
                        gradient += std::abs(lu_dir - last);
                        last = lu_dir > 0 ? lu_dir : last;
                    }
                    //int t = gradient;//(gradient + max_lu) / 2;
                    int t = max_lu;
                    rgb *col = &col4;

                    img.pixels[4 * (pv * width + pu) + 0] = col->x * t;
                    img.pixels[4 * (pv * width + pu) + 1] = col->y * t;
                    img.pixels[4 * (pv * width + pu) + 2] = col->z * t;
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
    if (!save_png(std::string("picture/chunk_hit_ctrb" + std::to_string(iter++)) + ".png", img))
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
    float time;
    int chunk_size;
    int data_size;
    int cur_chk;
    int iter;
    
    void reset() {
        memset(data, 0, sizeof(int) * data_size * chunk_size);
        time = 0;
    }

    PassRecord(int chunk_size) : chunk_size(chunk_size) {
        data_size = chunk_size + 1;
        data = new int[data_size * chunk_size];
        reset();
        iter = 0;
    }

    ~PassRecord(){
        delete[] data;
    }

    int get_send(int src, int dst) { return data[data_size * src + dst];}

    int get_recv(int src) { return data[data_size * src + chunk_size]; }

    void write_send(float * rays, int size,  bool primary, int chunk) {
        int* cur_chk_send = &data[data_size * chunk]; 
        
        int width = primary ? PRIMARY_WIDTH : SECONDARY_WIDTH;
        int *iptr = (int*) rays;
        for(int i = 0; i < size; i++) {
            int next_chk = iptr[i * width + 9];// & 0xFF;;
            cur_chk_send[next_chk] ++;
        }
    }

    void write_recv(int n, int chunk) {
        data[chunk * data_size + chunk_size] += n;
    }

    void write_time(float t) { time += t; }

    void print(int rank) {
        std::ofstream os = std::ofstream("out/passrecord" + std::to_string(iter++));
        int col = chunk_size + 1;
        for(int i = 0; i < chunk_size; i++) {
            for(int j = 0; j < data_size; j++) {
                os<< data[i * col + j]<<" | ";
            }
            os<<"\n";
        }
    }
    
    // gather pass record, get ray process speed rays/ms 
    float gather() {
        int *reduce_data = new int[data_size * chunk_size];
        MPI_Allreduce(data, reduce_data, data_size * chunk_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
        delete[] data;
        data = reduce_data;
    
        float total_rays = 0;
        for(int i = 0; i < chunk_size; i++) {
            total_rays += data[i * data_size + chunk_size];
        }
        float speed = total_rays / time;
        printf("speed %f\n", speed);
        return speed;
    }
};

int64_t get_chunk_storage_size(int chunk) {
    struct stat statbuf;
    std::string path = "data/"; 
                if (chunk < 100) path += "0";
                if (chunk < 10) path += "0";
                path += std::to_string(chunk) + "/" + "bvh_avx.bin";
                std::cout<<path<<"\n";
    long size = 0;  
    if(stat(path.data(),&statbuf)==0)
        return statbuf.st_size; 
    else
        error("no simple mesh provided\n");
} 

struct LocalChunks {
    struct Chunk {
        std::vector<float> recv_rate; // possible rays from n rays in chunk n 
        int64_t load_time;          
        int64_t storge;
        int id, other_send_left;  //chunk id, init order, possible recv from other proc rays
        float weight;
        int  priority;
        bool loaded;

       // Chunk(int id, float w) :id(id), weight(w) { }
       // Chunk(int id, float w, int64_t ltime) :id(id), weight(w), load_time(ltime) { }
        Chunk(int id, float weight) :id(id), weight(weight) { }
        
        Chunk(int id, float weight, int64_t storge, int64_t load_time) 
            :id(id), weight(weight), storge(storge), load_time(load_time) 
        {
            printf("chunk %d load_time %ld storge %ld\n", id, load_time, storge);
        }
    };
    
    int max_cache_num; 
    std::vector<Chunk> chunks;
    long load_speed;  //b / ms
    float render_speed; //ray msgs / ms
    std::queue<int> history;
    
    bool first_frame;
    int current, next;
    bool new_loaded;
    bool frame_schedule;
    
    LocalChunks() {
        //get load speed;
        new_loaded= true;
        next = -1;
        load_speed = 0;
        render_speed = 0;
        first_frame = true;
        frame_schedule = false;
        max_cache_num = 1;
    }

    void reset(int simple_chunk) {
        chunks.clear(); 
        new_loaded= true;
        next = -1;
        
        if(!SIMPLE_TRACE)
            return;

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

    void insert(int chk, float weight, int64_t chunk_data_size) { 
        //get chunk storage size
        if(load_speed != 0) 
            chunks.push_back(Chunk(chk, weight, chunk_data_size, chunk_data_size / load_speed)); 
        else
            chunks.push_back(Chunk(chk, weight)); 
    }

    int size() { return chunks.size(); }
    void set_new_loaded(bool a) { new_loaded = a; }
    bool get_new_loaded() { return new_loaded; }
    Chunk& get_chunk(int i) { return chunks[i]; }
    
    int find(int chk) {
        for(int i = 0; i < chunks.size(); i++) 
            if(chk == chunks[i].id) return i;
        return -1;
    }
   
    void recv_rays(int size) {
        int i = 0;
        for(i; i < chunks.size(); i++) 
            if(current == chunks[i].id)
                break;
        printf("recv_rays other send left %d %d %d\n", current, chunks[i].other_send_left, size);
        chunks[i].other_send_left -= size;
    }

    float max_load_time() {
        float max = 0;
        for(auto& chk : chunks)
            max = chk.load_time > max ? chk.load_time : max;
    }

    //get priority chunk, new chunk can't be previous chunk
    int get_max_priority_chunk(int cur_chk, RayStreamList * outlist) {
        float min_priority = - max_load_time() - 1;  
        float max = min_priority; 
        int next_chk = -1;
        int local_size = chunks.size(); 
        for(auto& chk : chunks) {
            if(chk.id == cur_chk || outlist[chk.id].ray_size() == 0) { 
                chk.priority = min_priority; 
            } else {
                printf("chk %d : ", chk.id);
                if(SIMPLE_TRACE && !first_frame && frame_schedule) {
                    int predict_ray_size = outlist[chk.id].ray_size() 
                                         + chk.other_send_left / local_size;
                                        // - other local chunk ray size * chk.recv_rate
                   // if(cur_chk >= 0) {
                   //     printf("cur chk >= 0 %d recv rate %f \n", outlist[cur_chk].ray_size(), chk.recv_rate[cur_chk]);
                   //     predict_ray_size += outlist[cur_chk].ray_size() * chk.recv_rate[cur_chk];
                   // }
                    printf("render predict size %d raysize %d local size %d left %d speed %f  load time %ld "
                                , predict_ray_size, outlist[chk.id].ray_size(), local_size, chk.other_send_left, render_speed, chk.load_time);
                    if(chk.loaded) 
                        chk.priority = predict_ray_size / render_speed - min_priority; /*always > 0*/
                    else 
                        chk.priority = predict_ray_size / render_speed - chk.load_time - min_priority; /*always > 0*/

                 //   printf("chunk ray size %d, next size %d recv_rate %f, other left %d  "
                 //           , outlist[chk.id].ray_size(), outlist[cur_chk].ray_size(), chk.recv_rate[cur_chk], chk.recv_rate[cur_chk], chunks[i].other_send_left);
                 //   printf(" chunk %d get priority %f\n", chk.id, priority);
                    chk.priority *= chk.weight;
                } else {
                    int cur_size = outlist[chk.id].ray_size();
                    chk.priority = cur_size;
                }
            }
            if(chk.priority > max) {
                next_chk = chk.id;
                max = chk.priority;
            }
        }
        //if(next_chk == -1)
        //    warn("cur ", cur_chk, " next ", next_chk);
        if(cur_chk == next_chk) 
            return -1;
        
        if(next_chk >= 0) 
            chunks[find(next_chk)].loaded = true;
            
        return next_chk;
    }
   
    void loaded_chunk_size(int& num, int64_t& size) {
        num = 0; size = 0;
        for(auto& chk : chunks) {
            if(chk.loaded) {
                size += chk.storge;
                num++;
            }
        }
    }

    int min_priority_loaded() {
        int min = INT32_MAX; //max_priority; 
        int id = -1;
        for(int i = 0; i < chunks.size(); i++) {
            if(chunks[i].id == current || chunks[i].id == next)
                continue;
            if(chunks[i].loaded && chunks[i].priority < min) {
                min = chunks[i].priority;
                id = i;
            }
        }
        return id;
    }

    int select_new_chunk(RayStreamList * outlist) {
        if(next == -1) { 
            printf("before schedule1  c %d n %d\n", current, next);
            current = get_max_priority_chunk(current, outlist); 
            next = PRELOAD ? get_max_priority_chunk(current, outlist) : -1;
            printf("schedule1  c %d n %d\n", current, next);
        } else {
            if(outlist[next].size() > 0)
                current = next;
            else 
                current = get_max_priority_chunk(current, outlist);
            next = PRELOAD ? get_max_priority_chunk(current, outlist) : -1;
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

        int loaded_num;
        int64_t loaded_size;
        loaded_chunk_size(loaded_num, loaded_size);
        printf("loaded num %d\n", loaded_num);
        while(loaded_num > max_cache_num) {
            int local_idx = min_priority_loaded();
            if(local_idx >= 0) {
                int chk_id = chunks[local_idx].id; 
                //for(int i = 0; i < devNum; i++) 
                rodent_unload_chunk_data(chk_id, 0); 
                printf("unload chunk %d loaded chunk num %d\n", chk_id, loaded_num);
                chunks[local_idx].loaded = false;
            } else {
                break;
            }
            loaded_chunk_size(loaded_num, loaded_size);
        }
    }

    void FIFO_clear() {
        int queue_size = history.size();
        printf("clear history :");
        for(int i = 0; i < queue_size; i++) {
            printf("%d ", history.front()); 
            rodent_unload_chunk_data(history.front(), 0); 
            history.pop();
        }
        printf("\n");
    }

    void FIFO() {
        int queue_size = history.size();
        printf("fifo history size %d :", history.size());
        for(int i = 0; i < queue_size; i++) {
            printf("%d ", history.front()); 
            if(history.front() == current) 
                return;
            history.push(history.front());
            history.pop();
        }
        
        if(queue_size > max_cache_num) {
            int unload_chunk = history.front();
            history.pop(); 
            
            rodent_unload_chunk_data(unload_chunk, 0); 
            printf("unload %d ", unload_chunk);
        }
        history.push(current);
        printf("\n");
    }

    void init(bool sync) {
        current = chunks[0].id;
        chunks[0].loaded = true;
        if(chunks.size() > 1) {
            next    = SIMPLE_TRACE ? chunks[1].id : -1;
            chunks[1].loaded = true;
        }
        new_loaded = true;
        if(sync) 
            FIFO();
    }
    
    bool new_chunk() { return next != current && next != -1; }
    void clear() {
        for(auto& chunk: chunks) {
            printf("clear chunk %d %d\n", chunks.size(), chunk.id);
            //for(int i = 0; i < devNum; i++) 
            rodent_unload_chunk_data(chunk.id, 0); 
        }
        chunks.clear(); 
    }
};

struct ChunkProcs {
    std::vector<int> proc_list;
    int64_t storage;
    int recv, send;
    int iter;
   
    ChunkProcs () { iter = 0; }
    ChunkProcs (ChunkProcs & a) {
        storage = a.storage;
        recv = a.recv;
        send = a.send;
        iter = a.iter;
    }
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

    ChunkManager(int size, int proc_size) : chunk_size(size) {
        chunk_list = new ChunkProcs[chunk_size + 1]; //+1 for simple mesh
        for(int i = 0; i < chunk_size; i++)
            chunk_list[i].storage = get_chunk_storage_size(i);
        //cache n% assigned chunks 
        local_chunks.max_cache_num = CACHE_RATE * chunk_size / proc_size;
        if(SIMPLE_TRACE) 
            local_chunks.max_cache_num++;
        printf("cache size = %d\n", local_chunks.max_cache_num);
    }
    ~ChunkManager() { }

    int switch_current_chunk(RayStreamList * outlist) { 
        local_chunks.select_new_chunk(outlist); 
        scheduler_history.emplace_back(local_chunks.current);
        scheduler_history.emplace_back(local_chunks.next);
    }
    
    bool is_local_chunk(int c) { return local_chunks.find(c) >= 0; }

    int get_proc(int c) {
        if(c >= chunk_size || chunk_list[c].size() == 0) 
            return -1;

        printf("before get chunk %d\n ", chunk_list[c].proc_list.size());
        printf("get chunk %d  proc %d \n", c, chunk_list[c].get_proc());
        return chunk_list[c].get_proc();
    }

    void set_render_speed(float speed) { 
        printf("set render speed %f \n", speed );
        local_chunks.render_speed = speed; 
    }

    bool update_chunk(int chunk, int proc) {
        chunk_list[chunk].set_proc(proc);
    }

    void update_local_chunks(int proc_rank) {
        float weight = 1; 
        // set proc load speed
        for(int i = 0; i < chunk_size + 1; i++) {
            ChunkProcs &chk = chunk_list[i];
            for(auto& proc : chk.proc_list) {
                printf("chunk manager insert %d %d\n", i, proc);
                if(proc == proc_rank) {
                    local_chunks.insert(i, weight, chunk_list[i].storage);
                    weight *= 0.8;
                }
            }
        }
        if(local_chunks.size() == 0)            
            local_chunks.insert(proc_rank, 1, chunk_list[proc_rank].storage);
    }

    void insert(int chk, int proc) { 
        printf("chunk_manager insert %d %d\n", chk, proc);
        chunk_list[chk].proc_list.emplace_back(proc); 
    }
    
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

    int get_next_chunk() {return local_chunks.next;}
    int get_current_chunk() {return local_chunks.current;}
    bool new_chunk(){return local_chunks.new_chunk();}
    void clear_local_chunks() { local_chunks.clear(); }
};

// Grid hierarchy for scheduling
class Scheduler {
public:
    ImageBlock project_cube_to_image(Camera *camera, BBox box, int chunk, bool Record, ImageBlock image);

    void get_neighbor(int chk_id); 
    void get_projected_chunk(int proc_rank, int proc_size, ImageBlock& image, bool sync); 
    void fill_local_chunks(int, bool);
    void data_distributed(bool sync);
    void allcopy(); 
    int  get_most_unloaded_chunk(int *block, int chunk_size);
    
    void preprocess(Camera *, int, bool);    

    int* get_render_block() { return &render_block[0]; }
    void set_render_block(int i) { render_block = blockList[i]; }
    size_t get_block_count() { return block_count; }
    int get_spp() { return spp; }
    void set_load_chunk_hit() { 
        load_chunk_hit = true; 
        chunk_manager->local_chunks.frame_schedule = true;
        chunk_manager->local_chunks.clear();
    }
    Scheduler( int, int, int, int, int);
    ~Scheduler();

    Camera* camera;
    ChunkManager * chunk_manager;
    PassRecord * pass_record;
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

};

Scheduler::Scheduler( int width, int height, int spp, int prank, int psize)
       : width(width), height(height), spp(spp), proc_rank(prank), proc_size(psize)
{
    domain.resize(width * height);
    std::fill(domain.begin(), domain.end(), -1);
    depth.resize(width * height);
    int chunk_size = chunks.size;
    load_chunk_hit = false;
    int all_chunk_size = SIMPLE_TRACE ? chunk_size + 1 : chunk_size;
    chunk_manager = new ChunkManager(all_chunk_size, proc_size);
    pass_record = new PassRecord(all_chunk_size);
}

Scheduler::~Scheduler() {
    delete pass_record;
    delete chunk_manager;
}

int Scheduler::get_most_unloaded_chunk(int *block, int chunk_size) {
    int* chunks = new int [chunk_size];
    for(int i = 0; i< chunk_size; i++)
        chunks[i] = 0;
    for(int y = block[1]; y < block[3]; y++) {
        for(int x = block[0]; x < block[2]; x++) {
            if(domain[y * width + x] >= 0)
                chunks[domain[y * width + x]] ++;
        }
    }
    int max = 0, c = -1;
    int unloaded = -1;
    for(int i = 0; i < chunk_size; i++) {
        if(chunks[i] > max ) {
            max = chunks[i];
            c = i;
        }
    }
    return c;
}

void Scheduler::allcopy() {
    for(int i = 0; i < proc_size; i++)
        chunk_manager->insert(0, i);

    ImageBlock image(width, height);
    ImageBlock gblock = project_cube_to_image(camera, BBox(get_bbox()), 0, false, image);
    int split_image = false;
    if(split_image) {
        
        float length[] = { float(gblock.xmax - gblock.xmin)
                         , float(gblock.ymax - gblock.ymin)} ;
        splat(block_count, &scale[0], &length[0], 2);
        //global available block 
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
    } else {
        render_block = gblock; 
        spp = spp / proc_size;
        for(int i = 0; i < chunks.size; i++)
            chunk_manager->insert(i, i);
    }
}

bool load_speed_cmp(const std::pair<int, int64_t> a, const std::pair<int, int64_t> b) {
    return a.second > b.second;
}

void Scheduler::preprocess(Camera *cam, int block, bool sync) {    
    chunk_manager->clear_local_chunks();
    
    if(load_chunk_hit)
        chunk_manager->set_render_speed(pass_record->gather());
    
    camera = cam;
    block_count = block;
    int chunk_size = chunks.size;
    if(chunk_size == 1)
        allcopy();
    else  
        data_distributed(sync);
    
    printf("generate chunk manager: ");
    LocalChunks &local_chunks = chunk_manager->local_chunks; 
    if(sync || !load_chunk_hit) {
        chunk_manager->update_local_chunks(proc_rank);
    } else {
        local_chunks.reset(SIMPLE_TRACE ? chunk_size-1 : chunk_size);
        int64_t load_speed = local_chunks.load_speed;
        //get local chunk scheduler info: other left, proc_recv..
        // remap chunk procs by load speed
        std::vector<int64_t> proc_load_speed(proc_size);
        MPI_Allgather(&load_speed, 1, MPI_LONG, proc_load_speed.data(), 1, MPI_LONG, MPI_COMM_WORLD);
        std::vector<std::pair<int, int64_t>> load_speed_map;
        for(int i = 0; i < proc_size; i ++) 
            load_speed_map.emplace_back(std::make_pair(i, proc_load_speed[i]));

        sort(load_speed_map.begin(), load_speed_map.end(), load_speed_cmp);
        if(proc_rank == 0)
            for(int i = 0; i < chunk_size + 1; i++) {
                printf(" chunk %d old loaded by: ", i);
                for(auto& proc: chunk_manager->chunk_list[i].proc_list)
                    printf("%d ", proc);
                printf("\n");
            }

        ChunkProcs* old_chunk_list = new ChunkProcs[chunk_size + 1];
        for(int i = 0; i < chunk_size + 1; i++) 
            old_chunk_list[i] = chunk_manager->chunk_list[i];

        for(int p = 0; p < proc_size; p++) {
            int old_proc = load_speed_map[p].first;
            for(int i = 0 ; i < chunk_size + 1; i++) {
                ChunkProcs &chk = old_chunk_list[i];
                for(int j = 0; j < chk.proc_list.size(); j++) {
                    if(chk.proc_list[j] == old_proc) 
                        chunk_manager->chunk_list[i].proc_list[j] = p;
                }
            }
        }
        delete[] old_chunk_list;
        local_chunks.chunks.clear();
        chunk_manager->update_local_chunks(proc_rank);
        chunk_manager->local_chunks.chunks.back().weight = 0.01;

        printf("\nlocal chunks: \n");
        for(int i = 0; i < chunk_size; i++) {
            chunk_manager->chunk_list[i].recv = pass_record->get_recv(i); //recv size
        }

        //
        printf("write local chunk recv rate: ");
        for(auto& chunk : chunk_manager->local_chunks.chunks) {
            int cur_chk = chunk.id;
            printf(" %d : \n", cur_chk);
            chunk.recv_rate.clear();
            chunk.recv_rate.resize(chunk_size + 2);
            chunk.other_send_left = 0;
            for(int i = 0; i < chunk_size; i++) {
                if(i == cur_chk) continue;
                printf("| i %d |\n", i);
                if(chunk_manager->local_chunks.find(i) < 0) {
                    chunk.other_send_left += pass_record->get_send(i, cur_chk); //chunk i send to cur_chk
                    printf("* %d (%d  %d) %f *\n", i, pass_record->get_send(i, cur_chk), chunk_manager->chunk_list[i].recv, chunk.recv_rate[i]);
                } else {
                    if(chunk_manager->chunk_list[i].recv > 0)
                        chunk.recv_rate[i] = pass_record->get_send(i, cur_chk) / chunk_manager->chunk_list[i].recv;
                    printf("| %d (%d  %d) %f |\n", i, pass_record->get_send(i, cur_chk), chunk_manager->chunk_list[i].recv, chunk.recv_rate[i]);
                }
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

    pass_record->reset();
    chunk_manager->local_chunks.init(false /*sync*/); 
}

void Scheduler::get_projected_chunk(int proc_rank, int proc_size, ImageBlock& image, bool sync) {
    int chunk_size = chunks.size; 
    for(int i = 0; i < width * height; i++)
        domain[i] = -1;
    for(int i = 0; i < chunk_size; i++) {
        printf("chunks list %f %f %f max %f %f %f\n", chunks.list[i].min.x, chunks.list[i].min.y, chunks.list[i].min.z
                                                    , chunks.list[i].max.x, chunks.list[i].max.y, chunks.list[i].max.z);
        project_cube_to_image(camera, chunks.list[i], i, true, image);
    }
   
    //global available block
    ImageBlock gblock = project_cube_to_image(camera, chunks.bbox, -1, false, image);
    float length[] = { float(gblock.xmax - gblock.xmin)
                     , float(gblock.ymax - gblock.ymin)} ;
    splat(proc_size, &scale[0], &length[0], 2);
    
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
            if(sync || SIMPLE_TRACE)
                chunk = get_most_unloaded_chunk(&block[0], chunks.size); 
            else 
                chunk = -1;
            printf("%d most unload chunk %d\n", proc_rank, chunk);
        } else if(chunk_size == 1) {
            chunk = 0;
        } else {
            error("chunk size ", chunk_size, " proc size ", proc_size, "invalid");
        }
      
        //if no chunk project to this image block, push (-1, proc_rank)  
        if(chunk > -1)
            chunk_manager->insert(chunk, i);

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

void Scheduler::fill_local_chunks(int proc_size, bool sync) {
    int chunk_size = chunks.size; 
    if(chunk_manager->all_assigned()) return;
    
    ChunkProcs * chunk_list = chunk_manager->chunk_list;
    if(proc_rank == 0)
        pass_record->print(1);
    //load_chunk_hit = false; 
    // put left chunk to proc
    if(!sync && load_chunk_hit) {
        printf(" process left chunk second frame\n");
        std::vector<ChunkHitInfo> chunk_hit_info;
        chunk_hit_info.resize(chunk_size);
        memset(chunk_hit_info.data(), 0, sizeof(ChunkHitInfo) * chunk_size);
        //get recv / send, 
        for(int i = 0; i < chunk_size; i++) {
            chunk_hit_info[i].id = i;
            chunk_hit_info[i].center = -1;
            for(int j = 0; j < chunk_size; j ++) {
                chunk_hit_info[i].send_count += pass_record->get_send(i, j); 
            }
            chunk_hit_info[i].recv_count = pass_record->get_recv(i);
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
        
        // Avoid center empty, push most send chunk to center
        assert(!center.empty());

        for(auto& chk : center) 
            printf("projected cluster center %d\n", chk);
        while(center.size() < proc_size) {
            int max_dis = 0, cand_chk = -1;
            for(int i = 0; i < chunk_size; i++) {
                // Find the most frequent data transfer with center, to fill centers
                std::vector<int>::iterator it = find(center.begin(), center.end(), i);
                if (it != center.end()) continue; // it is already in centers
                int dis = 0;
                for(auto& cent_chk : center) {
                     dis += pass_record->get_send(i, cent_chk); 
                     dis += pass_record->get_send(cent_chk, i); 
                }
                if(dis > max_dis) {
                    max_dis = dis;
                    cand_chk = i;
                }
                printf("%d d %d| ", i, dis);
            }
            printf("\n");
            center.emplace_back(cand_chk);
            chunk_hit_info[cand_chk].center = cand_chk;
        }
        for(auto& chk : center) 
            printf("cluster center %d\n", chk);


        std::vector<int> cluster_size(chunk_size);
        int max_cluster_size = chunk_size / proc_size - 1;
        
        int st = 0;
        int local_size = chunk_size / proc_size - 1;
        for(int i = 0; i < local_size; i++) {
            int ed = st + proc_size;
            for(st; st  < ed; st++) {
                int cent_chk = center[st % proc_size];
                int min_dist = INT32_MAX;
                int cand_chk = -1;
                for(int chk = 0; chk < chunk_size; chk ++) {
                    int dist = 0; 
                    if(chunk_hit_info[chk].center != -1) continue;
                    for(int cluster_chk = 0; cluster_chk < chunk_size; cluster_chk ++) {
                        if(chunk_hit_info[cluster_chk].center == cent_chk) {
                            dist += pass_record->get_send(chk, cluster_chk) + pass_record->get_send(cluster_chk, chk);
                        } 
                    }
                    //dist += pass_record->get_send(chk, cent_chk) + pass_record->get_send(cent_chk, chk);
                    if(dist < min_dist) {
                        min_dist = dist;
                        cand_chk = chk;
                    }
                }
                cluster_size[st % proc_size] ++;
                chunk_hit_info[cand_chk].center = cent_chk;
            }
            st++;
        }

        chunk_manager->reset_chunk_list(); 
        chunk_list = chunk_manager->chunk_list;
        for(int i = 0; i < center.size(); i++) {
            int cent_chk = center[i];
            for(auto& chk: chunk_hit_info) {
                if(chk.center == cent_chk) {
                    chunk_manager->insert(chk.id, i);
                }
            }
        }
    } else {
        // find unloaded chunk 
        if(!sync) {
            int proc = 0;
            for(int i = 0; i < chunk_size; i++) 
                if(chunk_list[i].size() > 0) 
                    proc++;
            for(int i = 0; i < chunk_size; i++) 
                if(chunk_list[i].size() == 0)
                    chunk_manager->insert(i, proc++ % proc_size);
        }
    }
    if(SIMPLE_TRACE) {
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

void Scheduler::data_distributed(bool sync) {
    ImageBlock image(width, height);
    int chunk_size = chunks.size;
    chunk_manager->reset_chunk_list();
  
    get_projected_chunk(proc_rank, proc_size, image, sync);
    
    //put other not projected chunks to procs 
    fill_local_chunks(proc_size, sync);

    //write_project_result(width, height, domain.data()); 
}

ImageBlock Scheduler::project_cube_to_image(Camera *camera, BBox box, int chunk, bool Record, ImageBlock image) {
    if(box.is_inside(camera->eye)) {
        int num_pixels = width * height;
        for(int i = 0; i < num_pixels; i++) {
            if(chunk >=0) {
                domain[i] = chunk;
                depth[i] = 0;
            }
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
    //top-left 0 0 
    ImageBlock block( (p_min.x + 1) * width  / 2
                    , (1 - p_max.y) * height / 2
                    , (p_max.x + 1) * width  / 2
                    , (1 - p_min.y) * height / 2);
   
    if(Record) { 
        for(int y = block.ymin; y < block.ymax; y++) {
            for(int x = block.xmin; x < block.xmax; x++) {
                int id = y * width + x;
                if((domain[id] == -1 || p_min.z < depth[id]) && chunk >= 0) {
                    domain[id] = chunk;
                    depth[id] = p_min.z;
                }
            }
        }
    }

    return block;
}

inline void record_send_num(int * send, int chk, int ctrb) {
    if((ctrb & 0xFF) > 1)        send[chk] += ((ctrb & 0xFF) >> 1); 
    if((ctrb >> 8) & 0xFF > 1)   send[chk] += (((ctrb >> 8) & 0xFF) >> 1); 
    if((ctrb >> 16) & 0xFF > 1)  send[chk] += (((ctrb >> 16) & 0xFF) >> 1); 
    if((ctrb >> 24) & 0xFF > 1)  send[chk] += (((ctrb >> 24) & 0xFF) >> 1); 
}


