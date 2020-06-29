#pragma once
#include "decomposition.h"
#include <iostream>
#include <thread>
#ifndef PI 
#define PI 3.14159265359f
#endif

class ProcStatus {
private:
    int proc_size, proc_rank;
    int dev_num;
    int chunk_size;
    int cpu_thread_num;
    int loaded_local_chunk_id;
    std::vector<bool> thread_idle;
    std::vector<bool> proc_idle;
    std::vector<int>  local_chunk;
    std::vector<int>  chunk_map;
    
    std::vector<int> global_rays;
    
    bool master, exit;
public:    
    float3 eye;
    float3 dir;
    float3 right;
    float3 up;
    float w, h;
    struct TileScheduler *tileScheduler;   
    int width, height, spp;
    bool image_decompose;
    bool rayQ;
    std::mutex mutex;

    void lock(){mutex.lock();}
    void unlock(){mutex.unlock();}
    
    
    void set_chunks(){
        //get chunk infomation from file
        if(chunk_size > proc_size) {
            printf("error, chunk_size > proc_size %d %d", chunk_size, proc_size);
            return;
        }
        for(int i = 0; i < chunk_size; i++) {
            chunk_map.emplace_back(i);
        }
        local_chunk.emplace_back(proc_rank);
        loaded_local_chunk_id = 0;
    }
  
    void thread_reset() {
        cpu_thread_num = std::thread::hardware_concurrency() - dev_num + 1;
        printf("\ncpu thread num %d\n", cpu_thread_num);
        thread_idle.resize(cpu_thread_num);
        for(int i = 0; i < cpu_thread_num; i++) 
            thread_idle[i] = false;
    }
    
    void proc_reset() {
        proc_idle.resize(proc_size);
        for(int i = 0; i < proc_size; i++) 
            proc_idle[i] = false;
        printf("proc reset %d proc idle %d %d\n", proc_size, proc_idle[0], proc_idle[1]);
        std::cout<< proc_idle[0] << proc_idle[1] <<std::endl;
    }
  
    void set_self_idle(){proc_idle[proc_rank] = true;} 
    
    void set_self_busy(){proc_idle[proc_rank] = false;}

    void set_proc_busy(int i){proc_idle[i] = false;}
    
    void set_proc_idle(int i){proc_idle[i] = true; }
    
    int get_dev_num() {return dev_num;}

    bool all_proc_idle(){
//        proc_wait[comm->rank] = proc_idle;
        for(int i = 0; i < proc_size; i++){
            if(!proc_idle[i]) 
                return false;
        } 
        return true;
    }

    int get_chunk_size(){return chunk_size;}

    int get_loaded_chunk(){ return local_chunk[loaded_local_chunk_id];}
    
    std::vector<int>& get_local_chunk(){ return local_chunk;}

    int get_new_chunk() {
        loaded_local_chunk_id = (loaded_local_chunk_id + 1) % local_chunk.size();
        return local_chunk[loaded_local_chunk_id];
    }

    int get_dst_proc(int chunk_id) { return chunk_map[chunk_id];}

    void set_thread_idle(int id, bool wait) { thread_idle[id] = wait;}
    
    bool is_thread_idle(int id) { return thread_idle[id]; }
    
    bool isRayQueuing() { return rayQ; }
    
    int get_rank() {return proc_rank;}
    
    int get_size() {return proc_size;}
    
//    int get_sent_ray_count(int n){ return sent_rays; }

//    int get_recv_ray_count(int n) { return recv_rays; }
    
    void accumulate_sent(int n) { global_rays[proc_rank] += n; }
    
    void accumulate_recv(int n) { global_rays[proc_rank + proc_size] += n;}
 
    int* get_status() {return global_rays.data();}

    bool all_rays_received() {
        int sent = 0, recv = 0;
        for(int i = 0; i < proc_size; i++) {
            sent += global_rays[i];
            recv += global_rays[i + proc_size];
        }
        printf("check %d %d\n", sent, recv);
        if (sent * 0.99 <= recv )
            return true;
        else 
            return false;
    }

    bool updata_global_rays(int* a) {
        int sent = 0, recv = 0;
        for(int i = 0; i < proc_size; i++) {
            global_rays[i] = std::max(global_rays[i], a[i]);
            global_rays[i + proc_size] = std::max(global_rays[i + proc_size], a[i + proc_size]);
            sent += global_rays[i];
            recv += global_rays[i + proc_size];
        }
        printf("updata %d %d\n", sent, recv);
        if (sent * 0.99 <= recv )
            return true;
        else 
            return false;
        return true;
    }


    bool set_exit() {return exit = true;}
    
    bool Exit() {return exit;}

    bool all_thread_waiting() {
        for(int i = 0; i < cpu_thread_num; i++) {
            if(!thread_idle[i]) 
                return false;
        }
        return true;
    }

    ProcStatus(const float3& e, const float3& d, const float3& u, float fov, int width, int height,
            int spp, int rank, int size, bool image, int cSize, bool rayQ, int dev, bool master) 
        : spp(spp), width(width), height(height), image_decompose(image), rayQ(rayQ), proc_rank(rank), 
          chunk_size(cSize), dev_num(dev), master(master) 
    {
        
        float ratio = (float) width / (float)height;
        eye = e;
        dir = normalize(d);
        right = normalize(cross(dir, u));
        up = normalize(cross(right, dir));

        w = std::tan(fov * PI / 360.0f);
        h = w / ratio;
        
        proc_size = master ? size - 1 : size;
        printf("before tile scheduler\n");
        tileScheduler  = new TileScheduler(width, height, proc_rank, proc_size);
        printf("after tile scheduler\n");
        thread_reset(); 
        proc_reset();
         
        printf("after tile scheduler\n");
        exit = false;

        for(int i = 0; i < chunk_size; i++) {
            chunk_map.emplace_back(i);
        }
        local_chunk.emplace_back(proc_rank);
        loaded_local_chunk_id = 0;
        
        global_rays.resize(proc_size * proc_size);
    }

    void camera_rotate(float yaw, float pitch) {
        dir = ::rotate(dir, right,  -pitch);
        dir = ::rotate(dir, up,     -yaw);
        dir = normalize(dir);
        right = normalize(cross(dir, up));
        up = normalize(cross(right, dir));
    }

    void camera_move(float x, float y, float z) {
        eye += right * x + up * y + dir * z;
    }
    
    int get_tile_info(int* region, float* processTime) {
        return tileScheduler->get_camera_info(region, spp, processTime, image_decompose); 
    }

};


