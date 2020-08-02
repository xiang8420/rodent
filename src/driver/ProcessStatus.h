#pragma once
#include "interface.h"
#include <iostream>
#include <thread>
#include "decomposition.h"

class ProcStatus {
private:
    int dev_num;
    int chunk_size;
    int cpu_thread_num;
    std::vector<bool> thread_idle;
    int  local_chunk;
    std::vector<int>  chunk_map;
    std::vector<int>  proc_chunk;
    int region[4], proc_spp; 
   
    int buffer_size, buffer_capacity; 
    bool master, exit, load_new_chunk;
public:    
    std::vector<bool> proc_idle;
    int proc_size, proc_rank;
    std::vector<int> global_rays;
    
    struct ImageDecomposition *camera;   
    float3 eye;
    float3 dir;
    float3 right;
    float3 up;
    float w, h;
    int width, height, spp;
    bool image_decompose;
    bool rayQ;
    std::mutex mutex;

    void lock(){mutex.lock();}
    void unlock(){mutex.unlock();}
    
    void set_chunk(int rank, int chunk) {
        proc_chunk[rank] = chunk;

        chunk_map[chunk] = rank;
        if (rank == proc_rank)
            local_chunk = chunk;
    }
  
    void thread_reset() {
        cpu_thread_num = 1;//std::thread::hardware_concurrency() - dev_num + 1;
        printf("\ncpu thread num %d\n", cpu_thread_num);
        thread_idle.resize(cpu_thread_num);
        for(int i = 0; i < cpu_thread_num; i++) 
            thread_idle[i] = false;
    }
    
    void proc_reset() {
        proc_idle.resize(proc_size);
        for(int i = 0; i < proc_size; i++) 
            proc_idle[i] = false;
        printf("proc reset %d proc idle", proc_size);
    //    std::cout<< proc_idle[0] <<"  "<< proc_idle[1] <<std::endl;
    }
 
    bool is_proc_idle(){return proc_idle[proc_rank]; }
    void set_proc_busy(int i){
       // for(int i = 0; i < cpu_thread_num; i++) 
       //     thread_idle[i] = false;
        proc_idle[i] = false;
    }
    void set_proc_idle(int i){proc_idle[i] = true; }
    
    int get_dev_num() {return dev_num;}

    bool all_proc_idle(){
        for(int i = 0; i < proc_size; i++){
            if(!proc_idle[i]) 
                return false;
        } 
        return true;
    }
    
    int get_idle_proc() {
        for(int i = 0; i < proc_size; i++){
            if(proc_idle[i]) 
                return i;
        } 
        return -1;
    }
    int* get_chunk_map(){return chunk_map.data();}
    
    int* get_proc_chunk(){return proc_chunk.data();}

    int get_chunk_size(){return chunk_size;}

    int get_local_chunk(){ return local_chunk;}
   
    int get_proc(int c) {return chunk_map[c];}
    
    int get_chunk(int p) {return proc_chunk[p];}
    
    int get_target_proc(int chunk_id) { return chunk_map[chunk_id];}
    
    void set_thread_idle(int id, bool wait) { thread_idle[id] = wait;}
    
    bool is_thread_idle(int id) { return thread_idle[id]; }
    
    bool isRayQueuing() { return rayQ; }
    
    int get_buffer_capacity(){return buffer_capacity;}

    void set_buffer_capacity(int n){buffer_capacity = n;}
    
    int get_buffer_size(){return buffer_size;}

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
        if (sent <= recv )
            return true;
        else 
            return false;
    }
    
    bool update_global_rays(int* a) {
        int sent = 0, recv = 0;
        for(int i = 0; i < proc_size; i++) {
            global_rays[i] = std::max(global_rays[i], a[i]);
            global_rays[i + proc_size] = std::max(global_rays[i + proc_size], a[i + proc_size]);
            sent += global_rays[i];
            recv += global_rays[i + proc_size];
        }
        printf("updata %d %d\n", sent, recv);
        if (sent <= recv )
            return true;
        else 
            return false;
        return true;
    }

    //update chunk map , return if this proc need load new chunk
    bool update_chunk(int* s) {
        for(int i = 0; i < chunk_size; i++) {
            printf("chunk_size %d chunk_map[i] %d  s[i] %d\n", chunk_size, chunk_map[i], s[i]);
            if(chunk_map[i] != s[i]) {
                chunk_map[i] = s[i];
                if(chunk_map[i] == proc_rank) {
                    local_chunk = i;
                    load_new_chunk = true;
                    printf("update chunk true\n");
                }  
            } 
        } 
        printf("update chunk false\n");
        return load_new_chunk;
    }
    
    bool update_chunk(int candidate_proc, int candidate_chunk) {
        for(int i = 0; i < chunk_size; i++) {
            if(chunk_map[i] == candidate_proc) {
                chunk_map[i] = -1;  //set previous chunk unloaded
            } 
        } 
        chunk_map[candidate_chunk] = candidate_proc;
        printf("candidate chunk %d canidate proc %d\n", candidate_chunk, candidate_proc);
        if(candidate_proc == proc_rank) {
             local_chunk = candidate_chunk;
             load_new_chunk = true;
             return true;
        } 
        return false;
    }
    
    void SetNewChunk(){load_new_chunk = true;}

    void chunk_loaded() {
        thread_reset(); 
        load_new_chunk = false;
    }

    bool has_new_chunk(){return load_new_chunk;}

    void set_exit() {exit = true;}
    
    bool Exit() {return exit;}

    bool all_thread_waiting() {
        for(int i = 0; i < cpu_thread_num; i++) {
            if(!thread_idle[i]) 
                return false;
        }
        return true;
    }

    ProcStatus(const float3 e, const float3 d, const float3 u, float fov, int width, int height,
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
        
        camera  = new ImageDecomposition(e, d, u, fov, width, height, proc_rank, proc_size);
        
        printf("after tile scheduler\n");
        thread_reset(); 
        proc_reset();
         
        printf("after tile scheduler\n");
        exit = false;
        load_new_chunk = true;
      
        //inital chunk map  
        chunk_map.resize(chunk_size);
        for(int i = 0; i < proc_size; i++) 
            chunk_map[i] = i;
        for(int i = proc_size; i < chunk_size; i++) 
            chunk_map[i] = -1;
        local_chunk = proc_rank;
        // 
        proc_chunk.resize(proc_size);
        proc_spp = spp;
        camera->get_camera_info(&region[0], chunk_map.data(), local_chunk, proc_spp, image_decompose); 
        printf("region %d %d %d %d \n chunk map ", region[0], region[1], region[2], region[3]);
        for(int i = 0; i < chunk_size; i++) {
            printf(" %d", chunk_map[i]);
        }
        //
        buffer_size =  1048576;
        buffer_capacity = (buffer_size & ~((1 << 5) - 1)) + 32; // round to 32
        
        global_rays.resize(proc_size * proc_size);
    }

    int get_tile_info(int* res, float* processTime) {
        for(int i = 0; i < proc_size ;i++)
            res[i] = region[i];
        return proc_spp; 
    }
    
    void camera_rotate(float yaw, float pitch) {
        camera->rotate(yaw, pitch);
    }

    void camera_move(float x, float y, float z) {
        eye += right * x + up * y + dir * z;
    }
};


