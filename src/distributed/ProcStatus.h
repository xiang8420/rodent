#pragma once
#include "../driver/interface.h"
#include <iostream>
#include <thread>
#include "decomposition.h"

class ProcStatus {
private:
    //thread 
    std::vector<bool> thread_idle;
    int cpu_thread_num;
    int dev_num;

    //domain settings
    int chunk_size, local_chunk;
    std::vector<int>  chunk_map;

    //decomposition
    ImageDecomposition *camera;   

    std::mutex mutex;
   
    int buffer_size, buffer_capacity; 
    bool master, exit, load_new_chunk;
    int size, rank;
    std::vector<bool> proc_idle;
public:    
    // camera setting 
    bool image_decompose;
    int region[4]; 
    int width, height, spp;
    
    std::vector<int>  global_rays;

public:
    void lock(){mutex.lock();}
    
    void unlock(){mutex.unlock();}
    
    int* get_chunk_map(){return chunk_map.data();}
    
    int get_chunk_size(){return chunk_size;}

    int get_local_chunk(){ return local_chunk;}
   
    int get_proc(int c) {return chunk_map[c];}
    
    void set_thread_idle(int id, bool wait) { thread_idle[id] = wait;}
    
    bool is_thread_idle(int id) { return thread_idle[id]; }
    
    int get_buffer_capacity(){return buffer_capacity;}

    void set_buffer_capacity(int n){buffer_capacity = n;}
    
    int get_buffer_size(){return buffer_size;}

    void accumulate_sent(int n) { global_rays[rank] += n; }
    
    void accumulate_recv(int n) { global_rays[rank + size] += n;}
 
    int* get_status() {return global_rays.data();}
    
    bool is_proc_idle(){return proc_idle[rank]; }
    
    void set_proc_busy(int i){ proc_idle[i] = false; }
    
    void set_proc_idle(int i){proc_idle[i] = true; }
    
    int get_dev_num() {return dev_num;}
    
    ImageDecomposition* get_camera() { return camera; } 
    
    void set_new_chunk(){load_new_chunk = true;}
    
    
    void set_chunk(int rank, int chunk) {
        chunk_map[chunk] = rank;
        if (rank == rank)
            local_chunk = chunk;
    }
  
    void thread_reset() {
        cpu_thread_num = 1;//std::thread::hardware_concurrency();
        printf("\ncpu thread num %d\n", cpu_thread_num);
        thread_idle.resize(cpu_thread_num);
        for(int i = 0; i < cpu_thread_num; i++) 
            thread_idle[i] = false;
    }
    
    void proc_reset() {
        proc_idle.resize(size);
        for(int i = 0; i < size; i++) 
            proc_idle[i] = false;
        printf("proc reset %d proc idle", size);
    //    std::cout<< proc_idle[0] <<"  "<< proc_idle[1] <<std::endl;
    }
 

    bool all_proc_idle(){
        for(int i = 0; i < size; i++){
            if(!proc_idle[i]) 
                return false;
        } 
        return true;
    }
    
    int get_idle_proc() {
        for(int i = 0; i < size; i++){
            if(proc_idle[i]) 
                return i;
        } 
        return -1;
    }

    bool all_rays_received() {
        int sent = 0, recv = 0;
        for(int i = 0; i < size; i++) {
            sent += global_rays[i];
            recv += global_rays[i + size];
        }
        printf("check %d %d\n", sent, recv);
        if (sent <= recv )
            return true;
        else 
            return false;
    }
    
    bool update_global_rays(int* a) {
        int sent = 0, recv = 0;
        for(int i = 0; i < size; i++) {
            global_rays[i] = std::max(global_rays[i], a[i]);
            global_rays[i + size] = std::max(global_rays[i + size], a[i + size]);
            sent += global_rays[i];
            recv += global_rays[i + size];
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
                if(chunk_map[i] == rank) {
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
        if(candidate_proc == rank) {
             local_chunk = candidate_chunk;
             load_new_chunk = true;
             return true;
        } 
        return false;
    }
    

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
    
    ~ProcStatus(){
        printf("procstatus delete\n");
    }

    ProcStatus(const float3 e, const float3 d, const float3 u, float fov, int width, int height,
            int spp_global, int comm_rank, int comm_size, int cSize, int dev, bool master) 
        : width(width), height(height), rank(comm_rank), chunk_size(cSize), dev_num(dev), master(master) 
    {
        
        size = master ? comm_size - 1 : comm_size;
        printf("before tile scheduler\n");
        
        camera  = new ImageDecomposition(e, d, u, fov, width, height, rank, size);
        
        printf("after tile scheduler\n");
        thread_reset(); 
        proc_reset();
         
        printf("after tile scheduler\n");
        exit = false;
      
        //inital chunk map  
        chunk_map.resize(chunk_size);
        // 
        spp = spp_global;
        camera->get_camera_info(&region[0], chunk_map.data(), local_chunk, spp, true); 
        load_new_chunk = true;
        
        printf("proc status region %d %d %d %d \n chunk map ", region[0], region[1], region[2], region[3]);
        for(int i = 0; i < chunk_size; i++) {
            printf(" %d", chunk_map[i]);
        }
        //
        buffer_size = 1024 * 1024 / cpu_thread_num;
        buffer_capacity = (buffer_size & ~((1 << 5) - 1)) + 32; // round to 32
        
        global_rays.resize(size * size);
    }

    int get_tile_info(int* res, float* processTime) {
        
        printf("proc status get tile info region %d %d %d %d \n chunk map ", region[0], region[1], region[2], region[3]);
        for(int i = 0; i < 4; i++)
            res[i] = region[i];
        printf("proc status get tile info region %d %d %d %d \n chunk map ", res[0], res[1], res[2], res[3]);
        return spp; 
    }
    
    void camera_rotate(float yaw, float pitch) {
        camera->rotate(yaw, pitch);
    }

    void camera_move(float x, float y, float z) {
        camera->move(x, y, z);
    }
};


