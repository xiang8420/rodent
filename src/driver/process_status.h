#pragma once
#include "decomposition.h"
#include <iostream>
#ifndef PI 
#define PI 3.14159265359f
#endif

class ProcStatus {
private:
    int mpi_rank, mpi_size;
    int chunk_size;
    int proc_num, thread_num;
    int loaded_local_chunk_id;
    std::vector<bool> thread_idle;
    std::vector<bool> proc_idle;
    std::vector<int>  local_chunk;
    std::vector<int>  chunk_map;
    bool exit;
    int sent_rays, recv_rays;

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
        if(chunk_size > mpi_size) {
            printf("error, chunk_size > mpi_size %d %d", chunk_size, mpi_size);
            return;
        }
        for(int i = 0; i < chunk_size; i++) {
            chunk_map.emplace_back(i);
        }
        local_chunk.emplace_back(mpi_rank);
        loaded_local_chunk_id = 0;
    }
  
    void thread_reset() {
        thread_idle.resize(thread_num);
        for(int i = 0; i < thread_num; i++) 
            thread_idle[i] = false;
    }
    
    void proc_reset() {
        proc_idle.resize(proc_num);
        for(int i = 0; i < proc_num; i++) 
            proc_idle[i] = false;
        printf("proc reset %d proc idle %d %d\n", proc_num, proc_idle[0], proc_idle[1]);
        std::cout<< proc_idle[0] << proc_idle[1] <<std::endl;
    }
  
    void set_self_idle(){proc_idle[mpi_rank] = true;} 
    
    void set_self_busy(){proc_idle[mpi_rank] = false;}

    void set_proc_busy(int i){proc_idle[i] = false;}
    
    void set_proc_idle(int i){proc_idle[i] = true; }
    
    bool all_proc_idle(){
//        proc_wait[comm->rank] = proc_idle;
        for(int i = 0; i < proc_num; i++){
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
    
    int get_mpi_rank() {return mpi_rank;}
    int get_mpi_size() {return mpi_size;}
    
    int get_sent_ray_count(){return sent_rays;}
    void accumulate_sent(int n){ sent_rays += n;}
    
    int get_recv_ray_count(){return recv_rays;}
    void accumulate_recv(int n){ recv_rays += n;}
   
    bool set_exit() {return exit = true;}
    bool Exit() {return exit;}

    bool all_thread_waiting() {
        for(int i = 0; i < thread_num; i++) {
            if(!thread_idle[i]) 
                return false;
        }
        return true;
    }

    ProcStatus(const float3& e, const float3& d, const float3& u, float fov, int width, int height,
            int spp, int mRank, int mSize, bool image, int cSize, bool rayQ, int pn, int tn) 
        : spp(spp), width(width), height(height), image_decompose(image), rayQ(rayQ), mpi_rank(mRank), 
          mpi_size(mSize), chunk_size(cSize), proc_num(pn), thread_num(tn) 
    {
        
        float ratio = (float) width / (float)height;
        eye = e;
        dir = normalize(d);
        right = normalize(cross(dir, u));
        up = normalize(cross(right, dir));

        w = std::tan(fov * PI / 360.0f);
        h = w / ratio;

        printf("before tile scheduler\n");
        tileScheduler  = new TileScheduler(width, height, mpi_rank, mpi_size);
        printf("after tile scheduler\n");
        thread_reset(); 
        proc_reset();
         
        printf("after tile scheduler\n");
        exit = false;

        for(int i = 0; i < chunk_size; i++) {
            chunk_map.emplace_back(i);
        }
        local_chunk.emplace_back(mpi_rank);
        loaded_local_chunk_id = 0;
        
        sent_rays = 0;
        recv_rays = 0;
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


