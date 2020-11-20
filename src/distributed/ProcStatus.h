#pragma once
#include "../driver/interface.h"
#include <iostream>
#include <thread>
#include "decomposition.h"

#define MAX_CHUNK 3
#define MAX_PROC  128

class ProcStatus {
private:
    //thread 
    std::vector<bool> thread_idle;
    int work_thread_num;
    int dev_num;

    //domain settings
    int chunk_proc[MAX_CHUNK * MAX_PROC];
    std::vector<int> local_chunks;
    int chunk_size, current_chunk;

    int stream_logic_capacity, stream_store_capacity;
    int out_stream_capacity; 
    bool exit, load_new_chunk;
    int proc_size, proc_rank;
    std::vector<bool> proc_idle;

    bool rough_trace;
public:    
    ImageDecomposition *camera;   
    std::vector<int>  global_rays;
    
    ~ProcStatus(); 

    ProcStatus(int proc_rank, int proc_size, int cSize, int dev, bool); 

    void thread_reset(); 
    
    void proc_reset(); 

    bool all_proc_idle();

    int get_idle_proc(); 

    bool all_rays_received(); 
    
    bool update_global_rays(int* a); 

    //update chunk map , return if this proc need load new chunk
    bool update_chunk(int* s); 
    
    bool update_chunk(int candidate_proc, int candidate_chunk); 

    bool all_thread_waiting(); 
    
    void updata_local_chunk(); 
    
    bool is_local_chunk(int c); 

    void chunk_loaded() { 
        thread_reset(); 
        set_proc_busy(proc_rank); 
        load_new_chunk = false; 
    }

    bool has_new_chunk(){return load_new_chunk;}

    int* get_chunk_proc(){return &chunk_proc[0];}
    
    int get_chunk_size(){return chunk_size;}

    int get_current_chunk(){ return current_chunk;}
 
    bool generate_rays(){ return !(rough_trace && current_chunk == chunk_size - 1);}

    int get_proc(int c); 

    void set_thread_idle(int id, bool wait) { thread_idle[id] = wait;}
    
    bool is_thread_idle(int id) { return thread_idle[id]; }

    void set_stream_store_capacity(int n){stream_store_capacity = n;}
    
    int get_stream_logic_capacity(){return stream_logic_capacity; }
    int get_stream_store_capacity(){return stream_store_capacity; }
    int get_out_stream_capacity(){return out_stream_capacity; }
    int get_thread_num() {return work_thread_num; }

    void accumulate_sent(int n) { global_rays[proc_rank] += n; }
    
    void accumulate_recv(int n) { global_rays[proc_rank + proc_size] += n;}
 
    int* get_status() {return global_rays.data();}
    
    bool is_proc_idle(){return proc_idle[proc_rank]; }
    
    void set_proc_busy(int i){ proc_idle[i] = false; }
    
    void set_proc_idle(int i){proc_idle[i] = true; }
    
    int get_dev_num() {return dev_num;}
    
    ImageDecomposition* get_camera() { return camera; } 
    
    void set_new_chunk(){load_new_chunk = true;}
    
    void switch_current_chunk();

    void set_exit() {exit = true;}
    
    void set_working() {exit = false;}
    
    bool get_rough_trace() {return rough_trace;}

    bool Exit() {return exit;}

    void reset(); 
};

ProcStatus::ProcStatus(int proc_rank, int proc_size, int cSize, int dev, bool rough_trace) 
    : proc_rank(proc_rank), proc_size(proc_size), dev_num(dev), rough_trace(rough_trace)
{
    
    work_thread_num = 8;//std::thread::hardware_concurrency();
    thread_idle.resize(work_thread_num);
    thread_reset(); 
    proc_reset();
     
    exit = false;
    //including simple mesh
    chunk_size = rough_trace ? cSize + 1 : cSize;
    printf("chunk_size %d cSize %d\n", chunk_size, cSize);
    
    stream_logic_capacity = 256 * 256;//1024 * 1024 / work_thread_num;
    stream_store_capacity = (stream_logic_capacity & ~((1 << 5) - 1)) + 32; // round to 32
    
    out_stream_capacity = 1024*1024;//stream_logic_capacity; 
    
    global_rays.resize(proc_size * proc_size);
}

void ProcStatus::reset() {
    thread_idle.resize(work_thread_num);
    thread_reset(); 
    proc_reset();
    exit = false;
}

ProcStatus::~ProcStatus() {
    proc_idle.clear();
    global_rays.clear();

    printf("procstatus delete\n");
}

void ProcStatus::switch_current_chunk() {
    int local_size = local_chunks.size();
    for(int i = 0; i < local_size; i++) {
        if(local_chunks[i] == current_chunk) {
            printf("load new chunk old chunk %d || ", current_chunk);
            current_chunk = local_chunks[(i + 1) % local_size];
            load_new_chunk = true;
            printf("new chunk %d\n", current_chunk);
            return;
        }    
    }
    error("load invalid chunk");
}

void ProcStatus::updata_local_chunk() {

    current_chunk = -1;
    int i = 0;
    while(chunk_proc[i * 2] != -1) {
        printf("update chunk proc %d %d \n", chunk_proc[i * 2], chunk_proc[i * 2 + 1]);
        if(chunk_proc[i * 2 + 1] == proc_rank) {
            current_chunk = chunk_proc[i * 2];  
            local_chunks.emplace_back(chunk_proc[i * 2]);
            load_new_chunk = true;
        }
        i++;
    }
    assert(load_new_chunk);
}

bool ProcStatus::all_thread_waiting() {
    for(int i = 0; i < work_thread_num; i++) {
        if(!thread_idle[i]) 
            return false;
    }
    return true;
}

void ProcStatus::thread_reset() {
    printf("\ncpu thread num %d\n", work_thread_num);
    for(int i = 0; i < work_thread_num; i++) 
        thread_idle[i] = false;
}

void ProcStatus::proc_reset() {
    printf("proc reset %d proc idle", proc_size);
    proc_idle.resize(proc_size);
    for(int i = 0; i < proc_size; i++) 
        proc_idle[i] = false;
//    std::cout<< proc_idle[0] <<"  "<< proc_idle[1] <<std::endl;
}
 
bool ProcStatus::all_proc_idle() {
    for(int i = 0; i < proc_size; i++){
        if(!proc_idle[i]) 
            return false;
    } 
    return true;
}
    
int ProcStatus::get_proc(int c) {
    int i = 0;
    while(chunk_proc[i * 2] != -1) {
        //printf("get proc %d %d\n", );
        if(chunk_proc[i * 2] == c) {
            return chunk_proc[i * 2 + 1];
        }
        i++;
    }
    return -1;
}

bool ProcStatus::is_local_chunk(int c) {
    int local_size = local_chunks.size();
    for(int i = 0; i < local_size; i++) {
        if(local_chunks[i] == c) {
            return true;
        }    
    }
    return false;
}


int ProcStatus::get_idle_proc() {
    for(int i = 0; i < proc_size; i++){
        if(proc_idle[i]) 
            return i;
    } 
    return -1;
}

bool ProcStatus::all_rays_received() {
    int sent = 0, recv = 0;
    for(int i = 0; i < proc_size; i++) {
        sent += global_rays[i];
        recv += global_rays[i + proc_size];
    }
    printf("check all rays recv %d %d\n", sent, recv);
    if (sent <= recv )
        return true;
    else 
        return false;
}

bool ProcStatus::update_global_rays(int* a) {
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

bool ProcStatus::update_chunk(int* schedule) {
    error("update chunk 1");
    for(int i = 0; i < chunk_size; i++) {
        printf("chunk_size %d chunk_proc[i] %d  schedule[i] %d\n", chunk_size, chunk_proc[i], schedule[i]);
        if(chunk_proc[i] != schedule[i]) {
            chunk_proc[i] = schedule[i];
            if(chunk_proc[i] == proc_rank) {
                current_chunk = i;
                load_new_chunk = true;
                printf("update chunk true\n");
            }  
        } 
    } 
    printf("update chunk false\n");
    return load_new_chunk;
}

bool ProcStatus::update_chunk(int candidate_proc, int candidate_chunk) {
    //error("update chunk 2");
    int i = 0;
    while(chunk_proc[i * 2] != -1) {
        printf("update chunk proc %d %d \n", chunk_proc[i * 2], chunk_proc[i * 2 + 1]);
        if(chunk_proc[i * 2 + 1] == candidate_proc) 
            chunk_proc[i * 2 + 1] = -1;  
        if(chunk_proc[i * 2] == candidate_chunk) 
            chunk_proc[i * 2 + 1] = candidate_proc;  
        i++;
    }
    for(int i = 0; i < chunk_size; i++)
        printf("%d %d |", chunk_proc[i*2], chunk_proc[i*2+1]);
    printf("\n");
    printf("candidate chunk %d canidate proc %d\n", candidate_chunk, candidate_proc);
    if(candidate_proc == proc_rank) {
         current_chunk = candidate_chunk;
         load_new_chunk = true;
         return true;
    } 
    return false;
}

