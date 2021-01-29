#ifndef DEF_PROC_STATUS
#define DEF_PROC_STATUS

#include "../driver/interface.h"
#include <iostream>
#include <thread>

class ProcStatus {
private:
    //proc
    bool exit;
    int proc_size, proc_rank;
    std::vector<bool> proc_idle;

    //thread 
    std::vector<bool> thread_idle;
    int work_thread_num;
    int dev_num;

    //global capacity settings
    int stream_logic_capacity;
    int stream_store_capacity;
    int out_stream_capacity; 

    bool simple_trace;
    int chunk_size;
public:    
    ChunkManager * chunk_manager;

    std::mutex thread_mutex;
    std::vector<int>  global_rays;
    
    ~ProcStatus(); 

    ProcStatus(int proc_rank, int proc_size, int cSize, int dev); 

    void thread_reset(); 
    void proc_reset(); 
    bool all_proc_idle();
    int get_idle_proc(); 
    bool all_rays_received(); 
    bool update_global_rays(int* a); 
    bool all_thread_waiting(); 

    void set_thread_idle(int id, bool wait) { thread_idle[id] = wait;}
    bool is_thread_idle(int id) { return thread_idle[id]; }
    int get_stream_logic_capacity(){return stream_logic_capacity; }
    int get_stream_store_capacity(){return stream_store_capacity; }
    int get_out_stream_capacity(){return out_stream_capacity; }
    int get_thread_num() {return work_thread_num; }
    void accumulate_sent(int n) { global_rays[proc_rank] += n; }
    void accumulate_recv(int n) { global_rays[proc_rank + proc_size] += n;}
    int* get_status() {return global_rays.data();}
    bool is_proc_idle(){return proc_idle[proc_rank]; }
    void set_proc_idle(int i){proc_idle[i] = true; }
    void set_proc_idle(){proc_idle[proc_rank] = true; }
    int get_dev_num() {return dev_num;}
    int  get_chunk_size(){ return chunk_size; }
    void set_exit() {exit = true;}
    void set_working() {exit = false;}
    bool get_simple_trace() {return simple_trace;}
    bool Exit() {return exit;}
    void reset(); 
    
    void set_proc_busy(int i){ proc_idle[i] = false; }
    void chunk_loaded() { 
        thread_reset(); 
        set_proc_busy(proc_rank); 
        chunk_manager->local_chunks.set_new_loaded(false); 
    }
    int  switch_current_chunk(RayStreamList *outlist) { chunk_manager->switch_current_chunk(outlist); } 
    int  get_next_chunk() { return chunk_manager->get_next_chunk();  }
    bool has_new_chunk(){ return chunk_manager->local_chunks.get_new_loaded(); }
    int  get_current_chunk(){ return chunk_manager->get_current_chunk(); }
    bool generate_rays(){ return !(simple_trace && chunk_manager->get_current_chunk() == chunk_manager->chunk_size - 1); }
    int  get_proc(int c) { return chunk_manager->get_proc(c); } 
    bool update_chunk(int* schedule) { chunk_manager->update_chunk(schedule, proc_rank); } 
    bool update_chunk(int cand_proc, int cand_chunk){ chunk_manager->update_chunk(cand_proc, cand_chunk, proc_rank); } 
    bool is_local_chunk(int c) { return chunk_manager->is_local_chunk(c); }
};

ProcStatus::ProcStatus(int proc_rank, int proc_size, int csize, int dev) 
    : proc_rank(proc_rank), proc_size(proc_size), chunk_size(csize), dev_num(dev)
{
    work_thread_num = 16;//std::thread::hardware_concurrency();
    thread_idle.resize(work_thread_num);
    thread_reset();
    proc_reset();
     
    exit = false;

    //including simple mesh
    
    stream_logic_capacity = STREAM_CAPACITY;
    stream_store_capacity = (stream_logic_capacity & ~((1 << 5) - 1)) + 32; // round to 32
    out_stream_capacity = OUT_BUFFER_CAPACITY;
    global_rays.resize(proc_size * proc_size);

    simple_trace = SIMPLE_TRACE;
    chunk_manager = NULL;//new ChunkManager(simple_trace ? cSize + 1 : cSize);
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
    thread_idle.clear();
    printf("%d procstatus delete\n", proc_rank);
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
    if (sent <= recv )
        return true;
    else 
        return false;
    return true;
}

#endif
