#include "interface.h"
#include <thread>
#include <condition_variable>
#include <string.h>
#include <sys/file.h>
#include <dirent.h>
#include <unistd.h>
#define PRIMARY_WIDTH 21
#define SECONDARY_WIDTH 14
#define MAX_CLIENT_CHUNK 3

#include "raylist.h"
#include "process_status.h"

struct Master {
    struct Communicator *comm;
    struct RayList **raylists;
    struct RayList *buffer;
    int worker_local_chunks[128][3];
    int recv_capacity, buffer_capacity, comm_count;
    int chunk_size, worker_size; 
    bool *mpi_worker_wait;
    struct ProcStatus *rs;
    double recv_t, wait_t, send_t, total_t, read_t, write_t, st, ed;

    Master(struct Communicator *comm, struct ProcStatus *rs);
    
    ~Master();
    
    int get_dst_worker_id();

    int get_max_rays_chunk(bool unloaded); 

    bool all_queue_empty(){ 
        for(int i = 0; i < chunk_size; i++) {
            if(!raylists[i]->empty()) {return false;}
        }
        return true;
    }

    bool all_worker_finish(); 
    // ray queue empty and recv all worker end msg
    bool shutdown(){ return false; }

    void asyn_master();

    void rayQueuing(float * buffer, int size, int capacity); 
    
    void save_ray_batches(float *buffer, size_t size, size_t capacity, size_t thread_id); 
    
    void ray_batching_master(int sppTask, int film_width, int film_height); 

    void run(); 
};

class Worker {

protected:
    struct Communicator *comm;
    struct ProcStatus *rs;
    
    int buffer_capacity, buffer_size;
    int master, worker_size; 
    bool proc_idle;
    
    int recv_loop_count, master_loop_count;

    std::condition_variable buffer_not_full; // primary, secondary buffer size < max
    std::condition_variable buffer_not_empty; // primary + secondary > 0 

public: 
    Worker(struct Communicator *comm, struct ProcStatus *rs);

    ~Worker(){}
    
    virtual void write_rays_buffer() = 0;

    virtual int worker_load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) = 0;
    
    virtual void run(float* rProcessTime, int deviceNum) = 0;

    virtual void worker_save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary) = 0;
  
};

class StupidWorker : public Worker {
    
protected:    
    std::vector<int> local_chunks;
    struct RayList * inList;
    struct RayList * buffer;
    struct RayList *outList;
    bool current_chunk_empty;
public:    

    StupidWorker(struct Communicator *comm, struct ProcStatus *rs);

    ~StupidWorker();

    bool incoming_empty(){ return inList->empty() && buffer->empty(); }

    bool all_queue_empty(){ return inList->empty() && outList->empty() && buffer->empty();}

    // send primary and secondary
    void worker_save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary);
  
    void write_rays_buffer();
    
    int worker_load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait);
    
    void mpi_thread();

    static void message_thread(void* tmp); 
    
    static void work_thread(struct ProcStatus *rs, int region[4], int sppTask, int iter, int dev, int chunk, bool valid_camera);

    void run(float* rProcessTime, int deviceNum); 
}; 


class SmartWorker : public Worker {

protected:    
    int * statistic;
    struct RayList **List;
    struct RayList **outList;
    struct RayList * inList;
    struct RayList * buffer;
    
    MessageQ *inMsgQ;
    MessageQ *outMsgQ;

    std::mutex   out_mutex, buffer_mutex, in_mutex;
    
public:    
    SmartWorker(struct Communicator *comm, struct ProcStatus *rs);

    ~SmartWorker();

    bool incoming_empty(){ return inList->empty() && buffer->empty();}

    bool all_queue_empty() { 
        bool res = true;
        out_mutex.lock();
        for(int i = 0;i < worker_size; i++)
            if(!outList[i]->empty()) { res = false; break; }

        out_mutex.unlock();
        in_mutex.lock();
        res &= inList->empty() && buffer->empty();
        in_mutex.unlock();
        return res;
    }
    
    bool check_rendering_status(); 
   
    void set_distributed_buffer(); 
    // send primary and secondary
    void worker_save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary);
  
    void write_rays_buffer();
    
    int worker_load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait);
   
    void mpi_thread();
    
    void count_rays();

    static void message_thread(void* tmp);
    
    static void work_thread(void* tmp, float* processTime, int devId, int devNum, bool preprocess);

    void run(float* processTime, int deviceNum);

}; 

