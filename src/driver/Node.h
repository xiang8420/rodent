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
#include "decomposition.h"
#include "ProcessStatus.h"


class Node {

protected:
    Communicator *comm;
    ProcStatus *ps;
     
    int master, worker_size; 
    bool proc_idle;
    
    int recv_loop_count, master_loop_count;

    std::condition_variable buffer_not_full; // primary, secondary buffer size < max
    std::condition_variable buffer_not_empty; // primary + secondary > 0 

public: 
    Node(struct Communicator *comm, struct ProcStatus *ps)
        : comm(comm), ps(ps) 
    {
        master = comm->size - 1;
        worker_size = master;
        
        recv_loop_count = 0;
        master_loop_count = 0;    
        proc_idle = false;
    }

    ~Node(){}
    
    ProcStatus *proc_status(){return ps;}

    virtual void write_rays_buffer() = 0;

    virtual int load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) = 0;
    
    virtual void save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary) = 0;

    virtual void run(float* rProcessTime) = 0;
};




