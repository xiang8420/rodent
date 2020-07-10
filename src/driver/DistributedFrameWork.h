#include "interface.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <string.h>
#include <sys/file.h>
#include <dirent.h>
#include <unistd.h>

#include "Node.h"
#include "P2PNode.h"
#include "MasterWorker.h"

//
struct DistributedFrameWork {
    Node *node;
    std::string type;
    
    Communicator *comm;
    ProcStatus *ps;

    DistributedFrameWork(std::string type, Communicator *comm, ProcStatus *ps)
         :type(type), ps(ps), comm(comm) 
    {
        if(type == "P2P") {
            node = new P2PNode(comm, ps);
        } else if (type == "MasterWorker") {
            if(comm->rank == comm->size - 1)
                node = new Master(comm, ps);
            else    
                node = new Worker(comm, ps);
        } 
    }

    ~DistributedFrameWork() {
        delete node;
    }
    
    void run(float* time) {
        if(type == "P2P") {
            node->run(time);
        } else if (type == "MasterWorker") {
            if(comm->rank == comm->size - 1)
                node->run(time);
            else    
                node->run(time);
        }
    }
};

static std::unique_ptr<DistributedFrameWork> dfw;

void setup_distributed_framework(std::string type, struct Communicator *comm, struct ProcStatus *ps) {
    dfw.reset(new DistributedFrameWork(type, comm, ps));
}

void cleanup_distributed_framework() {
    dfw.reset();
}

void dfw_run(float* processTime) {
    dfw->run(processTime);
}

void send_rays(float *rays, size_t size, size_t capacity, bool isPrimary){
    printf("worker send\n");
    dfw->node->save_outgoing_buffer(rays, size, capacity, isPrimary);
}

int recv_rays(float **rays, size_t size, bool isPrimary, int thread_id, bool thread_wait){
    return dfw->node->load_incoming_buffer(rays, size, isPrimary, thread_id, thread_wait);
}

void master_save_ray_batches(float *rays, size_t size, size_t capacity, size_t thread_id) {
//    dfw->worker->save_ray_batches(rays, size, capacity, thread_id);
}

int32_t worker_buffer_size() {
    return dfw->node->proc_status()->get_buffer_size();
}
