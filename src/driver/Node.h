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
    struct RayList **rayList;
    struct RayList * outList;
    struct RayList * inList;
    struct RayList * buffer;
    Communicator *comm;
    ProcStatus *ps;
     
    int master, worker_size; 
    bool proc_idle;
    
    int recv_loop_count, master_loop_count;
    int renderer_get_rays, renderer_save_rays; 

    std::mutex  out_mutex, buffer_mutex, in_mutex, thread_mutex;
    
    std::condition_variable buffer_not_full, buffer_not_empty; //
    std::condition_variable render_start; //

public: 
    Node(Communicator *comm, ProcStatus *ps);

    ~Node() {}

    ProcStatus *proc_status(){return ps;}

    void write_rays_buffer();

    int load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait);
    
    static void work_thread(void* tmp, float* processTime, int devId, int devNum, bool preprocess, bool shoot_rays);
   
    virtual void save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary) = 0;

    virtual void run(float* rProcessTime) = 0;

    bool rayList_empty();
    
    bool inList_empty(); 
    
    bool outList_empty();
    
    bool inout_list_empty();
    
};

bool Node::outList_empty() {  
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++) {
        if(rayList[i]->type == "out" && !rayList[i] -> empty() ) {
            return false;
        }
    }
    return true;
}

bool Node::inList_empty(){ 
    return inList->empty() && buffer->empty(); 
}

bool Node::rayList_empty() {  
    int chunk_size = ps->get_chunk_size();
    bool res = true;
    out_mutex.lock();
    for(int i = 0;i < chunk_size; i++)
        if(!rayList[i]->empty()) { res = false; break; }
    out_mutex.unlock();
    in_mutex.lock();
    res &= inList->empty() && buffer->empty();
    in_mutex.unlock();
    return res;
}

bool Node::inout_list_empty() {
    return inList->empty() && outList->empty() && buffer->empty();
}

Node::Node(Communicator *comm, ProcStatus *ps) : comm(comm), ps(ps) {
    master = comm->size - 1;
    worker_size = master;
    
    recv_loop_count = 0;
    master_loop_count = 0;    
    proc_idle = false;
    
    renderer_get_rays = 0;
    renderer_save_rays = 0;
}

void Node::write_rays_buffer() {
    if(inList->empty()) return;
    
    std::unique_lock <std::mutex> lock(buffer->mutex); 
    while(!buffer->empty()) {
        buffer_not_full.wait(lock);
    }
    comm->os<<"mthread get inlist lock\n";
    
    inList->write_to_device_buffer(buffer, comm->rank);
    
    buffer_not_empty.notify_one();
    comm->os<<"mthread release inlist lock\n";
    lock.unlock();
}

int Node::load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) {
    if(comm->size == 1 || ps->Exit() || ps->has_new_chunk()) {
        comm->os << "rthread exit new chunk \n";
        return -1;
    }
    
    comm->os<<"rthread buffer size" <<buffer->size()<<"\n";
    if(buffer->empty()) { 
        comm->os<<"rthread notify" <<buffer->size()<<"\n";
        buffer_not_full.notify_all();
        comm->os<<"rthread notify" <<buffer->size()<<"\n";
    }

    struct Rays *queue = primary ? buffer->get_primary() : buffer->get_secondary();
    comm->os <<"rthread "<<thread_id<<"read incoming buffer"<<thread_wait<< "size "<<queue->size<<"\n";
    int width = primary ? 21 : 14;
    std::unique_lock <std::mutex> lock(buffer->mutex); 
    
    ps->set_thread_idle(thread_id, thread_wait);
    comm->os <<"rthread idle "<< ps->is_thread_idle(thread_id) 
             <<" width "<<width
             <<" queue->size"<< queue->size
             <<" exit"<<ps->Exit() 
             <<" buffersize"<<buffer->size()
             <<" thread wait"<<thread_wait
             << std::endl;
    while (thread_wait && buffer->empty() && !ps->Exit() && !ps->has_new_chunk()) {
        comm->os<<"rthread wait for incoming lock"<<ps->is_thread_idle(thread_id)<<" width "<<width<<"\n";
        buffer_not_empty.wait(lock);
        comm->os<<"rthread get condition" <<thread_wait<<" "<<buffer->empty()<<ps->Exit()<<"\n";
    }
    if(ps->Exit() || ps->has_new_chunk()) {  
        return -1;
    }
    comm->os<<"rthread after wait\n";
    if(!queue->empty()) {
        ps->set_thread_idle(thread_id, false);
        int copy_size = queue->size; 
        memcpy(*rays, queue->get_data(), ps->get_buffer_capacity() * width * sizeof(float)); 
        queue->size = 0;

        //printf("%d queue size %d %d %d\n", comm->rank, queue->size, primary, queue->store_width);
        int* ids = (int*)(*rays);
        comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
        for(int i = 0; i < 10; i ++) {
            comm->os<<"#" <<ids[i + rays_size]<<" "<<ids[i + rays_size + ps->get_buffer_capacity() * 9];
        }
        comm->os<<"\n";
        renderer_get_rays += copy_size;
        return copy_size + rays_size;
    }
//    if(buffer->empty()) { 
//        comm->os<<"rthread notify" <<buffer->size()<<"\n";
//        buffer_not_full.notify_all();
//        comm->os<<"rthread notify" <<buffer->size()<<"\n";
//        
//    }
    lock.unlock(); 

    comm->os <<"rthread idle "<< ps->is_thread_idle(thread_id) 
             <<" width "<<width
             <<" queue->size "<< queue->size
             <<" exit "<<ps->Exit() 
             <<" buffersize "<<buffer->size()
             <<" new chunk "<<ps->has_new_chunk()
             <<" thread wait "<<thread_wait
             << std::endl;
    if(ps->Exit() || ps->has_new_chunk()) {  
        return -1;
    }
    return rays_size;
}

void Node::work_thread(void* tmp, float *process_time, int devId, int devNum, bool preRendering, bool generate_rays) {
    printf("work threadstart\n");

    Node * wk = (Node*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus *ps = wk->ps;
    
    int region[4]; 
    int sppProc = ps->get_tile_info(&region[0], process_time); 
    int sppDev = sppProc / devNum;
    int seed = comm->rank * devNum + devId;
    printf("width %d height%d spp %d dev id %d local chunk %d\n", ps->width, ps->height, sppDev, devId, ps->get_local_chunk() );
    printf("region %d %d %d %d\n", region[0], region[1], region[2], region[3]);
    Settings settings {
        Vec3 { ps->eye.x, ps->eye.y, ps->eye.z },
        Vec3 { ps->dir.x, ps->dir.y, ps->dir.z },
        Vec3 { ps->up.x, ps->up.y, ps->up.z },
        Vec3 { ps->right.x, ps->right.y, ps->right.z },
        ps->w, ps->h,
        Vec4_i32 { region[0], region[1], region[2], region[3]},
        sppDev
    };
    if(preRendering) {
        prerender(&settings);
    } else {
        render(&settings, seed, devId, ps->get_local_chunk(), generate_rays);
    }
}

// get chunk retun chunk id
int get_sent_list(RayList ** raylist, ProcStatus *ps) {
    int n = ps->get_chunk_size();
    int t = -1; 
    int max = 0; 
    for(int i = 0; i < n; i++) {
        printf(" get sent list %d %d  %d\n", i, raylist[i] -> size(), ps->get_proc(i));
        if(raylist[i] -> size() > max && ps->get_proc(i) != -1 && i != ps->get_local_chunk()) {
            max = raylist[i] -> size();
            t = i;
        }
    }
    printf(" %d %d\n", t, max);
    if(max <= 0 || t < 0) 
        return -1;
    
    return t;
}
