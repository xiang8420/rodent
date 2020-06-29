#include "interface.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <string.h>
#include <sys/file.h>
#include <dirent.h>
#include <unistd.h>

#include "decomposition.h"
#include "communicator.h"
#include "multiproc.h"

int32_t * get_prerender_result();

Master::Master(struct Communicator *comm, struct ProcStatus *rs)
    :comm(comm), rs(rs)
{
    int chunk_size = rs->get_chunk_size();
    worker_size    = comm->size - 1;
    comm_count     = 0;
    raylists       = new RayList *[chunk_size];
    comm->os<< "master chunk size"<<chunk_size<<"\n";
    mpi_worker_wait = new bool[worker_size];
    for(int i = 0; i < worker_size; i++) {
        mpi_worker_wait[i] = false;
    } 
    
    recv_capacity = 1048608;
    buffer_capacity = 8 * recv_capacity;
    for(int i = 0; i < chunk_size; i++) {
        raylists[i] = new RayList(buffer_capacity, "out");
        worker_local_chunks[i][0] = i;
    }
    buffer = new RayList(recv_capacity, "buffer"); 
    
    st = clock();
    recv_t = 0; wait_t = 0;
    send_t = 0; total_t = 0;
    write_t = 0; read_t = 0;
}

Master::~Master(){
    int chunk_size = rs->get_chunk_size();
    for(int i = 0; i < chunk_size; i++){
        delete raylists[i];
    }
    delete buffer;
}; 

int Master::get_dst_worker_id(){
    int dst = comm_count++ % worker_size; 
    return dst; 
}

int Master::get_max_rays_chunk(bool unloaded) {
    int chunk_size = rs->get_chunk_size();
    int max_ray = 0, id_max = 0;
    comm->os<<"get max rays:";
    for(int i = 0; i < chunk_size; i++) {
        int ray_size = raylists[i]->size(); 
        comm->os<<"| "<<ray_size;
        if(ray_size > max_ray) {
            if(unloaded) {
                int j;
                for(j = 0; j < worker_size; j++) {
                    if(worker_local_chunks[j][0] == i) {
                        break;
                    }
                }
                if (j == worker_size) {
                    max_ray = ray_size;
                    id_max = i; 
                } 
            } else {
                max_ray = ray_size;
                id_max = i; 
            }
            
        }
    }
    comm->os<<"\n";
    if(max_ray == 0) return -1;
    return id_max;
}


bool Master::all_worker_finish() {
    bool finish = true;
    for(int i = 0; i < worker_size; i++){
        if(!mpi_worker_wait[i]) {
            finish = false; break;
        }
    }
    return finish;
}

void Master::asyn_master() {
    comm->os<<"asyn master\n";
    int chunk_size = rs->get_chunk_size();
    int msg[MSG_SIZE];
    int stopped = 0; 
//        return;
    for(;;) {
        int workerId = get_dst_worker_id();
        int workerChunk = worker_local_chunks[workerId][0]; 
        comm->os<<"master "<<workerId<<" "<<workerChunk<<"\n"; 
        comm->recv_msg(workerId, &msg[0]);
        comm->os<<"master recv msg"<< msg[0]<<" "<< msg[1]<<" "<<msg[2]<<" "<<msg[3]<<" "<<msg[4]<<"\n";
        mpi_worker_wait[workerId] = (msg[0] == -1);
       
        // scene chunk schedule
        if(mpi_worker_wait[workerId] && raylists[workerChunk]->empty()) {
            comm->os<<"has queues ";
            for(int i = 0; i < chunk_size;i++) {
                comm->os<<i<<"| "<<raylists[i]->primary->size<<" "<<raylists[i]->secondary->size;
            }
            comm->os<<"\n";
            if(!all_queue_empty()) {
                int tmp = workerChunk;
                workerChunk = get_max_rays_chunk(true);
                if(workerChunk != -1) {
                    msg[0] = workerChunk; 
                    msg[3] = std::min(msg[3], raylists[workerChunk]->primary->size); 
                    msg[4] = std::min(msg[4], raylists[workerChunk]->secondary->size);
                    comm->os<<"msg[3]"<<msg[3]<<" raylist primary"<<raylists[workerChunk]->primary->size
                            <<"msg[4]"<<msg[4]<<"  secondary "    <<raylists[workerChunk]->secondary->size<<"\n";
                    worker_local_chunks[workerId][0] = workerChunk;
                    comm->send_msg(workerId, &msg[0]);
                    
                    comm->recv_rays(workerId, msg[1], buffer->get_primary());
                    RayList::classification(raylists, buffer->get_primary());
                    comm->recv_rays(workerId, msg[2], buffer->get_secondary());
                    RayList::classification(raylists, buffer->get_secondary());
                    comm->os<<"master send"<<workerId<<"need new chunk"<<workerChunk;
                    mpi_worker_wait[workerId] = false; 
                    continue;     
                } else {
                    comm->os<<"no chunk need to be loaded, use old chunk \n";
                    workerChunk = tmp; 
                    // Nothing happened, the worker still idle
                }
            } else if(all_worker_finish()) {
                comm->os<<"send end"<<workerId;
                comm->send_end(workerId);
                stopped ++;
                comm->os<<"master stop"<<stopped<<"\n";
                if(stopped == worker_size) {break;}
                continue;
            } 
        } 
        comm->os<<"msg[3]"<<msg[3]<<" raylist primary"<<raylists[workerChunk]->primary->size
                <<"msg[4]"<<msg[4]<<"  secondary "    <<raylists[workerChunk]->secondary->size<<"\n";
        
        msg[0] = workerChunk;
        msg[3] = std::min(msg[3], raylists[workerChunk]->primary->size); 
        msg[4] = std::min(msg[4], raylists[workerChunk]->secondary->size);
        
        comm->send_msg(workerId, &msg[0]);
        //if workerId idle send new rays and new scene id 
        comm->recv_rays(workerId, msg[1], buffer->get_primary());
        RayList::classification(raylists, buffer->get_primary());
        comm->os<<"master classify primary\n";
        comm->recv_rays(workerId, msg[2], buffer->get_secondary());
        RayList::classification(raylists, buffer->get_secondary());
        
        comm->os<<"master send msg"<<msg[0]<<" "<<msg[1]<<" "<<msg[2]<<" "<<msg[3]<<" "<<msg[4];
        comm->send_rays(workerId, msg[3], raylists[workerChunk]->primary);
        comm->send_rays(workerId, msg[4], raylists[workerChunk]->secondary);
        if(msg[3] > 1 || msg[4] > 0) {
            mpi_worker_wait[workerId] = false;
        } 
    }
    comm->os<<"master send all over\n";
    for(int i = 0; i < worker_size;i++){
        printf("%d %d \n", i, raylists[i]->primary->get_size());
    }
    printf("\n");
    for(int i = 0; i < worker_size;i++){
        printf("%d %d \n", i, raylists[i]->secondary->get_size());
    }
    printf("\n");
    total_t = (double)(clock() - st) / CLOCKS_PER_SEC;
    printf( "\n Master recv time%f send time %f  total %f \n write %f read %fseconds\n", recv_t, send_t, total_t, write_t, read_t );
}

void Master::rayQueuing(float * buffer, int size, int capacity) {
    int *grid  = (int *)&buffer[9 * capacity];
    int st = 0, ed = 0;
    while(st < size) {
        ed = st;
        int st_gId = grid[st];
        int ed_gId = grid[ed];
        while(st_gId == ed_gId && ed < size) {
            ed ++;
            ed_gId = grid[ed];
        }
        printf("|%d %d %d|", st_gId, st, ed);
        raylists[st_gId]->primary->read_dev_rays(buffer, st, ed - st, capacity);
        st = ed;
    } 
//        printf("\n");
}

void Master::save_ray_batches(float *buffer, size_t size, size_t capacity, size_t thread_id) {
    rayQueuing(buffer, size, capacity);   
    int chunk_size = rs->get_chunk_size();
    printf("%ld size %ld batch status \n", thread_id, size);
    for(int i = 0; i< chunk_size;i++) {
        printf("%d ", raylists[i]->primary->size); 
    }
    printf("\n");
}

void Master::ray_batching_master(int sppTask, int film_width, int film_height) {
    Settings settings {
        Vec3 { rs->eye.x, rs->eye.y, rs->eye.z },
        Vec3 { rs->dir.x, rs->dir.y, rs->dir.z },
        Vec3 { rs->up.x, rs->up.y, rs->up.z },
        Vec3 { rs->right.x, rs->right.y, rs->right.z },
        rs->w, rs->h,
        Vec4_i32 { 0, 0, film_width, film_height},
        sppTask
    };
    raybatching(&settings, 0, 0);
    
    int msg[MSG_SIZE];
    int stopped = 0; 
    for(int i = 0; i < worker_size; i++) {
        worker_local_chunks[i][0] = -2;
    }
    
    for(;;) { 
        printf("master wait for msg\n");
        int workerId = get_dst_worker_id();
        
        //comm->Gather_worker_info();
        comm->recv_msg(workerId, &msg[0]);
        printf("master recv msg %d %d %d %d %d\n", msg[0], msg[1], msg[2], msg[3], msg[4]);
        
        int workerChunk = worker_local_chunks[workerId][0]; 
        mpi_worker_wait[workerId] = (msg[0] == -1);
        if(workerChunk < 0 || (mpi_worker_wait[workerId] && raylists[workerChunk]->empty())) {
            int chunk_size = rs->get_chunk_size();
            for(int i = 0; i < chunk_size;i++) {
                printf("%d %d %d|", i, raylists[i]->primary->size, raylists[i]->secondary->size);
            }
            printf("\n");
            if(!all_queue_empty()) {
                int tmp = workerChunk;
                workerChunk = get_max_rays_chunk(true);
                if(workerChunk != -1) {
                    msg[0] = workerChunk; 
                    worker_local_chunks[workerId][0] = workerChunk;
                    comm->send_msg(workerId, &msg[0]);
                    comm->recv_rays(workerId, msg[1], buffer->get_primary());
                    RayList::classification(raylists, buffer->get_primary());
                    comm->recv_rays(workerId, msg[2], buffer->get_secondary());
                    RayList::classification(raylists, buffer->get_secondary());
                    printf("master send %d need new chunk%d", workerId, workerChunk);
                    mpi_worker_wait[workerId] = false; 
                    continue;     
                } else {
                    printf("no chunk need to be loaded, use old chunk \n");
                    workerChunk = tmp; 
                    // Nothing happened, the worker still idle
                }
            } else if(all_worker_finish()) {
                printf("send end %d", workerId);
                comm->send_end(workerId);
                stopped ++;
                printf("master stop %d\n", stopped);
                if(stopped == worker_size) {break;}
                continue;
            } 
        } 
        msg[0] = workerChunk;
        printf("master before send msg %d %d %d %d %d\n", msg[0], msg[1], msg[2], msg[3], msg[4]);
        msg[3] = std::min(msg[3], raylists[workerChunk]->primary->size); 
        msg[4] = std::min(msg[4], raylists[workerChunk]->secondary->size);
        
        comm->send_msg(workerId, &msg[0]);
        printf("master send msg %d %d %d %d %d\n", msg[0], msg[1], msg[2], msg[3], msg[4]);
        //if workerId idle send new raylists and new scene id 
        comm->recv_rays(workerId, msg[1], buffer->get_primary());
        RayList::classification(raylists, buffer->get_primary());
        comm->recv_rays(workerId, msg[2], buffer->get_secondary());
        RayList::classification(raylists, buffer->get_secondary());
        
        comm->send_rays(workerId, msg[3], raylists[workerChunk]->primary);
        comm->send_rays(workerId, msg[4], raylists[workerChunk]->secondary);
        if(msg[3] > 1 || msg[4] > 0) {
            mpi_worker_wait[workerId] = false;
        } 
        //breadcast control worker recv or send to other worker?
        //if all done return
        //send to worker
    }
    total_t = (double)(clock() - st) / CLOCKS_PER_SEC;
    printf( "\n Master recv time%f send time %f  total %f \n write %f read %fseconds\n", recv_t, send_t, total_t, write_t, read_t );
}

void Master::run() {
    printf("master->run()"); 
    if(rs->isRayQueuing()) {
        printf("batching\n");
        ray_batching_master(rs->spp, rs->width, rs->height);
    } else {
        printf("asyn master\n");
        asyn_master(); 
    }
}

Worker::Worker(struct Communicator *comm, struct ProcStatus *rs) 
        : comm(comm), rs(rs)
{
    master = comm->size - 1;
    worker_size = master;
    buffer_capacity = 1048608;
    buffer_size = 1048576;
    
    recv_loop_count = 0;
    master_loop_count = 0;    
    proc_idle = false;
    
}


//
StupidWorker::StupidWorker(struct Communicator *comm, struct ProcStatus *rs) 
    : Worker(comm, rs) {
    inList  = new RayList(buffer_capacity * 4, "in");
    outList = new RayList(buffer_capacity * 4, "out");
    buffer   = new RayList(buffer_capacity, "buffer");
    current_chunk_empty = false; 
}

StupidWorker::~StupidWorker(){
    delete outList;
    delete inList;
    delete buffer; 
}

void StupidWorker::write_rays_buffer() {
    if(inList->empty()) return;
    comm->os<<"mthread before copy primaryary  "<<inList->get_primary()->size<<"|"<<inList->get_primary()->empty() 
                     <<" buffer primary: "<<buffer->get_primary()->size<<"|"<< buffer->get_primary()->empty()<<"\n";
    
    if(buffer->get_primary()->empty() && !inList->get_primary()->empty()) {
        struct Rays *rays = buffer->get_primary();
        rays->size = inList->get_primary()->write_dev_rays(rays->get_data(), buffer_size, buffer_capacity);
    }
    if(buffer->get_secondary()->empty() && !inList->get_secondary()->empty()) {
        struct Rays *rays = buffer->get_secondary(); 
        rays->size = inList->get_secondary()->write_dev_rays(rays->get_data(), buffer_size, buffer_capacity);
    }
    if(!buffer->empty()) {
         buffer_not_empty.notify_one();
    }
}

int StupidWorker::worker_load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) {
    struct Rays *queue = primary ? buffer->get_primary() : buffer->get_secondary();
    comm->os <<"rthread "<<thread_id<<"read incoming buffer"<<thread_wait<< "size "<<queue->size<<"\n";
    int width = primary ? 21 : 14;
    std::unique_lock <std::mutex> inList_lock(inList->mutex); 
    
    rs->set_thread_idle(thread_id, thread_wait);
    comm->os <<"rthread idle "<< rs->is_thread_idle(thread_id) 
             <<" width "<<width
             <<" queue->size"<< queue->size
             <<" exit"<<rs->Exit() 
             <<" buffersize"<<buffer->size()
             <<" thread wait"<<thread_wait
             <<std::endl;
    while (thread_wait && buffer->empty() && !current_chunk_empty) {
        comm->os<<"rthread wait for incoming lock"<<rs->is_thread_idle(thread_id)<<"\n";
        buffer_not_empty.wait(inList_lock);
        comm->os<<"rthread get not empty condition" <<thread_wait<<" "<<buffer->empty()<<rs->Exit()<<"\n";
    }
    if(!queue->empty() && rays_size < buffer_size) {
        rs->set_thread_idle(thread_id, false);
        int copy_size;
        if(rays_size == 0) {
            copy_size = queue->size; 
            queue->size = 0;
            memcpy(*rays, queue->data, buffer_capacity * width * sizeof(float)); 
        } else {
            copy_size = std::min(queue->size , buffer_size - (int)rays_size); 
            queue->size = queue->size - copy_size;
            int src = queue->size; 
            int dst = rays_size;
            for(int i = 0; i < width; i++) {
               memcpy( &rays[dst + i * buffer_capacity], &queue->data[src + i * buffer_capacity], copy_size * sizeof(float));
            }
        }
        int* ids = (int*)(*rays);
        comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
        for(int i = 0; i < 10; i ++) {
            comm->os<<"#" <<ids[i + rays_size]<<" "<<ids[i + rays_size + 1048608 * 9];
        }
        comm->os<<"\n";
        return copy_size + rays_size;
    }
        
    inList_lock.unlock(); 
    if(current_chunk_empty){  
        printf("ray size %ld   queue primary  size %d queue secondary size %d\n", 
                rays_size, inList->get_primary()->size, inList->get_secondary()->size);
        printf("recv stop mpi %d thread %d %ld\n", comm->rank, thread_id, rays_size);
        return -1;
    }
    return rays_size;
}

void StupidWorker::worker_save_outgoing_buffer(float *retired_rays, size_t size, size_t capacity, bool primary){
    int *id = (int *) retired_rays;
    int *chunk = (int *) &retired_rays[capacity * 9];
    struct Rays *list = primary ? outList->get_primary() : outList->get_secondary();
    list->lock();
    list->read_dev_rays(retired_rays, 0, size, capacity);
    list->unlock();
}


void StupidWorker::mpi_thread() {
    outList->lock();
    int msg[MSG_SIZE];
    proc_idle =  rs->all_thread_waiting() && all_queue_empty();
    msg[0] = proc_idle ? -1 : 1;
    msg[1] = outList -> get_primary() -> size; 
    msg[2] = outList -> get_secondary() -> size;
    msg[3] = buffer_size - inList -> get_primary() -> size; 
    msg[4] = buffer_size - inList -> get_secondary() -> size;
    // worker send information to master
    comm->send_msg(master, &msg[0]);
    comm->os<<"worker chunk"<<local_chunks[0]<<"thread idle"<<rs->is_thread_idle(0)<<"all t wait"<<rs->all_thread_waiting()<<" "
           << "queue empty"<<all_queue_empty()<<" msg "<<msg[0]<<" "<<msg[1]<<" "<<msg[2]<<"inlist primaryi "<<inList->get_primary()->size 
           << " secondary "<<inList->get_secondary()->size<<"buffer primary "<< buffer->get_primary()->size<<"secondary"<<buffer->get_secondary()->size;
    comm->recv_msg(master, &msg[0]);
    // finish or need new chunk
    if(msg[0] != local_chunks[0] ) {
        local_chunks[0] = msg[0];
        current_chunk_empty = true;
        comm->os<<"recv stop\n";
        buffer_not_empty.notify_all();
        return;
    }
    comm->os<<comm->rank<<"recv msg"<<msg[0]<<" "<< msg[1]<<" "<<msg[2]<<" "<<msg[3];
    comm->send_rays(master, msg[1], outList->get_primary());
    comm->send_rays(master, msg[2], outList->get_secondary());
    outList->unlock();
    comm->recv_rays(master, msg[3], inList->get_primary());
    comm->recv_rays(master, msg[4], inList->get_secondary());
    inList->lock(); 
    write_rays_buffer(); 
    inList->unlock(); 
    
    master_loop_count ++;
    usleep(100);
}

void StupidWorker::message_thread(void* tmp) {
    struct StupidWorker *p = (struct StupidWorker*)tmp;
    //start recv thread and render
//        return;
    for(;;) {
        p->mpi_thread();
        if (p->current_chunk_empty){
            printf("worker comm proc return\n ");
            return;
        }
    }   
}

void StupidWorker::work_thread(struct ProcStatus *rs, int region[4], int sppTask, int iter, int dev, int chunk, bool camera_ray){
    printf("image region %d %d %d %d", region[0], region[1], region[2], region[3]);
    Settings settings {
        Vec3 { rs->eye.x, rs->eye.y, rs->eye.z },
        Vec3 { rs->dir.x, rs->dir.y, rs->dir.z },
        Vec3 { rs->up.x, rs->up.y, rs->up.z },
        Vec3 { rs->right.x, rs->right.y, rs->right.z },
        rs->w, rs->h,
        Vec4_i32 { region[0], region[1], region[2], region[3]},
        sppTask
    };
    render(&settings, iter + dev, dev, chunk, camera_ray);
}

void StupidWorker::run(float* frame_time) {
    
    //multi thread
    int deviceNum = rs->get_dev_num();
    int region[4];
    int spp = rs->spp; 
    int sppProc = rs->get_tile_info(&region[0], frame_time); 
    printf("stupid worker %d\n", comm->rank); 
    int sppDev = sppProc / deviceNum;
    int iter = 0; 
    if(rs->isRayQueuing()) {
        local_chunks.emplace_back(-2);
        iter = 1;
    } else {
        local_chunks.emplace_back(comm->rank);
    }
    while(true) {
        current_chunk_empty = false;
        int chunk = local_chunks[0];
        printf("%d worker %d\n", comm->rank, chunk);
        if( chunk == -1) {
            return;
        }
        //Device thread cpu thread num
        std::vector<std::thread> threadPool;
        for(int i = 0; i < deviceNum; i++) {
            printf("worker %d chunk %d image region %d %d %d %d", comm->rank, chunk, region[0], region[1], region[2], region[3]);
            threadPool.emplace_back(std::thread(work_thread, rs, region, sppDev, comm->rank, i, chunk, iter == 0));
            iter ++;
        }
        if( master != -1) {
            threadPool.emplace_back(std::thread(message_thread, (void*)this));
        }
        for( auto &thread: threadPool) {
            thread.join();
        }
    }
}

// Asynchronous p2p worker without master

void SmartWorker::set_distributed_buffer() {
    int chunk_size = rs->get_chunk_size();
    List = new RayList *[chunk_size];
    for(int i = 0; i < chunk_size; i++) {
        List[i] = new RayList(buffer_capacity * 4, "out");
    }
    for(auto c : rs->get_local_chunk()) {
        List[c]->type = "in";
    }
    outList = List;
    comm->os << "out list size"<<outList[0]->size() << "capacity" << outList[0]->get_primary()->get_capacity() << std::endl; 
    comm->os << "loaded chunk()"<<rs->get_loaded_chunk()<< std::endl; 
    //when we need to load a new chunk updata inList pointer.
    inList  = List[rs->get_loaded_chunk()];
    buffer  = new RayList(buffer_capacity, "buffer");

    rs->set_chunks();
}

SmartWorker::SmartWorker(struct Communicator *comm, struct ProcStatus *rs) 
    : Worker(comm, rs) 
{
    int chunk_size = rs->get_chunk_size();
    statistic = new int[chunk_size * 4/*dep*/ * 2]; 
    std::fill(statistic, statistic + chunk_size * 4/*dep*/ * 2, 0);
}

SmartWorker::~SmartWorker() {
    delete buffer;
    for(int i = 0; i < rs->get_chunk_size(); i++){
        delete List[i];
    }
}

void SmartWorker::write_rays_buffer() {
    comm->os<<"mthread write ray buffer"<<inList->empty()<<"\n";
    if(inList->empty()) return;
    comm->os<<"mthread inlist not empty\n";
    inList->lock();
    comm->os<<"mthread get inlist lock\n";
    if(buffer->get_primary()->empty() && !inList->get_primary()->empty()) {
        comm->os<<"mthread before copy primaryary  "<<inList->get_primary()->size<<"|"<<inList->get_primary()->empty() 
                         <<" buffer primary: "<<buffer->get_primary()->size<<"|"<< buffer->get_primary()->empty()<<"\n";
        struct Rays *rays = buffer->get_primary();
        rays->size = inList->get_primary()->write_dev_rays(rays->get_data(), buffer_size, buffer_capacity);
        comm->os<<"mthread after copy primaryary  "<<inList->get_primary()->size<<"|"<<inList->get_primary()->empty() 
                         <<" buffer primary: "<<buffer->get_primary()->size<<"|"<< buffer->get_primary()->empty()<<"\n";
    }
    if(buffer->get_secondary()->empty() && !inList->get_secondary()->empty()) {
        comm->os<<"mthread before copy secondaryary  "<<inList->get_secondary()->size<<"|"<<inList->get_secondary()->empty() 
                         <<" buffer secondary: "<<buffer->get_secondary()->size<<"|"<< buffer->get_secondary()->empty()<<"\n";
        comm->os<<"copy secondary\n";
        struct Rays *rays = buffer->get_secondary(); 
        rays->size = inList->get_secondary()->write_dev_rays(rays->get_data(), buffer_size, buffer_capacity);
        comm->os<<"mthread after copy secondaryary  "<<inList->get_secondary()->size<<"|"<<inList->get_secondary()->empty() 
                         <<" buffer secondary: "<<buffer->get_secondary()->size<<"|"<< buffer->get_secondary()->empty()<<"\n";
    }
    buffer_not_empty.notify_one();
    inList->unlock();
    comm->os<<"mthread release inlist lock\n";
}

int SmartWorker::worker_load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) {
    if(comm->size == 1 || rs->Exit()) {
        comm->os << "rthread exit\n";
        return -1;
    }    
    struct Rays *queue = primary ? buffer->get_primary() : buffer->get_secondary();
    comm->os <<"rthread "<<thread_id<<"read incoming buffer"<<thread_wait<< "size "<<queue->size<<"\n";
    int width = primary ? 21 : 14;
    std::unique_lock <std::mutex> inList_lock(inList->mutex); 
    
    rs->set_thread_idle(thread_id, thread_wait);
    comm->os <<"rthread idle "<< rs->is_thread_idle(thread_id) 
             <<" width "<<width
             <<" queue->size"<< queue->size
             <<" exit"<<rs->Exit() 
             <<" buffersize"<<buffer->size()
             <<" thread wait"<<thread_wait
             <<std::endl;
    while (thread_wait && buffer->empty() && !rs->Exit()) {
        comm->os<<"rthread wait for incoming lock"<<rs->is_thread_idle(thread_id)<<"\n";
        buffer_not_empty.wait(inList_lock);
        comm->os<<"rthread get not empty condition" <<thread_wait<<" "<<buffer->empty()<<rs->Exit()<<"\n";
    }

    if(!queue->empty()) {
        rs->set_thread_idle(thread_id, false);
        int copy_size = queue->size; 
        queue->size = 0;
        memcpy(*rays, queue->data, buffer_capacity * width * sizeof(float)); 

        //printf("%d queue size %d %d %d\n", comm->rank, queue->size, primary, queue->store_width);
        int* ids = (int*)(*rays);
        comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
        for(int i = 0; i < 10; i ++) {
            comm->os<<"#" <<ids[i + rays_size]<<" "<<ids[i + rays_size + 1048608 * 9];
        }
        comm->os<<"\n";
        return copy_size + rays_size;
    }
        
    inList_lock.unlock(); 

    if(rs->Exit()){  
//        printf("ray size %ld   queue primary  size %d queue secondary size %d\n", 
//                rays_size, inList->get_primary()->size, inList->get_secondary()->size);
//        printf("recv stop mpi %d thread %d %ld\n", comm->rank, thread_id, rays_size);
        return -1;
    }
    return rays_size;

}

void SmartWorker::worker_save_outgoing_buffer(float *retired_rays, size_t size, size_t capacity, bool primary){
    comm->os<<"rthread save outgoing buffer"<<size<<"\n";
    out_mutex.lock(); 
    int width = primary?21:14; 
    int* ids = (int*)(retired_rays);
    for(int i = 0; i < 5; i ++) {
        comm->os<<"| "<< ids[i] <<" "<< ids[i + 1048608 * 9] << " ";
    }
    comm->os<<"\n";
    RayList::classification(outList, retired_rays, size, capacity, primary);
    out_mutex.unlock(); 
}

void SmartWorker::work_thread(void* tmp, float *process_time, int devId, int devNum, bool preRendering) {
  
    SmartWorker * wk = (struct SmartWorker*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus *rs = wk->rs;
   
    int region[4]; 
    int sppProc = rs->get_tile_info(&region[0], process_time); 
    int sppDev = sppProc / devNum;
    printf("width %d height%d spp %d\n", rs->width, rs->height, sppDev);
    Settings settings {
        Vec3 { rs->eye.x, rs->eye.y, rs->eye.z },
        Vec3 { rs->dir.x, rs->dir.y, rs->dir.z },
        Vec3 { rs->up.x, rs->up.y, rs->up.z },
        Vec3 { rs->right.x, rs->right.y, rs->right.z },
        rs->w, rs->h,
        Vec4_i32 { region[0], region[1], region[2], region[3]},
        sppDev
    };
    if(preRendering) {
        prerender(&settings);
    } else {
        render(&settings, devId, devId, wk->rs->get_loaded_chunk(), true);
    }
}

int get_sent_list(RayList ** raylist, int n) {
    int t = -1, max = 0;
    for(int i = 0; i < n; i++) {
        if(raylist[i] -> size() > max && raylist[i]->type == "out") {
            max = raylist[i] -> size();
            t = i;
        }
    }
    if(max <= 0 || t < 0) 
        return -1;
    
    return t;
}

//check proc status, return if proc need to wait 
bool SmartWorker::check_rendering_status() {
    if(rs->all_thread_waiting()) {
        if(all_queue_empty()) {
            comm->os<< "mthread if all thread idle all queue empty, set itself idle\n";
            rs->set_self_idle();
            if(rs->all_proc_idle() && rs->all_rays_received()) {
                //if all worker end synchronous 
                comm->os<<"mthread all proc idle send collective\n";
                QuitMsg coll_msg(comm->rank); 
                comm->send_message(&coll_msg, rs);
                rs->set_exit();
                buffer_not_empty.notify_all();
            } else {
                comm->os<<"mthread proc idle send status\n";
                StatusMsg status(comm->rank, rs->get_status(), comm->size); 
                comm->send_message(&status, rs);
            }
            return true; 
        } else if(!incoming_empty()){
            comm->os<<"mthread incoming rays left inlist:"<<inList->size()<<" buffer:"<<buffer->size()<<"\n";
            rs->set_proc_busy(comm->rank);
            return false;
        } else {
            //load new chunk
            comm->os<<"mthread load another chunk\n";
            return true; 
        }
    } else {
        comm->os<<"mthread proc busy\n";
        rs->set_proc_busy(comm->rank);
        return false;
    }
} 

void SmartWorker::message_thread(void* tmp) {
  
    SmartWorker *wk = (struct SmartWorker*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus * rs = wk->rs;

    comm->os<<"message thread\n";
    int recv_count = 0;
    bool wait = false;
    while (!rs->Exit()) {
       
        wait = wk->check_rendering_status(); 
        if(rs->Exit()) break;
        
        while(wait && !rs->Exit()) {
            wait = !comm->recv_message(wk->List, rs);
        }
        
        wk->write_rays_buffer();
        
        if(rs->Exit()) break;

        int cId = get_sent_list(wk->outList, rs->get_chunk_size());
        comm->os<<"mthread get send list"<<cId<<"\n";
        
        if(cId >= 0) {
            wk->out_mutex.lock();
            comm->os<<"mthread new RayMsg"<<cId<<"\n";
            RayMsg *ray_msg = new RayMsg(wk->outList[cId], comm->rank, rs->get_dst_proc(cId), cId, false); 
            comm->os<<"mthread RayMsg"<<ray_msg->get_chunk()<<"\n";
            wk->out_mutex.unlock(); 
            comm->send_message(ray_msg, rs);
        }
        comm->purge_completed_mpi_buffers();
        comm->os<<"mthread single loop end \n";
	}
    comm->os<<"end message thread"<<rs->all_thread_waiting()<<"\n";
    wk->buffer_not_empty.notify_all();
    return;
} 

void SmartWorker::count_rays() {
    int width = rs->width; 
    int height = rs->height;
    int chunk_size = rs->get_chunk_size();
    int *data = get_prerender_result();
    int pixel_size = width * height;
    for(int dep = 0; dep < 8; dep++) {
        for(int i = 0; i < height; i ++) {
            for(int j = 0; j < width; j++) {
                int chunk = data[j + i * width + dep * pixel_size];
                if(chunk != -1)
                    statistic[dep * chunk_size + chunk] ++;
            }
        }
    }
    printf("statistic: ");
    for(int dep = 0; dep < 8; dep++) {
        for(int i = 0; i < chunk_size; i ++) {
            printf("%d ", statistic[dep * chunk_size + i]);
        }
        printf("\n");
    }
}
void SmartWorker::run(float* frame_time) {
    int deviceNum = rs->get_dev_num();
    bool PreRendering = false;
    if(PreRendering && comm->rank == 0) {
        work_thread(this, frame_time, 0, 1, true);
        count_rays(); 
    }
    // set domain and image distribution
    set_distributed_buffer();
    
    std::vector<std::thread> workThread;
    if(comm->size == 1) {
        printf("only one worker %d  ", comm->rank);
        for(int i = 0; i < deviceNum; i++) 
            workThread.emplace_back(std::thread(work_thread, this, frame_time, i, deviceNum, false));
        
        for( auto &thread: workThread) 
            thread.join();
    	
    } else {
        std::thread mthread(message_thread, this);
        //    //all proc start proc 0 send schedule? 
        while(!rs->Exit()) {
            int chunk = rs->get_new_chunk(); 
            for(int i = 0; i < deviceNum; i++) 
                workThread.emplace_back(std::thread(work_thread, this, frame_time, i, deviceNum, false));
            
            for( auto &thread: workThread) 
                thread.join();
            
        }
        printf("%d worker exit\n", comm->rank);
        mthread.join();
    }
    return;
} 


static std::unique_ptr<Master> master;
static std::unique_ptr<Worker> worker;

void setup_master(struct Communicator *comm, struct ProcStatus *rs) {
    master.reset(new Master(comm, rs));
}
void master_run() {
    master->run();
}

void cleanup_master(){
    master.reset();
}

void cleanup_worker() {
    worker.reset();
}

void setup_worker(struct Communicator *comm, struct ProcStatus *rs, bool with_master) {
    if(with_master) {
        worker.reset(new StupidWorker(comm, rs));
    } else {
        worker.reset(new SmartWorker(comm, rs));
    }
}

void worker_run(float* processTime) {
    worker->run(processTime);
}

void worker_send_rays(float *rays, size_t size, size_t capacity, bool isPrimary){
    printf("worker send\n");
    worker->worker_save_outgoing_buffer(rays, size, capacity, isPrimary);
}
int  worker_recv_rays(float **rays, size_t size, bool isPrimary, int thread_id, bool thread_wait){
    return worker->worker_load_incoming_buffer(rays, size, isPrimary, thread_id, thread_wait);
}
void master_save_ray_batches(float *rays, size_t size, size_t capacity, size_t thread_id) {
    master->save_ray_batches(rays, size, capacity, thread_id);
}


