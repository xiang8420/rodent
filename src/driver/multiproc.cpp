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

Master::Master(struct Communicator *comm, struct ProcStatus *ps)
    :comm(comm), ps(ps)
{
    int chunk_size = ps->get_chunk_size();
    worker_size    = comm->size - 1;
    comm_count     = 0;
    raylists       = new RayList *[chunk_size];
    comm->os<< "master chunk size"<<chunk_size<<"\n";
    mpi_worker_wait = new bool[worker_size];
    for(int i = 0; i < worker_size; i++) {
        mpi_worker_wait[i] = false;
    } 
    
    for(int i = 0; i < chunk_size; i++) {
        raylists[i] = new RayList(8 * ps->get_buffer_size(), "out", false);
        worker_local_chunks[i][0] = i;
    }
    buffer = new RayList(ps->get_buffer_size(), "buffer", false); 
    
    st = clock();
    recv_t = 0; wait_t = 0;
    send_t = 0; total_t = 0;
    write_t = 0; read_t = 0;
}

Master::~Master(){
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++){
        delete raylists[i];
    }
    delete buffer;
    delete comm;
}; 

int Master::get_dst_worker_id(){
    int dst = comm_count++ % worker_size; 
    return dst; 
}

int Master::get_max_rays_chunk(bool unloaded) {
    int chunk_size = ps->get_chunk_size();
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
    int chunk_size = ps->get_chunk_size();
    int msg[MSG_SIZE];
    int stopped = 0; 
//        return;
    int recv_size = 0;
    int sent_size = 0;
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
                    recv_size = recv_size + msg[1] + msg[2];
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
        recv_size = recv_size + msg[1] + msg[2];
        comm->recv_rays(workerId, msg[1], buffer->get_primary());
        RayList::classification(raylists, buffer->get_primary());
        comm->os<<"master classify primary\n";
        comm->recv_rays(workerId, msg[2], buffer->get_secondary());
        RayList::classification(raylists, buffer->get_secondary());
        
        comm->os<<"master send msg"<<msg[0]<<" "<<msg[1]<<" "<<msg[2]<<" "<<msg[3]<<" "<<msg[4];
        sent_size = sent_size + msg[3] + msg[4];
        comm->send_rays(workerId, msg[3], raylists[workerChunk]->primary);
        comm->send_rays(workerId, msg[4], raylists[workerChunk]->secondary);
        if(msg[3] > 1 || msg[4] > 0) {
            mpi_worker_wait[workerId] = false;
        } 
    }
    comm->os<<"master send all over\n";
    comm->os<<"master sent "<<sent_size<<" recv "<<recv_size<<"\n ";
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
        raylists[st_gId]->primary->read_dev_rays(buffer, st, ed - st, capacity, 0);
        st = ed;
    } 
//        printf("\n");
}

void Master::save_ray_batches(float *buffer, size_t size, size_t capacity, size_t thread_id) {
    rayQueuing(buffer, size, capacity);   
    int chunk_size = ps->get_chunk_size();
    printf("%ld size %ld batch status \n", thread_id, size);
    for(int i = 0; i< chunk_size;i++) {
        printf("%d ", raylists[i]->primary->size); 
    }
    printf("\n");
}

void Master::ray_batching_master(int sppTask, int film_width, int film_height) {
    Settings settings {
        Vec3 { ps->eye.x, ps->eye.y, ps->eye.z },
        Vec3 { ps->dir.x, ps->dir.y, ps->dir.z },
        Vec3 { ps->up.x, ps->up.y, ps->up.z },
        Vec3 { ps->right.x, ps->right.y, ps->right.z },
        ps->w, ps->h,
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
            int chunk_size = ps->get_chunk_size();
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
    if(ps->isRayQueuing()) {
        printf("batching\n");
        ray_batching_master(ps->spp, ps->width, ps->height);
    } else {
        printf("asyn master\n");
        asyn_master(); 
    }
}

Worker::Worker(struct Communicator *comm, struct ProcStatus *ps) 
        : comm(comm), ps(ps)
{
    master = comm->size - 1;
    worker_size = master;
    
    recv_loop_count = 0;
    master_loop_count = 0;    
    proc_idle = false;
}


//
StupidWorker::StupidWorker(struct Communicator *comm, struct ProcStatus *ps) 
    : Worker(comm, ps) {
    inList  = new RayList(ps->get_buffer_size() * 4, "in", false);
    outList = new RayList(ps->get_buffer_size() * 4, "out", false);
    buffer  = new RayList(ps->get_buffer_size(),     "buffer", false);
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
        rays->size = inList->get_primary()->copy_to_buffer(rays->get_data(), ps->get_buffer_size(), ps->get_buffer_capacity(), comm->rank, true);
    }
    if(buffer->get_secondary()->empty() && !inList->get_secondary()->empty()) {
        struct Rays *rays = buffer->get_secondary(); 
        rays->size = inList->get_secondary()->copy_to_buffer(rays->get_data(), ps->get_buffer_size(), ps->get_buffer_capacity(), comm->rank, false);
    }
    if(!buffer->empty()) {
         buffer_not_empty.notify_one();
    }
}

int StupidWorker::worker_load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) {
    struct Rays *queue = primary ? buffer->get_primary() : buffer->get_secondary();
    comm->os <<"rthread "<<thread_id<<"read incoming buffer"<<thread_wait<< "size "<<queue->size<<"\n";
    int width = primary ? 21 : 14;
    int buffer_capacity = ps->get_buffer_capacity();
    int buffer_size = ps->get_buffer_size();
    std::unique_lock <std::mutex> inList_lock(inList->mutex); 
    
    ps->set_thread_idle(thread_id, thread_wait);
    comm->os <<"rthread idle "<< ps->is_thread_idle(thread_id) 
             <<" width "<<width
             <<" queue->size"<< queue->size
             <<" exit"<<ps->Exit() 
             <<" buffersize"<<buffer->size()
             <<" thread wait"<<thread_wait
             <<std::endl;
    while (thread_wait && buffer->empty() && !current_chunk_empty) {
        comm->os<<"rthread wait for incoming lock"<<ps->is_thread_idle(thread_id)<<"\n";
        buffer_not_empty.wait(inList_lock);
        comm->os<<"rthread get not empty condition" <<thread_wait<<" "<<buffer->empty()<<ps->Exit()<<"\n";
    }
    if(!queue->empty() && rays_size < buffer_size) {
        ps->set_thread_idle(thread_id, false);
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
            comm->os<<"#" <<ids[i + rays_size]<<" "<<ids[i + rays_size + buffer_capacity * 9];
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
    outList->lock();
    list->read_dev_rays(retired_rays, 0, size, capacity, comm->rank);
    outList->unlock();
}


void StupidWorker::mpi_thread() {
    outList->lock();
    int msg[MSG_SIZE];
    proc_idle =  ps->all_thread_waiting() && all_queue_empty();
    msg[0] = proc_idle ? -1 : 1;
    msg[1] = outList -> get_primary() -> size; 
    msg[2] = outList -> get_secondary() -> size;
    msg[3] = ps->get_buffer_size() - inList -> get_primary() -> size; 
    msg[4] = ps->get_buffer_size() - inList -> get_secondary() -> size;
    // worker send information to master
    comm->send_msg(master, &msg[0]);
    comm->os<<"worker chunk"<<local_chunks[0]<<"thread idle"<<ps->is_thread_idle(0)<<"all t wait"<<ps->all_thread_waiting()<<" "
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

void StupidWorker::work_thread(struct ProcStatus *ps, int region[4], int sppTask, int iter, int dev, int chunk, bool camera_ray){
    printf("image region %d %d %d %d %d\n", region[0], region[1], region[2], region[3], sppTask);
//    int sppProc = ps->get_tile_info(&region[0], process_time); 
//    int sppDev = sppProc;
    printf("width %d height%d spp %d iter + dev %d\n", ps->width, ps->height, sppTask, iter + dev);
    Settings settings {
        Vec3 { ps->eye.x, ps->eye.y, ps->eye.z },
        Vec3 { ps->dir.x, ps->dir.y, ps->dir.z },
        Vec3 { ps->up.x, ps->up.y, ps->up.z },
        Vec3 { ps->right.x, ps->right.y, ps->right.z },
        ps->w, ps->h,
        Vec4_i32 { region[0], region[1], region[2], region[3]},
        sppTask
    };
    render(&settings, iter + dev, dev, chunk, camera_ray);
}

void StupidWorker::run(float* frame_time) {
    
    //multi thread
    int deviceNum = ps->get_dev_num();
    int region[4];
    int spp = ps->spp; 
    int sppProc = ps->get_tile_info(&region[0], frame_time); 
    printf("stupid worker %d\n", comm->rank); 
    int sppDev = sppProc / deviceNum;
    int iter = 0; 
    if(ps->isRayQueuing()) {
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
            threadPool.emplace_back(std::thread(work_thread, ps, region, sppDev, comm->rank, i, chunk, iter == 0));
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
    int chunk_size = ps->get_chunk_size();
    List = new RayList *[chunk_size];
    for(int i = 0; i < chunk_size; i++) {
        List[i] = new RayList(ps->get_buffer_size() * 16, "out", false);
    }
    for(auto c : ps->get_local_chunk()) {
        List[c]->type = "in";
    }
    outList = List;
    comm->os << "out list size"<<outList[0]->size() << "capacity" << outList[0]->get_primary()->get_capacity() << std::endl; 
    comm->os << "loaded chunk()"<<ps->get_loaded_chunk()<< std::endl; 
    //when we need to load a new chunk updata inList pointer.
    inList  = List[ps->get_loaded_chunk()];
    buffer  = new RayList(ps->get_buffer_size(), "buffer", false);

    ps->set_chunks();
}

SmartWorker::SmartWorker(struct Communicator *comm, struct ProcStatus *ps) 
    : Worker(comm, ps) 
{
    int chunk_size = ps->get_chunk_size();
    statistic = new int[chunk_size * 4/*dep*/ * 2]; 
    std::fill(statistic, statistic + chunk_size * 4/*dep*/ * 2, 0);
    renderer_get_rays = 0;
    renderer_save_rays = 0;
    write_rays = 0; 
    sent_rays = 0;
}

SmartWorker::~SmartWorker() {
    delete buffer;
    for(int i = 0; i < ps->get_chunk_size(); i++){
        delete List[i];
    }
}

void SmartWorker::write_rays_buffer() {
//    comm->os<<"mthread write ray buffer"<<inList->size()<<"\n";
//    comm->os<<"mthread buffer size "<<buffer->size()<<"\n";
    if(inList->empty()) return;
    
    std::unique_lock <std::mutex> lock(buffer->mutex); 
    while(!buffer->empty()) {
        buffer_not_full.wait(lock);
    }
    comm->os<<"mthread get inlist lock\n";
    
    buffer->copy(inList, comm->rank);
    
//    if(buffer->get_primary()->empty() && !inList->get_primary()->empty()) {
//        comm->os<<"mthread before copy primaryary  "<<inList->get_primary()->size<<"|"<<inList->get_primary()->empty() 
//                         <<" buffer primary: "<<buffer->get_primary()->size<<"|"<< buffer->get_primary()->empty()<<"\n";
//        struct Rays *rays = buffer->get_primary();
//        rays->size = inList->get_primary()->copy_to_buffer(rays->get_data(), ps->get_buffer_size(), ps->get_buffer_capacity(), comm->rank, true);
//        comm->os<<"mthread after copy primaryary  "<<inList->get_primary()->size<<"|"<<inList->get_primary()->empty() 
//                         <<" buffer primary: "<<buffer->get_primary()->size<<"|"<< buffer->get_primary()->empty()<<"\n";
//        write_rays += buffer->get_primary()->size;
//    }
//    if(buffer->get_secondary()->empty() && !inList->get_secondary()->empty()) {
//        comm->os<<"mthread before copy secondaryary  "<<inList->get_secondary()->size<<"|"<<inList->get_secondary()->empty() 
//                         <<" buffer secondary: "<<buffer->get_secondary()->size<<"|"<< buffer->get_secondary()->empty()<<"\n";
//        comm->os<<"copy secondary\n";
//        struct Rays *rays = buffer->get_secondary(); 
//        rays->size = inList->get_secondary()->copy_to_buffer(rays->get_data(), ps->get_buffer_size(), ps->get_buffer_capacity(), comm->rank, false);
//        comm->os<<"mthread after copy secondaryary  "<<inList->get_secondary()->size<<"|"<<inList->get_secondary()->empty() 
//                         <<" buffer secondary: "<<buffer->get_secondary()->size<<"|"<< buffer->get_secondary()->empty()<<"\n";
//        write_rays += buffer->get_secondary()->size;
//    }
    
    buffer_not_empty.notify_one();
    lock.unlock();
    comm->os<<"mthread release inlist lock\n";
}

int SmartWorker::worker_load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) {
    if(comm->size == 1 || ps->Exit()) {
        comm->os << "rthread exit\n";
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
    while (thread_wait && buffer->empty() && !ps->Exit()) {
        comm->os<<"rthread wait for incoming lock"<<ps->is_thread_idle(thread_id)<<"\n";
        buffer_not_empty.wait(lock);
        comm->os<<"rthread get not empty condition" <<thread_wait<<" "<<buffer->empty()<<ps->Exit()<<"\n";
    }

    if(!queue->empty()) {
        ps->set_thread_idle(thread_id, false);
        int copy_size = queue->size; 
        memcpy(*rays, queue->data, ps->get_buffer_capacity() * width * sizeof(float)); 
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
    if(buffer->empty()) { 
        comm->os<<"rthread notify" <<buffer->size()<<"\n";
        buffer_not_full.notify_all();
        comm->os<<"rthread notify" <<buffer->size()<<"\n";
        
    }
    lock.unlock(); 

    if(ps->Exit()) {  
//        printf("ray size %ld   queue primary  size %d queue secondary size %d\n", 
//                rays_size, inList->get_primary()->size, inList->get_secondary()->size);
//        printf("recv stop mpi %d thread %d %ld\n", comm->rank, thread_id, rays_size);
        return -1;
    }
    return rays_size;
}

void SmartWorker::worker_save_outgoing_buffer(float *retired_rays, size_t size, size_t capacity, bool primary){
    comm->os<<"rthread save outgoing buffer"<<size<<"\n";
    renderer_save_rays += size;
    out_mutex.lock(); 
    int width = primary?21:14; 
    int* ids = (int*)(retired_rays);
    for(int i = 0; i < 5; i ++) {
        comm->os<<"| "<< ids[i] <<" "<< ids[i + ps->get_buffer_capacity() * 9] << " ";
    }
    comm->os<<"\n";
    RayList::classification(outList, retired_rays, size, capacity, primary, comm->rank);
    out_mutex.unlock(); 
}

void SmartWorker::work_thread(void* tmp, float *process_time, int devId, int devNum, bool preRendering) {
  
    SmartWorker * wk = (struct SmartWorker*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus *ps = wk->ps;
   
    int region[4]; 
    int sppProc = ps->get_tile_info(&region[0], process_time); 
    int sppDev = sppProc / devNum;
    int seed = comm->rank * devNum + devId;
    printf("width %d height%d spp %d dev id %d\n", ps->width, ps->height, sppDev, devId);
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
        render(&settings, seed, devId, wk->ps->get_loaded_chunk(), true);
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
    if(ps->all_thread_waiting() && buffer->empty() ) {
        if(all_queue_empty()) {
            comm->os<< "mthread if all thread idle all queue empty, set itself idle\n";
            ps->set_self_idle();
            if(ps->all_proc_idle() && ps->all_rays_received()) {
                int *s = ps->get_status();
                comm->os<< "exit all rays received "<< s[0] <<" "<< s[1] <<" "<<s[2]<<" "<<s[3]<<"\n";

                //if all worker end synchronous 
                comm->os<<"mthread all proc idle send collective\n";
                QuitMsg quit_msg(comm->rank); 
                comm->send_message(&quit_msg, ps);
                ps->set_exit();
                comm->os<<"set exit\n";
                buffer_not_empty.notify_all();
            } else {
                int *s = ps->get_status();
                comm->os<< "all rays received "<< s[0] <<" "<< s[1] <<" "<<s[2]<<" "<<s[3]<<"\n";
                comm->os<<"mthread proc idle send status\n";
                StatusMsg status(comm->rank, ps->get_status(), comm->size); 
                comm->send_message(&status, ps);
            }
            return true; 
        } else {
            //maybe load new chunk
            comm->os<<"need another chunk or need send rays\n";
            ps->set_proc_busy(comm->rank);
            return false; 
        }
    } else {
//        comm->os<<"mthread proc busy\n";
        ps->set_proc_busy(comm->rank);
        return false;
    }
} 

void SmartWorker::message_thread(void* tmp) {
  
    SmartWorker *wk = (struct SmartWorker*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus * ps = wk->ps;

    comm->os<<"message thread\n";
    int recv_count = 0;
    while (!ps->Exit()) {
       
        wk->check_rendering_status(); 
        
        if(ps->Exit()) break;
        bool recv = false;
//        if(comm->rank == 0 || wk->inList->get_secondary()->empty()) { 
            do {
                recv = comm->recv_message(wk->List, ps);
            } while(!recv && ps->is_proc_idle() && !ps->Exit()) ;
//        }

        wk->write_rays_buffer();
        
        if(ps->Exit()) break;

        int cId = get_sent_list(wk->outList, ps->get_chunk_size());
//        comm->os<<"mthread get send list"<<cId<<"\n";
        
        if(cId >= 0) {
//            ps->set_proc_busy(ps->get_dst_proc(cId));
            wk->out_mutex.lock();
            comm->os<<"mthread new RayMsg"<<cId<<"\n";
//            wk->sent_rays += wk->outList[cId]->size();
            RayMsg *ray_msg = new RayMsg(wk->outList[cId], comm->rank, ps->get_dst_proc(cId), cId, false); 
            comm->os<<"mthread RayMsg"<<ray_msg->get_chunk()<<"\n";
            wk->out_mutex.unlock(); 
            comm->send_message(ray_msg, ps);
        }
        comm->purge_completed_mpi_buffers();
//        comm->os<<"mthread single loop end \n";
//        usleep(4000);
	}
    comm->os <<" end message thread"<<ps->all_thread_waiting()<<"\n";
    comm->os <<" inlist "<< wk->inList->size()
             <<" buffersize"<<wk->buffer->size()
             <<" get render rays "<<wk->renderer_get_rays
             <<" save render rays" <<wk->renderer_save_rays
             <<" sent rays" << wk->sent_rays
             <<" writer_rays " << wk->write_rays
             <<" recv "<<ps->global_rays[ps->proc_rank + ps->proc_size]
             <<std::endl;
    wk->buffer_not_empty.notify_all();
    return;
} 

void SmartWorker::count_rays() {
    int width = ps->width; 
    int height = ps->height;
    int chunk_size = ps->get_chunk_size();
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
    int deviceNum = ps->get_dev_num();
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
        while(!ps->Exit()) {
            int chunk = ps->get_new_chunk(); 
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

void setup_master(struct Communicator *comm, struct ProcStatus *ps) {
    master.reset(new Master(comm, ps));
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

void setup_worker(struct Communicator *comm, struct ProcStatus *ps, bool with_master) {
    if(with_master) {
        worker.reset(new StupidWorker(comm, ps));
    } else {
        worker.reset(new SmartWorker(comm, ps));
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

int32_t worker_buffer_size() {
    return worker->proc_status()->get_buffer_size();
}
