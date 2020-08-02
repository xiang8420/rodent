
// Original Master-Worker mode 

struct Master : public Node{
    int worker_local_chunks[128][3];
    int comm_count, worker_size; 
    bool *mpi_worker_wait;
    double recv_t, wait_t, send_t, total_t, read_t, write_t, st, ed;

    Master(struct Communicator *comm, struct ProcStatus *ps);
    
    ~Master();
    
    int get_dst_worker_id();

    int get_max_rays_chunk(bool unloaded); 

    bool all_worker_finish(); 
    // ray queue empty and recv all worker end msg
    bool shutdown(){ return false; }

    void asyn_master();

    void rayQueuing(float * buffer, int size, int capacity); 
    
    void save_ray_batches(float *buffer, size_t size, size_t capacity, size_t thread_id); 
    
    void ray_batching_master(int sppTask, int film_width, int film_height); 

    void run(float* rProcessTime); 
    
    void save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary);
  
//    int load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait);
};


class Worker : public Node {
    
protected:    
    std::vector<int> local_chunks;
    bool current_chunk_empty;
public:    

    Worker(struct Communicator *comm, struct ProcStatus *ps);

    ~Worker();

    bool incoming_empty(){ return inList->empty() && buffer->empty(); }

    // send primary and secondary
    void save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary);
    
    void mpi_thread();

    static void message_thread(void* tmp); 
    
    void run(float* rProcessTime); 
}; 

Master::Master(struct Communicator *comm, struct ProcStatus *ps)
    :Node(comm, ps)
{
    int chunk_size = ps->get_chunk_size();
    worker_size    = comm->size - 1;
    comm_count     = 0;
    rayList       = new RayList *[chunk_size];
    comm->os<< "master chunk size"<<chunk_size<<"\n";
    mpi_worker_wait = new bool[worker_size];
    for(int i = 0; i < worker_size; i++) {
        mpi_worker_wait[i] = false;
    } 
    
    for(int i = 0; i < chunk_size; i++) {
        rayList[i] = new RayList(8 * ps->get_buffer_size(), "out", false);
        worker_local_chunks[i][0] = i;
    }
    buffer = new RayList(8 * ps->get_buffer_size(), "buffer", false); 
    
    st = clock();
    recv_t = 0; wait_t = 0;
    send_t = 0; total_t = 0;
    write_t = 0; read_t = 0;
}

Master::~Master(){
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++){
        delete rayList[i];
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
        int ray_size = rayList[i]->size(); 
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
        if(mpi_worker_wait[workerId] && rayList[workerChunk]->empty()) {
            comm->os<<"has queues ";
            for(int i = 0; i < chunk_size;i++) {
                comm->os<<i<<"| "<<rayList[i]->primary->size<<" "<<rayList[i]->secondary->size;
            }
            comm->os<<"\n";
            if(!rayList_empty()) {
                int tmp = workerChunk;
                workerChunk = get_max_rays_chunk(true);
                if(workerChunk != -1) {
                    msg[0] = workerChunk; 
                    msg[3] = std::min(msg[3], rayList[workerChunk]->primary->size); 
                    msg[4] = std::min(msg[4], rayList[workerChunk]->secondary->size);
                    comm->os<<"msg[3]"<<msg[3]<<" raylist primary"<<rayList[workerChunk]->primary->size
                            <<"msg[4]"<<msg[4]<<"  secondary "    <<rayList[workerChunk]->secondary->size<<"\n";
                    worker_local_chunks[workerId][0] = workerChunk;
                    comm->send_msg(workerId, &msg[0]);
                    recv_size = recv_size + msg[1] + msg[2];
                    comm->recv_rays(workerId, msg[1], buffer->get_primary());
                    RayList::classification(rayList, buffer->get_primary());
                    comm->recv_rays(workerId, msg[2], buffer->get_secondary());
                    RayList::classification(rayList, buffer->get_secondary());
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
        comm->os<<"msg[3]"<<msg[3]<<" raylist primary"<<rayList[workerChunk]->primary->size
                <<"msg[4]"<<msg[4]<<"  secondary "    <<rayList[workerChunk]->secondary->size<<"\n";
        
        msg[0] = workerChunk;
        msg[3] = std::min(msg[3], rayList[workerChunk]->primary->size); 
        msg[4] = std::min(msg[4], rayList[workerChunk]->secondary->size);
        
        comm->send_msg(workerId, &msg[0]);
        //if workerId idle send new rays and new scene id 
        recv_size = recv_size + msg[1] + msg[2];
        comm->recv_rays(workerId, msg[1], buffer->get_primary());
        RayList::classification(rayList, buffer->get_primary());
        comm->os<<"master classify primary\n";
        comm->recv_rays(workerId, msg[2], buffer->get_secondary());
        RayList::classification(rayList, buffer->get_secondary());
        
        comm->os<<"master send msg"<<msg[0]<<" "<<msg[1]<<" "<<msg[2]<<" "<<msg[3]<<" "<<msg[4];
        sent_size = sent_size + msg[3] + msg[4];
        comm->send_rays(workerId, msg[3], rayList[workerChunk]->primary);
        comm->send_rays(workerId, msg[4], rayList[workerChunk]->secondary);
        if(msg[3] > 1 || msg[4] > 0) {
            mpi_worker_wait[workerId] = false;
        } 
    }
    comm->os<<"master send all over\n";
    comm->os<<"master sent "<<sent_size<<" recv "<<recv_size<<"\n ";
    for(int i = 0; i < worker_size;i++){
        printf("%d %d \n", i, rayList[i]->primary->get_size());
    }
    printf("\n");
    for(int i = 0; i < worker_size;i++){
        printf("%d %d \n", i, rayList[i]->secondary->get_size());
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
        rayList[st_gId]->primary->read_device_buffer(buffer, st, ed - st, capacity, 0);
        st = ed;
    } 
}

void Master::save_ray_batches(float *buffer, size_t size, size_t capacity, size_t thread_id) {
    rayQueuing(buffer, size, capacity);   
    int chunk_size = ps->get_chunk_size();
    printf("%ld size %ld batch status \n", thread_id, size);
    for(int i = 0; i< chunk_size;i++) {
        printf("%d ", rayList[i]->primary->size); 
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
        if(workerChunk < 0 || (mpi_worker_wait[workerId] && rayList[workerChunk]->empty())) {
            int chunk_size = ps->get_chunk_size();
            for(int i = 0; i < chunk_size;i++) {
                printf("%d %d %d|", i, rayList[i]->primary->size, rayList[i]->secondary->size);
            }
            printf("\n");
            if(!rayList_empty()) {
                int tmp = workerChunk;
                workerChunk = get_max_rays_chunk(true);
                if(workerChunk != -1) {
                    msg[0] = workerChunk; 
                    worker_local_chunks[workerId][0] = workerChunk;
                    comm->send_msg(workerId, &msg[0]);
                    comm->recv_rays(workerId, msg[1], buffer->get_primary());
                    RayList::classification(rayList, buffer->get_primary());
                    comm->recv_rays(workerId, msg[2], buffer->get_secondary());
                    RayList::classification(rayList, buffer->get_secondary());
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
        msg[3] = std::min(msg[3], rayList[workerChunk]->primary->size); 
        msg[4] = std::min(msg[4], rayList[workerChunk]->secondary->size);
        
        comm->send_msg(workerId, &msg[0]);
        printf("master send msg %d %d %d %d %d\n", msg[0], msg[1], msg[2], msg[3], msg[4]);
        //if workerId idle send new rayList and new scene id 
        comm->recv_rays(workerId, msg[1], buffer->get_primary());
        RayList::classification(rayList, buffer->get_primary());
        comm->recv_rays(workerId, msg[2], buffer->get_secondary());
        RayList::classification(rayList, buffer->get_secondary());
        
        comm->send_rays(workerId, msg[3], rayList[workerChunk]->primary);
        comm->send_rays(workerId, msg[4], rayList[workerChunk]->secondary);
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

void Master::run(float * a) {
    printf("master->run()"); 
    if(ps->isRayQueuing()) {
        printf("batching\n");
        ray_batching_master(ps->spp, ps->width, ps->height);
    } else {
        printf("asyn master\n");
        asyn_master(); 
    }
}

void Master::save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary) { 
    std::cerr << "master can't save_outgoing_buffer\n";
}

//Worker

Worker::Worker(struct Communicator *comm, struct ProcStatus *ps) 
    : Node(comm, ps) {
    inList  = new RayList(ps->get_buffer_size() * 16, "in", false);
    outList = new RayList(ps->get_buffer_size() * 16, "out", false);
    buffer  = new RayList(ps->get_buffer_size(),     "buffer", false);
    current_chunk_empty = false; 
}

Worker::~Worker(){
    delete outList;
    delete inList;
    delete buffer; 
}

void Worker::save_outgoing_buffer(float *retired_rays, size_t size, size_t capacity, bool primary){
    int *id = (int *) retired_rays;
    int *chunk = (int *) &retired_rays[capacity * 9];
    struct Rays *list = primary ? outList->get_primary() : outList->get_secondary();
    outList->lock();
    list->read_device_buffer(retired_rays, 0, size, capacity, comm->rank);
    outList->unlock();
}


void Worker::mpi_thread() {
    outList->lock();
    int msg[MSG_SIZE];
    proc_idle =  ps->all_thread_waiting() && inout_list_empty();
    msg[0] = proc_idle ? -1 : 1;
    msg[1] = outList -> get_primary() -> size; 
    msg[2] = outList -> get_secondary() -> size;
    msg[3] = ps->get_buffer_size() - inList -> get_primary() -> size; 
    msg[4] = ps->get_buffer_size() - inList -> get_secondary() -> size;
    // worker send information to master
    comm->send_msg(master, &msg[0]);
    comm->os<<"worker chunk"<<local_chunks[0]<<"thread idle"<<ps->is_thread_idle(0)<<"all t wait"<<ps->all_thread_waiting()<<" "
           << "queue empty"<<inout_list_empty()<<" msg "<<msg[0]<<" "<<msg[1]<<" "<<msg[2]<<"inlist primaryi "<<inList->get_primary()->size 
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

void Worker::message_thread(void* tmp) {
    struct Worker *p = (struct Worker*)tmp;
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

//void Worker::work_thread(struct ProcStatus *ps, int region[4], int sppTask, int iter, int dev, int chunk, bool camera_ray){
//    printf("image region %d %d %d %d %d\n", region[0], region[1], region[2], region[3], sppTask);
////    int sppProc = ps->get_tile_info(&region[0], process_time); 
////    int sppDev = sppProc;
//    printf("width %d height%d spp %d iter + dev %d\n", ps->width, ps->height, sppTask, iter + dev);
//    Settings settings {
//        Vec3 { ps->eye.x, ps->eye.y, ps->eye.z },
//        Vec3 { ps->dir.x, ps->dir.y, ps->dir.z },
//        Vec3 { ps->up.x, ps->up.y, ps->up.z },
//        Vec3 { ps->right.x, ps->right.y, ps->right.z },
//        ps->w, ps->h,
//        Vec4_i32 { region[0], region[1], region[2], region[3]},
//        sppTask
//    };
//    render(&settings, iter + dev, dev, chunk, camera_ray);
//}

void Worker::run(float* frame_time) {
    
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
            threadPool.emplace_back(std::thread(work_thread, this, frame_time, i, deviceNum, false, iter == 0));
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
