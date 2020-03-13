struct ChunkRays {
    struct RayQueue *primary;
    struct RayQueue *secondary;
    std::vector<int> worker_allocated;
    int capacity;

    ChunkRays(int n):capacity(n) {
        primary = new struct RayQueue(capacity, 21);
        secondary = new struct RayQueue(capacity, 14);
    }
    
    ~ChunkRays() {
        delete primary;
        delete secondary;
    }

    int raySize() {
        return primary->size + secondary->size;
    }

    bool empty() {
        return primary->empty() && secondary->empty();
    } 
};

struct Master {
    struct Communicator *comm;
    struct ChunkRays **rays;
    struct RayQueue *raybuffer;
    int worker_local_chunks[128][3];
    int recv_capacity, buffer_capacity, comm_count;
    int chunk_size, worker_size; 
    bool *worker_wait;
    
    double recv_t, wait_t, send_t, total_t, read_t, write_t, st, ed;

    Master(struct Communicator *comm, int chunk_size):comm(comm), chunk_size(chunk_size) {
        worker_size    = comm->size - 1;
        comm_count     = 0;
        rays           = new struct ChunkRays *[chunk_size];
        worker_wait = new bool[worker_size];
        for(int i = 0; i < worker_size; i++) {
            worker_wait[i] = false;
        } 
        
        recv_capacity = 1048608;
        buffer_capacity = 8 * recv_capacity;
        for(int i = 0; i < chunk_size; i++) {
            rays[i] = new struct ChunkRays(buffer_capacity);
            worker_local_chunks[i][0] = i;
        }
        raybuffer = new struct RayQueue(recv_capacity, 21); 
        
        st = clock();
        recv_t = 0; wait_t = 0;
        send_t = 0; total_t = 0;
        write_t = 0; read_t = 0;
    }

    ~Master(){
        for(int i = 0; i < chunk_size; i++){
            delete rays[i];
        }
        delete raybuffer;
    }; 
    
    int get_dst_worker_id(){
        int dst = comm_count++ % worker_size; 
        return dst; 
    }

    int get_max_rays_chunk(bool unloaded) {
        int max_ray = 0, id_max = 0;
        printf("get max rays:");
        for(int i = 0; i < chunk_size; i++) {
            int ray_size = rays[i]->raySize(); 
            printf("|%d ", ray_size);
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
        printf("\n");
        if(max_ray == 0) return -1;
        return id_max;
    }

    bool all_queue_empty() {
        for(int i = 0; i < chunk_size; i++) {
            if(!rays[i]->empty()) {return false;}
        }
        return true;
    }

    bool all_worker_finish() {
        bool finish = true;
        for(int i = 0; i < worker_size; i++){
            if(!worker_wait[i]) {
                finish = false; break;
            }
        }
        return finish;
    }
    // ray queue empty and recv all worker end msg
    bool shutdown(){
        return false;
    }

    void raysClassify(bool isPrimary){
        int width         = raybuffer->width; 
        int *bufferdata   = (int*) raybuffer->rays();
        int recv_size     = raybuffer->size;
        if(isPrimary) {
            for(int i = 0; i < recv_size; i++){
                int chunk = bufferdata[i * width + 9]  >> 12;
                rays[chunk]->primary->copy(raybuffer, i);
                if(rays[chunk]->primary->full()){
                    printf("error queue is is full %d %d \n", chunk, rays[chunk]->primary->get_size());   
                }
            }
        } else {
            for(int i = 0; i < recv_size; i++){
                int chunk = bufferdata[i * width + 9]  >> 12;
                rays[chunk]->secondary->copy(raybuffer, i);
                if(rays[chunk]->secondary->full()){
                    printf("error queue is is full %d %d \n", chunk, rays[chunk]->secondary->get_size());   
                }
            }
        }
        raybuffer->clear();
    }

    void run() {
        int msg[MSG_SIZE];
        int  stopped = 0; 
//        return;
        for(;;) {
            int workerId = get_dst_worker_id();
            int workerChunk = worker_local_chunks[workerId][0]; 
            comm->Recv_msg(workerId, &msg[0]);
            printf("master recv msg %d %d %d %d %d\n", msg[0], msg[1], msg[2], msg[3], msg[4]);
            worker_wait[workerId] = (msg[0] == -1);
           
            for(int i = 0; i < chunk_size;i++) {
                printf("%d %d %d|", i, rays[i]->primary->size, rays[i]->secondary->size);
            }
            for(int i = 0; i< worker_size;i++) {
                printf("(%d %d)", i, worker_local_chunks[i][0]);
            }
            printf("\n");
            // scene chunk schedule
            if(worker_wait[workerId] && rays[workerChunk]->empty()) {
                printf("has queues ");
                for(int i = 0; i < chunk_size;i++) {
                    printf("%d %d %d|", i, rays[i]->primary->size, rays[i]->secondary->size);
                }
                printf("\n");
                if(!all_queue_empty()) {
                    int tmp = workerChunk;
                    workerChunk = get_max_rays_chunk(true);
                    if(workerChunk != -1) {
                        msg[0] = workerChunk; 
                        worker_local_chunks[workerId][0] = workerChunk;
                        comm->Send_msg(workerId, &msg[0]);
                        
                        comm->Recv_rays(workerId, true, msg[1], raybuffer);
                        raysClassify(true);
                        comm->Recv_rays(workerId, false, msg[2], raybuffer);
                        raysClassify(false);
                        printf("master send %d need new chunk%d", workerId, workerChunk);
                        worker_wait[workerId] = false; 
                        continue;     
                    } else {
                        printf("no chunk need to be loaded, use old chunk \n");
                        workerChunk = tmp; 
                        // Nothing happened, the worker still idle
                    }
                } else if(all_worker_finish()) {
                    printf("send end %d", workerId);
                    comm->Send_end(workerId);
                    stopped ++;
                    printf("master stop %d\n", stopped);
                    if(stopped == worker_size) {break;}
                    
                    continue;
                } 
            } 
            msg[0] = workerChunk;
            msg[3] = std::min(msg[3], rays[workerChunk]->primary->size); 
            msg[4] = std::min(msg[4], rays[workerChunk]->secondary->size);
            
            comm->Send_msg(workerId, &msg[0]);
//            printf("master send msg %d %d %d %d %d\n", msg[0], msg[1], msg[2], msg[3], msg[4]);
            //if workerId idle send new rays and new scene id 
            comm->Recv_rays(workerId, true, msg[1], raybuffer);
            raysClassify(true);
            comm->Recv_rays(workerId, false, msg[2], raybuffer);
            raysClassify(false);
            
            comm->Send_rays(workerId, true, msg[3], rays[workerChunk]->primary);
            comm->Send_rays(workerId, false, msg[4], rays[workerChunk]->secondary);
            if(msg[3] > 1 || msg[4] > 0) {
                worker_wait[workerId] = false;
            } 
        }
        printf("master send all over\n");
        for(int i = 0; i < worker_size;i++){
            printf("%d %d \n", i, rays[i]->primary->get_size());
        }
        printf("\n");
        for(int i = 0; i < worker_size;i++){
            printf("%d %d \n", i, rays[i]->secondary->get_size());
        }
        printf("\n");
        total_t = (double)(clock() - st) / CLOCKS_PER_SEC;
        printf( "\n Master recv time%f send time %f  total %f \n write %f read %fseconds\n", recv_t, send_t, total_t, write_t, read_t );
    }
};

