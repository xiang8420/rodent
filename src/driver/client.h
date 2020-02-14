
void render_spawn(struct Settings const* settings, int iter, int dev, int chunk, int *tag){
    render(settings, iter + dev, dev, chunk);
    *tag = 0;
}

struct Client{
    struct RayQueue * outcome_primary;
    struct RayQueue * income_primary;

    struct RayQueue * outcome_secondary;     
    struct RayQueue * income_secondary;           

    struct RayQueue * income_buffer_primary;
    struct RayQueue * income_buffer_secondary;

    int buffer_capacity, buffer_size, slave_num; 
    int mpi_rank, mpi_size, server;
    bool idle, stop_recv, first;

    std::mutex mutex, income_buffer_mutex;

    struct Communicator *comm;
//    struct Renderer *render; 
    
    Client(struct Communicator *comm, bool has_server): comm(comm){
        mpi_rank = comm->rank;
        mpi_size = comm->size; 
        buffer_capacity = 1048608;
        buffer_size = 1048576;
        if(has_server) {
            server = mpi_size - 1;
            slave_num = server;
        } else {
            server = -1;
        }
        first = true;
        idle = false;
        stop_recv = false;
        outcome_primary      = new struct RayQueue(buffer_capacity * 4, 21);
        income_primary     = new struct RayQueue(buffer_capacity * 4, 21); 
        
        outcome_secondary      = new struct RayQueue(buffer_capacity * 4, 14);
        income_secondary     = new struct RayQueue(buffer_capacity * 4, 14); 
        
        income_buffer_primary = new struct RayQueue(buffer_capacity, 21);
        income_buffer_secondary = new struct RayQueue(buffer_capacity, 14);
    }

    ~Client(){
        delete outcome_primary;
        delete income_primary;
        printf("client end\n");
    }

    bool all_empty(){ 
        return income_primary->empty() && income_secondary->empty()
            && outcome_primary->empty() && outcome_secondary->empty()
            && income_buffer_primary->empty() && income_buffer_secondary->empty();
    }
    
    bool all_queue_empty(){ 
        return income_primary->empty() 
            && income_secondary->empty();
    }

    // send primary and secondary
    void send_rays(float *rays, size_t size, size_t capacity, bool primary, bool send_all){
        int *id = (int *) rays;
        int *chunk = (int *) &rays[capacity * 9];
        printf("send size %d\n", size);
        struct RayQueue *buffer = primary ? outcome_primary : outcome_secondary;
        int width = primary? 21 : 14; 

        mutex.lock();
        for(int i = 0; i < size; i++) {
            buffer->put(rays, i, 1, capacity);
            if(buffer->full()) {
                printf("error : outcome buffer is full %d\n", width);
            }
        }
        int* d = (int*)buffer->rays();
        printf("%d client send %d %d:", mpi_rank, width, buffer->size);
        for(int j = 0; j < 10; j++){
            printf("%d %d|", d[j * width], d[j * width + 9]);
        }
        printf("send over\n");
        mutex.unlock();
    }
    
    void communicate_with_server(){
        int msg[5];
        comm->Recv_msg(server, &msg[0]);
//        printf("client recv msg %d %d %d %d %d\n", msg[0], msg[1], msg[2], msg[3], msg[4]);
        if(msg[0] == -1) {
            printf("slave recv stop queue->size%d \n", mpi_rank);
            stop_recv = true;
            return;
        }
        mutex.lock();
//        printf("proc %d %d %d %d %d\n", mpi_rank, idle, all_empty(),msg[3], msg[4]);
        idle = idle && all_empty() && msg[3] <= 0 &&msg[4] <= 0 ;
        msg[0] = idle ? -1 : 1;
        msg[1] = outcome_primary->size; 
        msg[2] = outcome_secondary->size;
        msg[3] = std::min(msg[3], buffer_size - income_primary->size); 
        msg[4] = std::min(msg[4], buffer_size - income_secondary->size); 
        
        comm->Send_msg(server, &msg[0]);
//        printf("client send msg %d %d %d %d %d\n", msg[0], msg[1], msg[2], msg[3], msg[4]);
//        if(outcome_primary->size > 0) {
//                int *bufferdata   = (int*)outcome_primary->data ;
//                printf("client buffer recv:");
//                for(int i = 0; i < 10; i++){
//                    printf("|%d :: %d", bufferdata[(0 + i) * 21], bufferdata[(0 + i) * 21 + 9]);
//                }
//                printf("\n");
//        }
        comm->Send_rays(server, true, msg[1], outcome_primary);
        comm->Send_rays(server, false, msg[2], outcome_secondary);
        comm->Recv_rays(server, true, msg[3], income_primary);
//        if(msg[3] > 0){
//                    int *bufferdata   = (int*)income_primary->rays() ;
//                    printf("client buffer recv:");
//                    for(int i = 0; i < 10; i++){
//                        printf("|%d :: %d", bufferdata[(0 + i) * 21], bufferdata[(0 + i) * 21 + 9]);
//                    }
//        }
        comm->Recv_rays(server, false, msg[4], income_secondary);
//        income_buffer_mutex.lock();
//        printf("income_buffer_primary %d %d %d %d %d\n", mpi_rank, income_buffer_primary->size, income_primary->size, income_buffer_secondary->size, income_secondary->size);
        if(income_buffer_primary->empty() && !income_primary->empty()) {
            struct RayQueue *queue = income_buffer_primary; 
            queue->size = income_primary->get(queue->rays(), buffer_size, buffer_capacity);
        }
        if(income_buffer_secondary->empty() && !income_secondary->empty()) {
            struct RayQueue *queue = income_buffer_secondary; 
            queue->size = income_secondary->get(queue->rays(), buffer_size, buffer_capacity);
        }
//        printf("income_buffer_primary %d %d %d %d %d\n", mpi_rank, income_buffer_primary->size, income_primary->size, income_buffer_secondary->size, income_secondary->size);
//        income_buffer_mutex.unlock();
        mutex.unlock();
    }


    static void Server_client(void* tmp) {
        struct Client *p = (struct Client*)tmp;
        struct RayQueue *queue = p->income_primary;
        //start recv thread and render
//        return;
        for(;;){
            usleep(2000);
            p->communicate_with_server();
            if (p->stop_recv){
                return;
            }
        }   
    }
    
    void run(struct Settings settings, int dev_num) {
              
        auto ticks = std::chrono::high_resolution_clock::now();
        int tag[] = {-1, -1, -1, -1, -1}; 
        int finished = 0;

        //multi thread
        std::thread *thread[5];
        int granularity = 1; 
        
        int chunk = mpi_rank;
        for(int i = 0; i < dev_num; i++){
            if(tag[i] <= 0){
                if(tag[i] == 0){
                    thread[i]->join();
                }
                tag[i] = 1;
                printf("work dev %d\n", i);
                thread[i] = new std::thread(render_spawn, &settings, mpi_rank, i, chunk, &tag[i]);
                finished ++;
            } 
        }
        if(server != -1){
            thread[dev_num] = new std::thread(Server_client, (void*)this);
            thread[dev_num]->join();
        }
        for(int i = 0; i < dev_num; i++){
            thread[i]->join();
        }
    }

    int recv_rays(float *rays, size_t rays_size, bool no_task, bool primary){
        struct RayQueue *queue = primary ? income_buffer_primary : income_buffer_secondary;
        int width = primary ? 21 : 14;
        do{
//            income_buffer_mutex.lock(); 
            mutex.lock();
//            printf("queue size %d %d\n", queue->size, primary);
            if(!queue->empty() && rays_size < buffer_size){
                idle = false;
                int copy_size = std::min(queue->size , buffer_size - (int)rays_size); 
                queue->size = queue->size - copy_size;
                int src = queue->size; 
                int dst = rays_size;
//                printf("ray size %d copy size %d queue st %d rays st%d\n", rays_size, copy_size, src, dst);
                for(int i = 0; i < width; i++) {
                    memcpy( &rays[dst + i * buffer_capacity], &queue->data[src + i * buffer_capacity], copy_size * sizeof(float));
                }
                    int* ids = (int*)rays;
                    printf("%d ray to renderi size %d copy size %d \n", mpi_rank, rays_size, copy_size);
                    for(int i = 0; i < 10; i ++) {
                        printf("| %d %d *", ids[i + rays_size], ids[i + rays_size + 1048608 * 9]);
                    }
                    printf("\n");
                    mutex.unlock();
//                income_buffer_mutex.unlock(); 
                return copy_size + rays_size;
            }
            idle = no_task && all_empty();
//            printf("proc %d rendering %d %d %d %d %d %d %d\n",mpi_rank, no_task, income_primary->size, income_secondary->size,
//                    outcome_primary->size, outcome_secondary->size, income_buffer_primary->size, income_buffer_secondary->size);
            income_buffer_mutex.unlock(); 
            mutex.unlock();
            if(stop_recv){  
                printf("ray size %d   queue primary  size %d queue secondary size %d\n", 
                        rays_size, income_primary->size, income_secondary->size);
                printf("recv stop %d %d\n", mpi_rank, rays_size);
                return -1;
            }
//            printf("%d %d\n", no_task, all_empty());
            if(!no_task || !all_empty() ){ return rays_size;}
            usleep(400);
        }while(no_task);
    }
}; 



