
struct Worker{
    
    struct RayQueue * out_primary_queue;
    struct RayQueue * out_secondary_queue;     

    struct RayQueue * in_primary_queue;
    struct RayQueue * in_secondary_queue;           

    struct RayQueue * exchange_primary;
    struct RayQueue * exchange_secondary;

    int buffer_capacity, buffer_size, slave_num; 
    int mpi_rank, mpi_size, master;
    std::vector<int> local_chunks;
    bool idle, first, stop_recv;

    std::mutex outcome_buffer_mutex, income_buffer_mutex;

    int recv_loop_count = 0, master_loop_count = 0;
    struct Communicator *comm;

    Worker(struct Communicator *comm, bool has_master): comm(comm){
        mpi_rank = comm->rank;
        mpi_size = comm->size;
        stop_recv = false; 
        buffer_capacity = 1048608;
        buffer_size = 1048576;
       
        local_chunks.emplace_back(mpi_rank);

        if(has_master) {
            master = mpi_size - 1;
            slave_num = master;
        } else {
            master = -1;
        }
        first = true;
        idle = false;
        out_primary_queue      = new struct RayQueue(buffer_capacity * 4, 21);
        in_primary_queue     = new struct RayQueue(buffer_capacity * 4, 21); 
        
        out_secondary_queue      = new struct RayQueue(buffer_capacity * 4, 14);
        in_secondary_queue     = new struct RayQueue(buffer_capacity * 4, 14); 
        
        exchange_primary = new struct RayQueue(buffer_capacity, 21);
        exchange_secondary = new struct RayQueue(buffer_capacity, 14);
    }

    ~Worker(){
        delete out_primary_queue;
        delete in_primary_queue;
        printf("worker end\n");
    }

    bool all_empty(){ 
        return in_primary_queue->empty() && in_secondary_queue->empty()
 //           && out_primary_queue->empty() && out_secondary_queue->empty()
            && exchange_primary->empty() && exchange_secondary->empty();
    }
    
    bool all_queue_empty(){ 
        return in_primary_queue->empty() 
            && in_secondary_queue->empty();
    }

    // send primary and secondary
    void Write_outcome_buffer(float *rays, size_t size, size_t capacity, bool primary, bool send_all){
        int *id = (int *) rays;
        int *chunk = (int *) &rays[capacity * 9];
        printf("send size %d\n", size);
        struct RayQueue *buffer = primary ? out_primary_queue : out_secondary_queue;
        int width = primary? 21 : 14; 

        outcome_buffer_mutex.lock();
        for(int i = 0; i < size; i++) {
            buffer->put(rays, i, 1, capacity);
            if(buffer->full()) {
                printf("error : outcome buffer is full %d\n", width);
            }
        }
//        int* d = (int*)buffer->rays();
//        printf("%d worker send %d %d:", mpi_rank, width, buffer->size);
//        for(int j = 0; j < 10; j++){
//            printf("%d %d|", d[j * width], d[j * width + 9]);
//        }
//        printf("send over\n");
        outcome_buffer_mutex.unlock();
    }
  
    int Read_income_buffer(float **rays, size_t rays_size, bool no_task, bool primary) {
        struct RayQueue *queue = primary ? exchange_primary : exchange_secondary;
//        printf("%d read income buffer\n", mpi_rank);
        int width = primary ? 21 : 14;
        do {
            income_buffer_mutex.lock(); 
            if(!queue->empty() && rays_size < buffer_size) {
                printf("%d queue size0 %d %d\n", mpi_rank, queue->size, primary);
                idle = false;
                int copy_size;
                if(rays_size == 0) {
                    copy_size = queue->size; 
                    queue->size = 0;
                    memcpy(*rays, queue->data, buffer_capacity * width * sizeof(float)); 
                    printf("%d queue size %d %d\n", mpi_rank, queue->size, primary);
                } else {
                    copy_size = std::min(queue->size , buffer_size - (int)rays_size); 
                    queue->size = queue->size - copy_size;
                    int src = queue->size; 
                    int dst = rays_size;
                    for(int i = 0; i < width; i++) {
                        memcpy( &rays[dst + i * buffer_capacity], &queue->data[src + i * buffer_capacity], copy_size * sizeof(float));
                    }
                }
                income_buffer_mutex.unlock(); 
                return copy_size + rays_size;
            }
            idle = no_task && all_empty();
            income_buffer_mutex.unlock(); 
            if(stop_recv){  
                printf("ray size %d   queue primary  size %d queue secondary size %d\n", 
                        rays_size, in_primary_queue->size, in_secondary_queue->size);
                printf("recv stop %d %d\n", mpi_rank, rays_size);
                return -1;
            }
            if(!idle){ return rays_size;}
            usleep(300);
            recv_loop_count ++;
        } while(no_task);
    }
    
    void Communicate_with_master(){
        outcome_buffer_mutex.lock();
        int msg[MSG_SIZE];
        printf("idle %d", idle);
        idle = idle && all_empty();
        msg[0] = idle ? -1 : 1;
        msg[1] = out_primary_queue->size; 
        msg[2] = out_secondary_queue->size;
        msg[3] = buffer_size - in_primary_queue->size; 
        msg[4] = buffer_size - in_secondary_queue->size;
        // worker send information to master 
        comm->Send_msg(master, &msg[0]);
        printf("worker send msg %d %d %d %d %di %d %d\n", msg[0], msg[1], msg[2], 
                            in_primary_queue->size, in_secondary_queue->size, exchange_primary->size, exchange_secondary->size);
        comm->Recv_msg(master, &msg[0]);
        
        // finish or need new chunk
        if(msg[0] == -1 ) {
            printf("slave recv stop queue->size%d %d %d new chunk %d\n", mpi_rank, master_loop_count, recv_loop_count, msg[0]);
            local_chunks[0] = -1;
            outcome_buffer_mutex.unlock();
            stop_recv = true;
            return;
        }
        if(msg[0] != local_chunks[0]) {
            printf("%d recv new chunk %d\n", mpi_rank, msg[0]);
            local_chunks[0] = msg[0];
            comm->Send_rays(master, true, msg[1], out_primary_queue);
            comm->Send_rays(master, false, msg[2], out_secondary_queue);
            outcome_buffer_mutex.unlock();
            printf("%d after recv new chunk %d\n", mpi_rank, msg[0]);
            stop_recv = true;
            return;
        } 
        comm->Send_rays(master, true, msg[1], out_primary_queue);
        comm->Send_rays(master, false, msg[2], out_secondary_queue);
        outcome_buffer_mutex.unlock();
        comm->Recv_rays(master, true, msg[3], in_primary_queue);
        comm->Recv_rays(master, false, msg[4], in_secondary_queue);
        income_buffer_mutex.lock();
//        printf("exchange_primary %d %d %d %d %d\n", mpi_rank, exchange_primary->size, in_primary_queue->size, exchange_secondary->size, in_secondary_queue->size);
        if(exchange_primary->empty() && !in_primary_queue->empty()) {
            struct RayQueue *queue = exchange_primary;
            queue->size = in_primary_queue->get(queue->rays(), buffer_size, buffer_capacity);
        }
        if(exchange_secondary->empty() && !in_secondary_queue->empty()) {
            struct RayQueue *queue = exchange_secondary; 
            queue->size = in_secondary_queue->get(queue->rays(), buffer_size, buffer_capacity);
        }
        income_buffer_mutex.unlock();
        master_loop_count ++;
        usleep(100);
    }

    static void Comm_proc(void* tmp) {
        struct Worker *p = (struct Worker*)tmp;
        struct RayQueue *queue = p->in_primary_queue;
        //start recv thread and render
//        return;
        for(;;) {
            p->Communicate_with_master();
            if (p->stop_recv){
                printf("worker comm proc return\n ");
                return;
            }
        }   
    }
    
    static void render_proc(struct Settings const* settings, int iter, int dev, int chunk, bool valid_camera, int *tag){
        render(settings, iter + dev, dev, chunk, valid_camera);
        *tag = 0;
    }

    void run(struct Settings settings, int dev_num) {
        auto ticks = std::chrono::high_resolution_clock::now();
        int tag[] = {-1, -1, -1, -1, -1}; 
        int finished = 0;
        bool generate_rays = true;
        //multi thread
        std::thread *thread[5];
        while(true){
            stop_recv = false;
            int chunk = local_chunks[0];
            printf("%d load chunk %d", mpi_rank, chunk);
            if(chunk == -1) {
                return;
            } 
            for(int i = 0; i < dev_num; i++){
                if(tag[i] <= 0){
                    if(tag[i] == 0) {
                        thread[i]->join();
                    }
                    tag[i] = 1;
                    thread[i] = new std::thread(render_proc, &settings, mpi_rank, i, chunk, generate_rays, &tag[i]);
                    generate_rays = false;
                    finished ++;
                } 
            }
            if(master != -1){
                thread[dev_num] = new std::thread(Comm_proc, (void*)this);
            }
            for(int i = 0; i <= dev_num; i++){
                thread[i]->join();
                printf("%d join %d \n", mpi_rank, i);
            }
            for(int i = 0; i <= dev_num; i++) {tag[i] = -1;}
        }
    }
}; 



