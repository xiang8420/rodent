
void render_spawn(struct Settings const* settings, int iter, int dev, int chunk, int *tag){
    render(settings, iter + dev, dev, chunk);
    *tag = 0;
}

struct Client{
    struct RayQueue *send_buffer;
    struct RayQueue *wait_send_buffer;
    struct Ray *wait_buffer;
    struct RayQueue *queue;

    int buffer_capacity, buffer_size, slave_num; 
    int width;
    int mpi_rank, mpi_size, server;
    bool has_work, stop_recv, first;

    std::mutex mutex;

    struct Communicator *comm;
//    struct Renderer *render; 
    
    Client(struct Communicator *comm, bool has_server): comm(comm){
        mpi_rank = comm->rank;
        mpi_size = comm->size; 
        width = comm->width; 
        buffer_capacity = 1048608;
        buffer_size = 1048576;
        if(has_server) {
            server = mpi_size - 1;
            slave_num = server;
        } else {
            server = -1;
        }
        first = true;
        has_work = true;
        stop_recv = false;
        wait_buffer      = new struct Ray(buffer_capacity, width);  
        wait_send_buffer = new struct RayQueue(buffer_capacity, width);
        send_buffer      = new struct RayQueue(buffer_capacity, width);
        queue            = new struct RayQueue(buffer_capacity * 4, width); 
    }

    ~Client(){
        delete wait_buffer;
        delete wait_send_buffer;
        delete send_buffer;
        delete queue;
        printf("client end\n");
    }

    static void request(void* tmp) {
        struct Client *p = (struct Client*)tmp;
        struct RayQueue *queue = p->queue;
        int queue_capacity = queue->get_capacity();
        
        struct Ray *wait_buffer = p->wait_buffer;
        int buffer_capacity = p->buffer_capacity;
        int buffer_size     = p->buffer_size;
        struct Communicator *comm = p->comm;
        //start recv thread and render
//        return;
        for(;;){
            usleep(100);
            p->mutex.lock();
            if(!queue->isempty() && wait_buffer->isempty()){
                int queue_size = queue->get_size();
                int st = queue_size > buffer_size ? queue_size - buffer_size : 0;
                int n  = st == 0 ? queue_size : buffer_size;
                wait_buffer->put(queue->data, st, n);
                queue->size -= wait_buffer->get_size();
                p->has_work = true;
            }
            bool work = p->has_work; 
            p->mutex.unlock();

//                printf("client send request %d %d %d\n", queue->get_size(), buffer_size, queue_capacity);
            if(queue->get_size() + buffer_size < queue_capacity && !p->stop_recv) {
                comm->Send_request(p->mpi_rank, work, p->server, queue->capacity - queue->size);
                int msg = comm->Recv_rays(queue, p->server);
                if(msg == 0){
                    printf("slave recv stop queue->size%d\n", queue->size);
                    p->stop_recv = true;
                } else if(msg == 2) {
//                    printf("%d get rays ", p->mpi_rank);
//                    for(int i = 0; i < 10; i ++) {
//                        printf("%d %f |", (int)(*queue)[i * width], (*queue)[i * width + 1]);
//                    }
//                    p->recv_ray_num += buffer_size;
                    printf("\n");
                }
            } 
            if (queue->isempty() && p->stop_recv){
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
            thread[dev_num] = new std::thread(request, (void*)this);
            thread[dev_num]->join();
        }
        for(int i = 0; i < dev_num; i++){
            thread[i]->join();
        }
    }

    void send_rays(float *rays, size_t size, size_t capacity, bool send_all){
        int *id = (int *) rays;
        int *chunk = (int *) &rays[capacity * 9];
        printf("send size %d\n", size);
        for(int i = 0; i < size; i++) {
            wait_send_buffer->put(rays, i, 1, capacity);
            int* d = (int*)wait_send_buffer->data;
            if(wait_send_buffer->isfull()) {
//                send_ray_num += wait_send_buffer->get_size(); 
                printf("send :");
                for(int j = 0; j < 10; j++){
                    printf("%d|", d[j * 21 + 9]);
                }
                printf("send buffer_size%d \n\n", wait_send_buffer->get_size());
//                if(!first) {comm->Wait(0);} else {first = false;}
                swap(send_buffer, wait_send_buffer); 
                comm->Send_rays(send_buffer, send_buffer->get_size(), server);
            }
        }
        if(send_all && !wait_send_buffer->isempty()){
//                if(!first) {comm->Wait(0);} else {first = false;}
            printf("send :");
            int* d = (int*)wait_send_buffer->data;
//            for(int j = 0; j < 10; j++){
//                printf(" %d %d * %d %d|", id[j], chunk[j], d[j * 21], d[j * 21 + 9]);
//            }
            swap(send_buffer, wait_send_buffer); 
            comm->Send_rays(send_buffer, send_buffer->get_size(), server);
        }
        printf("send over\n");
    }

    int recv_rays(float *rays, size_t rays_size, bool idle){
//        printf("%d client need rays\n", mpi_rank);
        do{
    //        printf("%d waiting  ", mpi_rank);
            mutex.lock();
    //        printf("%d lock  ", mpi_rank);
            if(!wait_buffer->isempty()){
   //             swap(rays, wait_buffer->data);
                int copy_size = wait_buffer->size < buffer_size - rays_size ? wait_buffer->size : buffer_size - rays_size; 
                int buffer_st = wait_buffer->size - copy_size;
                int rays_st   = rays_size;
                for(int i = 0; i < 21; i++ ){
                    for(int j = 0; j < copy_size; j++){
                        int src = buffer_st + j;
                        int dst = rays_st + j;
                        rays[i * 1048608 + dst] = wait_buffer->data[i * 1048608 + src];
                    }
                }
                
                if(mpi_rank == 0) {    
                    printf("ray size %d buffer size %d copy size %d buffer st %d rays st%d\n", rays_size, wait_buffer->size, copy_size, buffer_st, rays_st);
                    printf("copy size %d\n", copy_size);
                    int* ids = (int*)rays;
                    for(int i = 0; i < 10; i ++) {
                        printf("| %d %d *", ids[i + rays_st], ids[i + rays_st + 9 * 1048608]);
                    }
                    printf("\n");
                    printf("wait buffer %d\n", wait_buffer->size);
                }

                wait_buffer->size -= copy_size;
                has_work = true;
                mutex.unlock();
                return copy_size + rays_size;
            } 
            has_work = !idle || rays_size > 0 || !queue->isempty();
            mutex.unlock();
            if(stop_recv || !idle){
                if(mpi_rank == 0) {    
                    printf("ray size %d buffer size %d \n", rays_size, wait_buffer->size);
                }
                return rays_size;
            }
            usleep(4000);
        }while(idle);
    }
}; 



