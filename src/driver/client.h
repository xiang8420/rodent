
void render_spawn(struct Settings const* settings, int iter, int dev, int chunk, int *tag){
    render(settings, iter + dev, dev, chunk);
    *tag = 0;
}

struct Client{
    struct RayQueue *send_buffer;
    struct RayQueue *wait_send_buffer;
    struct Ray *wait_buffer;
    struct RayQueue *queue;

    int buffer_size, slave_num; 
    int width;
    int mpi_rank, mpi_size, master;
    bool has_work, stop_recv, first;

    std::mutex mutex;

    struct Communicator *comm;
//    struct Renderer *render; 
    
    Client(struct Communicator *comm, bool has_server): comm(comm){
        mpi_rank = comm->rank;
        mpi_size = comm->size; 
        width = comm->width; 
        buffer_size = 1048608;
        if(has_server) {
            master = mpi_size - 1;
            slave_num = master;
        }
        first = true;
        has_work = true;
        stop_recv = false;
        wait_buffer      = new struct Ray(buffer_size, width);  
        wait_send_buffer = new struct RayQueue(buffer_size, width);
        send_buffer      = new struct RayQueue(buffer_size, width);
        queue            = new struct RayQueue(buffer_size * 4, width); 
    }
    ~Client(){}

    static void request(void* tmp) {
        struct Client *p = (struct Client*)tmp;
        struct RayQueue *queue = p->queue;
        int queue_capacity = queue->get_capacity();
        
        struct Ray *wait_buffer = p->wait_buffer;
        int buffer_size = p->buffer_size;
        struct Communicator *comm = p->comm;
        printf("client start request thread\n");        
        //start recv thread and render
        for(;;){
            usleep(300);
            p->mutex.lock();
            if(!queue->isempty() && wait_buffer->isempty()){
                printf("fill wait buffer\n");
                int queue_size = queue->get_size();
                int st = queue_size > buffer_size ? queue_size - buffer_size : 0;
                int n  = st == 0 ? queue_size : buffer_size;
                printf("wait buffer %d %d %d\n", queue_size, buffer_size, st);
               
                int id = p->width * st;
                wait_buffer->put(queue->data, id, n);
                queue->size -= wait_buffer->get_size();
            }
            p->mutex.unlock();

//                printf("client send request %d %d %d\n", queue->get_size(), buffer_size, queue_capacity);
            if(queue->get_size() + buffer_size < queue_capacity && !p->stop_recv) {
                comm->Send_request(p->mpi_rank, p->has_work, p->master, queue->capacity - queue->size);
                int msg = comm->Recv_rays(queue, p->master);
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
        thread[dev_num] = new std::thread(request, (void*)this);
        for(int i = 0; i <= dev_num; i++){
            thread[i]->join();
        }
    }

    void send_rays(float *rays, size_t size, size_t capacity, bool send_all){
        int *id = (int *) rays;
        int *chunk = (int *) &rays[capacity * 10];
        for(int i = 0; i < size; i++) {
            wait_send_buffer->put(rays, i, 1, capacity);
            int* d = (int*)wait_send_buffer->data;
            if(wait_send_buffer->isfull()) {
//                send_ray_num += wait_send_buffer->get_size(); 
//                for(int j = 0; j < 10; j++){
//                    printf("%d|", d[j * 20 + 10]);
//                }
                printf("send buffer_size%d \n\n", wait_send_buffer->get_size());
//                if(!first) {comm->Wait(0);} else {first = false;}
                swap(send_buffer, wait_send_buffer); 
                comm->Send_rays(send_buffer, send_buffer->get_size(), master);
            }
        }
        if(send_all && !wait_send_buffer->isempty()){
            printf("send buffer_size%d \n\n", wait_send_buffer->get_size());
//                if(!first) {comm->Wait(0);} else {first = false;}
            swap(send_buffer, wait_send_buffer); 
            comm->Send_rays(send_buffer, send_buffer->get_size(), master);
        }
        printf("send over\n");
    }

    int recv_rays(float *rays, size_t size){
        printf("client need rays\n");
        while(1){
            mutex.lock();
            if(stop_recv) {
                return 0;
            }
            if(!wait_buffer->isempty()){
                int size_new = wait_buffer->size;
   //             swap(rays, wait_buffer->data);
               
                 
                for(int i = 0; i < 20; i++ ){
                    for(int j = 0; j < size_new; j++){
                        rays[i * 1048608 + j] = wait_buffer->data[i * 1048608 + j];
                    }
                }
                int *id = (int *) rays;
                for(int i = 0; i < 10; i++) {
                    printf("%d %d|", id[i], id[i + wait_buffer->capacity * 10]);
                }
                
                wait_buffer->clear();
                has_work = true;
//                printf("before unlock");
                mutex.unlock();
//                printf("after lock");
                return size_new;
            } 
            has_work = false;
            mutex.unlock();
            usleep(100);
        }
    }
}; 



