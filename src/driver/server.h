
struct Server {
    struct RayQueue **rayqueue;
    struct RayQueue *raybuffer;
    int recv_capacity, buffer_capacity, slave_num, width; 
    bool *slave_wait;
    
    struct Communicator *comm;

    int recv_num, send_num;

    double recv_t, wait_t, send_t, total_t, read_t, write_t, st, ed;

    Server(struct Communicator *comm):comm(comm) {
        slave_num = comm->size - 1;
        width = comm->width;
        recv_capacity = 1048608;
        buffer_capacity = 16 * recv_capacity;
        raybuffer = new struct RayQueue(recv_capacity, width); 
        rayqueue  = new struct RayQueue *[slave_num];
        slave_wait = new bool[slave_num];
        for(int i = 0; i < slave_num; i++){
            rayqueue[i] = new struct RayQueue(buffer_capacity, width);
            slave_wait[i] = false;
        }
        
        recv_num = 0; send_num = 0;
        
        st = clock();
        recv_t = 0; wait_t = 0;
        send_t = 0; total_t = 0;
        write_t = 0; read_t = 0;
    }

    ~Server(){
        for(int i = 0; i < slave_num; i++){
            delete rayqueue[i];
        }
        delete raybuffer;
    }; 

    // ray queue isempty and recv all slave end msg
    bool shutdown(){
        return false;
    }

    void run(){
        int msg[3];
        bool all_done = false;
//        return;
        for(;;){
            // master recv, true: recv a request  
            //             false: recv ray data
//            sleep(1);
            double t = clock();
            if(comm->Server_recv(raybuffer, &msg[0])){
                recv_t += (double)(clock() - t) / CLOCKS_PER_SEC;
                int id = msg[1];
                if(msg[0] == 0) {
                    // slave has no work
                    slave_wait[id] = true;
//                    printf("idle %d %d\n", slave_wait[0], slave_wait[1]);
                    if(rayqueue[id]->isempty()) {
                        // if slave wait and no rays in its queue
                        all_done = true;
//                        printf("%d rayqueue empty\n", id);
                        for(int i = 0; i < slave_num; i++){
                            if(!slave_wait[i] || !rayqueue[i]->isempty()) {
                                all_done = false; break;
                            }
                        }
                    }
                    if(all_done) {
                        break;
                    }
                }

                int dst = msg[1];  
                int request_size = msg[2];

                t = clock();
                if(!rayqueue[dst]->isempty()) {
                    slave_wait[dst] = false;
                    struct RayQueue *r = rayqueue[dst]; 
                    comm->Send_rays(rayqueue[dst], rayqueue[dst]->size, dst);
                } else {
                    comm->Send_noray(dst);
                }
                send_t += (double)(clock() - t) / CLOCKS_PER_SEC;
            } else {
                recv_t += (double)(clock() - t) / CLOCKS_PER_SEC;
                int recv_size = raybuffer->get_size();

                // arrange data
                int *bufferdata   = (int*) raybuffer->data;
                for(int i = 0; i < recv_size; i++){
                    int chunk_id = bufferdata[i * width + 9]  >> 12;
//                    if(chunk_id == 0 && i < 9){            
//                        printf("|%d : %d : %d |", bufferdata[i * width], chunk_id, bufferdata[i * width + 9]);
//                    }
                    rayqueue[chunk_id]->copy(raybuffer, i);
                    if(rayqueue[chunk_id]->isfull()){
                        printf("error queue is is full %d %d %d\n", chunk_id, rayqueue[chunk_id]->get_size(), rayqueue[chunk_id]->get_capacity());   
                        rayqueue[chunk_id]->clear();
                    }
                }
                for(int i = 0; i < slave_num;i++){
                    printf("%d %d \n", i, rayqueue[i]->get_size());
                }
                printf("\n");
                raybuffer->clear();
            }
 //           for(int i = 0; i < slave_num;i++){
 //               printf("%d %d \n", i, rayqueue[i]->get_size());
 //           }
 //           printf("\n");
        }
        printf("server send all over\n");
        for(int i = 0; i < slave_num;i++){
            printf("%d %d \n", i, rayqueue[i]->get_size());
        }
        printf("\n");
        for(int i = 0; i < slave_num; i++){
             comm->Send_end(i);
        }
        total_t = (double)(clock() - st) / CLOCKS_PER_SEC;
        printf( "\n Server recv time%f send time %f  total %f \n write %f read %fseconds\n", recv_t, send_t, total_t, write_t, read_t );
    }
};

