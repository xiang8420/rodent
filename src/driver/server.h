
struct Server {
    struct RayQueue **queue_primary;
    struct RayQueue **queue_secondary;
    struct RayQueue *raybuffer;
    int recv_capacity, buffer_capacity, comm_obj, slave_num; 
    bool *slave_wait;

    struct Communicator *comm;
    int recv_num, send_num;

    double recv_t, wait_t, send_t, total_t, read_t, write_t, st, ed;

    Server(struct Communicator *comm):comm(comm) {
        slave_num = comm->size - 1;
        comm_obj  = 0;
        recv_capacity = 1048608;
        buffer_capacity = 8 * recv_capacity;
        queue_primary  = new struct RayQueue *[slave_num];
        queue_secondary = new struct RayQueue *[slave_num];
        slave_wait = new bool[slave_num];
         
        for(int i = 0; i < slave_num; i++){
            queue_primary[i] = new struct RayQueue(buffer_capacity, 21);
            queue_secondary[i] = new struct RayQueue(buffer_capacity, 14);
            slave_wait[i] = false;
        }
        raybuffer = new struct RayQueue(recv_capacity, 21); 
        
        recv_num = 0; send_num = 0;
        
        st = clock();
        recv_t = 0; wait_t = 0;
        send_t = 0; total_t = 0;
        write_t = 0; read_t = 0;
    }

    ~Server(){
        for(int i = 0; i < slave_num; i++){
            delete queue_primary[i];
            delete queue_secondary[i];
        }
        delete raybuffer;
    }; 
    bool all_queue_empty(int grid){
         return queue_primary[grid]->isempty()
             && queue_secondary[grid]->isempty();
    }

    int get_dst(){
        int dst = comm_obj++ % slave_num; 
        return dst; 
    }


    // ray queue isempty and recv all slave end msg
    bool shutdown(){
        return false;
    }

    void sort(RayQueue **queue){
        int width         = raybuffer->width; 
        int *bufferdata   = (int*) raybuffer->data;
        int recv_size     = raybuffer->size;
        for(int i = 0; i < recv_size; i++){
            int chunk_id = bufferdata[i * width + 9]  >> 12;
            if( i < 9){            
                printf("|%d : %d|", bufferdata[i * width], bufferdata[i * width + 9]);
            }
            queue[chunk_id]->copy(raybuffer, i);
            if(queue[chunk_id]->isfull()){
                printf("error queue is is full %d %d %d\n", chunk_id, queue[chunk_id]->get_size(), queue[chunk_id]->get_capacity());   
                queue[chunk_id]->clear();
            }
        }
        raybuffer->clear();
    }

    void run(){
        int msg[5];
        bool all_done = false;
//        return;
        int dst = 0;
           
        for(;;){
            int dst = get_dst();
            msg[0] = 1;  //
            msg[3] = queue_primary[dst]->size; //send primary size
            msg[4] = queue_secondary[dst]->size;//send secondary size ; 
//            printf("server send msg %d %d %d %d %d\n", msg[0], msg[1], msg[2], msg[3], msg[4]);
            comm->Send_msg(dst, &msg[0]);
            comm->Recv_msg(dst, &msg[0]);
//            printf("server recv msg %d %d %d %d %d\n", msg[0], msg[1], msg[2], msg[3], msg[4]);
            // msg[2] recv primary size msg[3] recv secondary size 
            
            comm->Recv_rays(dst, true, msg[1], raybuffer);
//            for(int i = 0; i < slave_num;i++){
//                printf("%d %d \n", i, queue_primary[i]->get_size());
//            }
//            printf("\n");
//            for(int i = 0; i < slave_num;i++){
//                printf("%d %d \n", i, queue_secondary[i]->get_size());
//            }
//            printf("\n");
            sort(queue_primary);
//            printf("server 1\n");
            comm->Recv_rays(dst, false, msg[2], raybuffer);
//            printf("server 2\n");
            sort(queue_secondary); 
//            printf("server 3\n");
            comm->Send_rays(dst, true, msg[3], queue_primary[dst]);
//            printf("server 4\n");
            comm->Send_rays(dst, false, msg[4], queue_secondary[dst]);
//            printf("server 5\n");
            if(msg[3] > 0 || msg[4] > 0) {
                slave_wait[dst] = false;
            }
            if(msg[0] == -1 /*slave idle*/ && msg[3] ==0 && msg[4] == 0) {
                slave_wait[dst] = true;
                all_done = true;
                for(int i = 0; i < slave_num; i++){
                    if(!slave_wait[i] || !all_queue_empty(i)) {
                        all_done = false; break;
                    }
                }
                if(all_done) {
                    break;
                }
            }
        }
        printf("server send all over\n");
        for(int i = 0; i < slave_num;i++){
            printf("%d %d \n", i, queue_primary[i]->get_size());
        }
        printf("\n");
        for(int i = 0; i < slave_num;i++){
            printf("%d %d \n", i, queue_secondary[i]->get_size());
        }
        printf("\n");
        for(int i = 0; i < slave_num; i++){
             comm->Send_end(i);
        }
        total_t = (double)(clock() - st) / CLOCKS_PER_SEC;
        printf( "\n Server recv time%f send time %f  total %f \n write %f read %fseconds\n", recv_t, send_t, total_t, write_t, read_t );


//        for(;;){
//            // master recv, true: recv a request  
//            //             false: recv ray data
////            sleep(1);
//            double t = clock();
//            if(comm->Server_recv(raybuffer, &msg[0])){
//                recv_t += (double)(clock() - t) / CLOCKS_PER_SEC;
//                int id = msg[1];
//                if(msg[0] == -1) {
//                    // slave has no work
//                    slave_wait[id] = true;
//                    if(all_queue_empty(id)) {
//                        // if slave wait and no rays in its queue
//                        all_done = true;
////                        printf("%d queue_primary empty\n", id);
////                        printf("idle proc:");
//                        for(int i = 0; i < slave_num; i++){
//                            if(!slave_wait[i] || !all_queue_empty(i)) {
//                                all_done = false; break;
//                            }
////                            printf("%d ",i);
//                        }
////                        printf("\n");
//                    }
//                    if(all_done) {
//                        break;
//                    }
//                }
//                int dst = msg[1];  
//                int request_primary_size = msg[2];
//                int request_secondary_size = msg[3];
//
//                t = clock();
//                if(request_primary_size > 0 && !queue_primary[dst]->isempty()) {
//                    
//                    int *bufferdata   = (int*)queue_primary[dst]->data;
//                    printf("server send  21 :");
//                    for(int i = 0; i < 10; i++){
//                        printf("|%d : %d", bufferdata[i * 21], bufferdata[i * 21 + 9]);
//                    }
//                    printf("\n");
//                    
//                    comm->Send_rays(queue_primary[dst], queue_primary[dst]->size, dst, true, true);
//                    slave_wait[dst] = false;
//                    printf("slave wait %d %d\n", slave_wait[0], slave_wait[1]);
//                } else if(request_secondary_size > 0 && !queue_secondary[dst]->isempty()) {
//                    
//                    int *bufferdata   = (int*)queue_secondary[dst]->data;
//                    printf("server send 14:" );
//                    for(int i = 0; i < 10; i++){
//                        printf("|%d : %d", bufferdata[i * 14], bufferdata[i * 14 + 9]);
//                    }
//                    printf("\n");
//                    
//                    comm->Send_rays(queue_secondary[dst], queue_secondary[dst]->size, dst, true, false);
//                    slave_wait[dst] = false;
//                    printf("slave wait %d %d\n", slave_wait[0], slave_wait[1]);
//                } else {
//                    comm->Send_noray(dst);
//                }
//                send_t += (double)(clock() - t) / CLOCKS_PER_SEC;
//            } else {
//                recv_t += (double)(clock() - t) / CLOCKS_PER_SEC;
//                int recv_size = raybuffer->get_size();
//                struct RayQueue **queue = raybuffer->width == 21 ? queue_primary : queue_secondary;
//                printf("recv width %d \n", raybuffer->width);
//
//                slave_wait[msg[1]] = (msg[0] == 2);
//                // arrange data
//                int *bufferdata   = (int*) raybuffer->data;
//                int width = raybuffer->width; 
////                printf("server recv :");
//                for(int i = 0; i < recv_size; i++){
//                    int chunk_id = bufferdata[i * width + 9]  >> 12;
////                    if( i < 9){            
////                        printf("|%d : %d|", bufferdata[i * width], bufferdata[i * width + 9]);
////                    }
//                    queue[chunk_id]->copy(raybuffer, i);
//                    if(queue[chunk_id]->isfull()){
//                        printf("error queue is is full %d %d %d\n", chunk_id, queue[chunk_id]->get_size(), queue[chunk_id]->get_capacity());   
//                        queue[chunk_id]->clear();
//                    }
//                }
//                printf("\n");
//                for(int i = 0; i < slave_num;i++){
//                    printf("%d %d %d \n", raybuffer->width, i, queue[i]->get_size());
//                }
//                printf("\n");
//                raybuffer->clear();
//            }
//        }
//        printf("server send all over\n");
//        for(int i = 0; i < slave_num;i++){
//            printf("%d %d \n", i, queue_primary[i]->get_size());
//        }
//        printf("\n");
//        for(int i = 0; i < slave_num;i++){
//            printf("%d %d \n", i, queue_secondary[i]->get_size());
//        }
//        printf("\n");
//        for(int i = 0; i < slave_num; i++){
//             comm->Send_end(i);
//        }
//        total_t = (double)(clock() - st) / CLOCKS_PER_SEC;
//        printf( "\n Server recv time%f send time %f  total %f \n write %f read %fseconds\n", recv_t, send_t, total_t, write_t, read_t );
    }
};

