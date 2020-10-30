// Original Master-Worker mode 
struct AllCopyNode : public Node{
    double recv_t, wait_t, send_t, total_t, read_t, write_t, st, ed;

    std::mutex  block_mutex;
    std::condition_variable block_cond; //
    
    AllCopyNode(struct Communicator *comm, struct ProcStatus *ps);
    
    ~AllCopyNode();
    
    void send_message();

    static void message_thread(void* tmp);
    
    void run(ImageDecomposition * camera);
    
    void save_outgoing_buffer(float *rays, size_t size, bool primary);

};

AllCopyNode::AllCopyNode(struct Communicator *comm, struct ProcStatus *ps)
    :Node(comm, ps)
{
    
    printf("new AllCopyNode\n");
    int chunk_size = ps->get_chunk_size();
    assert(chunk_size == 1 && comm->get_size() == chunk_size); 
    
    comm->os<< "master chunk size"<<chunk_size<<"\n";
  
    comm->os<<"master set up\n";
    printf("Memory %ld kb \n", physical_memory_used_by_process());
    
    st = clock();
}

AllCopyNode::~AllCopyNode() {
    printf("delete AllCopyNode\n");
}

void AllCopyNode::save_outgoing_buffer(float *retired_rays, size_t size, bool primary){ 
    error("Allcopynode don't need save outgoing");
}

void AllCopyNode::message_thread(void* tmp) {
  
    AllCopyNode *wk = (struct AllCopyNode*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus * ps = wk->ps;
    ImageDecomposition * camera = ps->get_camera();

    printf("%d message thread\n", comm->get_rank());
    comm->os<<"message thread\n";
    
    size_t block_count = camera->get_block_count();
    int i = comm->get_size();
    printf("comm size %d block count %ld\n", comm->get_size(), block_count);
    while(i < block_count) {
        //printf("mthread schedule %d mthread\n", i);
        //master waiting
        if(ps->is_proc_idle()) {
            printf("master mthread self schedule %d\n", i);
            camera->set_render_block(i++);
            ps->set_proc_busy(comm->get_rank());
            wk->block_cond.notify_all();
            continue;
        } else {
            MPI_Status status, s;
            int recv_ready = -1;
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_ready, &status);
            if(recv_ready) {
                int n = -1;
                MPI_Recv(&n, 1, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &s);
                printf("master mthread recv %d\n", n);
                MPI_Send(&i, 1, MPI_INT, n, 0, MPI_COMM_WORLD);
                printf("master mthread after send %d block %d\n", n, i);
                i++;
            }
        }
        usleep(100);
    }
    printf("end schedule1\n");
    ps->set_exit();
    wk->block_cond.notify_all();
    printf("end schedule2\n");
    int exit = -1;
    for(int i = 0; i < comm->get_size(); i++){
        if(i != comm->get_rank())
            MPI_Send(&exit, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
    printf("end schedule3\n");
    return;
} 

void AllCopyNode::run(ImageDecomposition * camera) {
   
    comm->os <<" start image decomposition  \n";
    //start rendering with origin tile
    //recv camera
    int deviceNum = ps->get_dev_num();

    ps->camera = camera;

    std::vector<std::thread> workThread;
    if(comm->isMaster()) {
        std::thread mthread(message_thread, this);
        for(int i = 0; i < deviceNum; i++)
            workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, true));

        for(auto &thread: workThread)
            thread.join();
        workThread.clear();
        while(!ps->Exit()) {
            
            ps->set_proc_idle(comm->get_rank());
            printf("end rthread\n"); 
            std::unique_lock <std::mutex> lock(block_mutex); 
            while (!ps->Exit() && !ps->is_proc_idle()) {
                printf("wait for block\n");
                block_cond.wait(lock);
                printf("get block\n");
            }
            for(int i = 0; i < deviceNum; i++)
                workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, true));

            for(auto &thread: workThread)
                thread.join();
            workThread.clear();
        }
        mthread.join();

    } else {
        for(int i = 0; i < deviceNum; i++) 
            workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, true));

        for(auto &thread: workThread)
            thread.join();
        
        workThread.clear();
        while(!ps->Exit()) {

            int new_block;
            MPI_Status status;
            printf("mthread worker idle %d\n", comm->get_rank());
            int comm_rank = comm->get_rank();
            MPI_Send(&comm_rank, 1, MPI_INT, comm->get_master(), 0, MPI_COMM_WORLD);
            MPI_Recv(&new_block, 1, MPI_INT, comm->get_master(), 0, MPI_COMM_WORLD, &status);
            printf("mthread worker  %d recv %d\n", comm->get_rank(), new_block);
            if(new_block == -1)
                break;
            camera->set_render_block(new_block);
            for(int i = 0; i < deviceNum; i++) 
                workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, true));

            for(auto &thread: workThread)
                thread.join();
            workThread.clear();
        }
        printf("worker %d exit\n", comm->get_rank());
    }


    printf("run AllCopyNode\n");
    return;
}

