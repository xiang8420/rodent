// Original Master-Worker mode 
struct AllCopyNode : public Node{
    double recv_t, wait_t, send_t, total_t, read_t, write_t, st, ed;

    std::mutex  block_mutex;
    std::condition_variable block_cond; //
    
    AllCopyNode(Communicator *comm, ProcStatus *ps, Scheduler*);
    
    ~AllCopyNode();
    
    void send_message();
    
    void save_outgoing_buffer(float *, size_t, bool){} ; 

    static void message_thread(void* tmp);
    
    int load_incoming_buffer(float **, bool, int) { }; 
    
    void run(Camera *);
    
};

AllCopyNode::AllCopyNode(Communicator *comm, ProcStatus *ps, Scheduler* scheduler)
    :Node(comm, ps, scheduler)
{}

AllCopyNode::~AllCopyNode() {
    printf("delete AllCopyNode\n");
}

void AllCopyNode::message_thread(void* tmp) {
  
    AllCopyNode *wk = (struct AllCopyNode*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus * ps = wk->ps;

    printf("%d message thread\n", comm->get_rank());
    comm->os<<"message thread\n";
    
    size_t block_count = wk->scheduler->get_block_count();
    int i = comm->get_size();
    printf("comm size %d block count %ld\n", comm->get_size(), block_count);
    while(i < block_count) {
        //printf("mthread schedule %d mthread\n", i);
        //master waiting
        if(ps->is_proc_idle()) {
            printf("master mthread self schedule %d\n", i);
            wk->scheduler->set_render_block(i++);
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

void AllCopyNode::run(Camera *cam) {
   
    ps->reset();
    scheduler->preprocess(cam, 1/*block size*/, false, false);
    comm->os <<" start image decomposition  \n";
    //start rendering with origin tile
    //recv camera
    int deviceNum = ps->get_dev_num();

    std::vector<std::thread> workThread;
    if(comm->isMaster()) {
        std::thread mthread(message_thread, this);
        launch_rodent_render(deviceNum, true);
        while(!ps->Exit()) {
            
            ps->set_proc_idle(comm->get_rank());
            printf("end rthread\n"); 
            std::unique_lock <std::mutex> lock(block_mutex); 
            while (!ps->Exit() && !ps->is_proc_idle()) {
                printf("wait for block\n");
                block_cond.wait(lock);
                printf("get block\n");
            }
            launch_rodent_render(deviceNum, true);
        }
        mthread.join();

    } else {
        launch_rodent_render(deviceNum, true);
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
            scheduler->set_render_block(new_block);
            launch_rodent_render(deviceNum, true);
        }
        printf("worker %d exit\n", comm->get_rank());
    }


    printf("run AllCopyNode\n");
    return;
}

