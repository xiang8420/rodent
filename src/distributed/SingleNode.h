struct SingleNode : public Node{
    SingleNode(struct Communicator *comm, struct ProcStatus *ps);
    
    ~SingleNode();
    
    void run(ImageDecomposition * camera);

};

SingleNode::SingleNode(struct Communicator *comm, struct ProcStatus *ps)
    :Node(comm, ps)
{
    printf("new SingleNode\n");
    printf("Memory %ld kb \n", physical_memory_used_by_process());
}

SingleNode::~SingleNode() {
    printf("delete SingleNode\n");
}

void SingleNode::run(ImageDecomposition * camera) {
    
    comm->os <<" start run message thread \n";
    int deviceNum = ps->get_dev_num();
    int iter = 0; 
    std::vector<std::thread> workThread;
    printf("only one worker %d  ", comm->get_rank());
    for(int i = 0; i < deviceNum; i++) 
        workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, true));
    
    for( auto &thread: workThread) 
        thread.join();

    workThread.clear();
    printf("run SingleNode\n");
    return;
}

