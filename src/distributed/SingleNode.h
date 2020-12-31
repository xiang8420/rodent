struct SingleNode : public Node{
    SingleNode(Communicator *comm, ProcStatus *ps, Scheduler *scheduler); 
    
    ~SingleNode();
    
    int load_incoming_buffer(float **, bool, int) { return -1; }; 
    
    void save_outgoing_buffer(float *, size_t, bool){} ; 
    
    void run(Scheduler * camera);

};

SingleNode::SingleNode( Communicator *comm, ProcStatus *ps, Scheduler * scheduler)
    :Node(comm, ps, scheduler)
{
    printf("new SingleNode\n");
    printf("Memory %ld kb \n", physical_memory_used_by_process());
}

SingleNode::~SingleNode() {
    printf("delete SingleNode\n");
}

void SingleNode::run(Scheduler * camera) {
    
    comm->os <<" start run message thread \n";
    int deviceNum = ps->get_dev_num();
    int iter = 0; 
    launch_rodent_render(camera, deviceNum, iter==0);
    printf("run SingleNode\n");
    return;
}

