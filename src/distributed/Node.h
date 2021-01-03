
class Node {

protected:
    Communicator *comm;
    ProcStatus *ps;
    Scheduler *scheduler;

    std::mutex  thread_mutex;
    
    std::condition_variable inList_not_empty; //
    std::condition_variable render_start; //

public: 
    int max_used_mem, pre_mem;
    
    Node(Communicator *comm, ProcStatus *ps, Scheduler *scheduler);

    ~Node();

    ProcStatus * get_proc_status(){return ps;}

    Communicator * get_communicator(){ return comm; }

    void launch_rodent_render(int, bool); 

    virtual void run() = 0;

    virtual int load_incoming_buffer(float **rays, bool primary, int thread_id) = 0;
    
    virtual void save_outgoing_buffer(float *, size_t, bool) = 0; 
    
    void loop_check(float i); 
};

Node::Node(Communicator *comm, ProcStatus *ps, Scheduler *scheduler) 
    : comm(comm), ps(ps), scheduler(scheduler)
{
    pre_mem = 0;
    max_used_mem = 0;
}

Node::~Node() {
    printf("%d delete Node\n", comm->get_rank());
}

void Node::launch_rodent_render(int devNum, bool generate_rays) {

    Camera *camera = scheduler->camera;
    int sppProc = scheduler->get_spp(); 
    int sppDev = sppProc / devNum;
    int cur_chk = ps->get_current_chunk();

    // only used for async dynamic schedule
    ps->chunk_loaded(); 

    std::vector<std::thread> workThread;
    for(int i = 0; i < devNum; i++) {
        int* region = scheduler->get_render_block();
        Settings settings {
            Vec3 { camera->eye.x, camera->eye.y, camera->eye.z },
            Vec3 { camera->dir.x, camera->dir.y, camera->dir.z },
            Vec3 { camera->up.x, camera->up.y, camera->up.z },
            Vec3 { camera->right.x, camera->right.y, camera->right.z },
            camera->w, camera->h,
            Vec4_i32 { region[0], region[1], region[2], region[3]},
            sppDev,
            ps->get_simple_trace()
        };
        int rnd = camera->iter * (comm->get_rank() + 1) * devNum + i;

        workThread.emplace_back(render, &settings, rnd, i, cur_chk, ps->get_next_chunk(), generate_rays);
    }

    render_start.notify_all(); 

    for(auto &thread: workThread)
        thread.join();
    
    comm->os<<comm->get_rank()<<" before set render start"<<"\n";
    // dynamic schedule for multi chunks    
    for(int i = 0; i < devNum; i++) 
        rodent_unload_chunk_data(cur_chk, i); 
}


void Node::loop_check(float i) {
    if(1) {    
//        int mem = physical_memory_used_by_process();
//        if(/*i> 1000 ||*/mem != pre_mem) {
//            comm->os<<i<<" memory use "<<mem<<"\n"; 
//            max_used_mem = std::max(mem, max_used_mem);
//            pre_mem = mem;
//        }
        //printf("%d mark %f\n", comm->get_rank(), i);
    }
}

