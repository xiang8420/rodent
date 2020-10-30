class Node {

protected:
    //RayStreamList Memory pool
    Communicator *comm;
    ProcStatus *ps;
    
    RayStreamList inList;  
     
    std::mutex  out_mutex, buffer_mutex, in_mutex, thread_mutex;
    
    std::condition_variable inList_not_empty; //
    std::condition_variable render_start; //

    int msg_tag;
    bool sync;
public: 
    Node(Communicator *comm, ProcStatus *ps);

    ~Node();

    ProcStatus * get_proc_status(){return ps;}

    Communicator * get_communicator(){ return comm; }

    int load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait);
    
    static void work_thread(void* tmp, ImageDecomposition * camera, int devId, int devNum, bool preprocess, bool shoot_rays);

    virtual void run(ImageDecomposition * camera) = 0;

    virtual void save_outgoing_buffer(float *, size_t, bool) = 0;
    
    int get_tag();
    
};

Node::Node(Communicator *comm, ProcStatus *ps) : comm(comm), ps(ps) {
    msg_tag = 0;
    sync = false;
}

int Node::get_tag(){
    return msg_tag++;
}

Node::~Node() {
    printf("%d delete Node\n", comm->get_rank());
}

int Node::load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) {
    comm->os<<"rthread load incoming buffer inlist size "<<inList.size()<<"\n";
    if(inList.empty()) {
        if(sync)
            return -1;
        if(ps->Exit() || ps->has_new_chunk() || ps->get_chunk_size() == 1) {
            comm->os << "rthread exit new chunk \n";
            return -1;
        }
    } 
    statistics.start("run => work_thread => load_incoming_buffer-wait");
    if(!inList.empty() && ps->Exit())
        error("inlist not empty but prepared to exit\n");

    std::unique_lock <std::mutex> lock(inList.mutex); 
    if(!sync) {
        ps->set_thread_idle(thread_id, true);
        
        std::cout<<"rthread  inlist priamry " <<inList.primary_size()<<" secondary "<<inList.secondary_size()<<"\n";
        while (inList.empty() && !ps->Exit() && !ps->has_new_chunk()) {
            comm->os<<"rthread wait for incoming lock"<<ps->is_thread_idle(thread_id)<<"\n";
            inList_not_empty.wait(lock);
            comm->os<<"rthread get condition  "<<ps->Exit()<<"\n";
        }
        if(ps->Exit() || ps->has_new_chunk()) {  
            statistics.end("run => work_thread => load_incoming_buffer-wait");
            return -1;
        }
    }
    statistics.end("run => work_thread => load_incoming_buffer-wait");

    std::cout<<"primary ? "<<primary<<" primary size "<<inList.primary_size()<<" secondary size "<<inList.secondary_size()<<"\n";

    struct RaysStream *rays_stream;
    if(primary && inList.primary_size() > 0) {
        rays_stream = inList.get_primary();
    } else if (!primary && inList.secondary_size() > 0) {
        rays_stream = inList.get_secondary();
    } else {
        return rays_size;
    }

    lock.unlock();

    statistics.start("run => work_thread => load_incoming_buffer-copy");

    ps->set_thread_idle(thread_id, false);
    int copy_size = rays_stream->size;
    int width = rays_stream->width;
    printf("copy primary size %d\n", copy_size);
 //   comm->os<<"rthread width "<<width <<" logic width "<<rays_stream->logic_width<<"\n";
    memcpy(*rays, rays_stream->get_data(), ps->get_stream_capacity() * width * sizeof(float)); 
 
//     //printf("%d rays_stream size %d %d %d\n", comm->get_rank(), rays_stream->size, primary, rays_stream->store_width);
//    int* ids = (int*)(*rays);
//    comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
//    for(int i = 0; i < std::min(5, copy_size); i ++) {
//        comm->os<<"#" <<ids[i + rays_size]<<" "<<ids[i + rays_size + ps->get_buffer_capacity() * 9];
//    }
//    comm->os<<"\n";
    delete rays_stream;

    statistics.end("run => work_thread => load_incoming_buffer-copy");
    return copy_size + rays_size;
        
}

void Node::work_thread(void* tmp, ImageDecomposition * camera, int devId, int devNum, bool preRendering, bool generate_rays) {
    printf("work threadstart\n");

    Node * wk = (Node*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus *ps = wk->ps;
    
    int* region = camera->get_render_block();
    int sppProc = camera->get_spp(); 
    int sppDev = sppProc / devNum;
    
    printf("width %d height%d spp %d dev id %d local chunk %d\n", camera->width, camera->height, sppDev, devId, ps->get_local_chunk() );
    printf("work thread rank %d region %d %d %d %d local chunk %d\n", comm->get_rank(), region[0], region[1], region[2], region[3], ps->get_local_chunk());
    if(generate_rays)
        printf("generate_rays\n");
    else 
        printf("no ray generate\n");

    statistics.start("run => work_thread");
    Settings settings {
        Vec3 { camera->eye.x, camera->eye.y, camera->eye.z },
        Vec3 { camera->dir.x, camera->dir.y, camera->dir.z },
        Vec3 { camera->up.x, camera->up.y, camera->up.z },
        Vec3 { camera->right.x, camera->right.y, camera->right.z },
        camera->w, camera->h,
        Vec4_i32 { region[0], region[1], region[2], region[3]},
        sppDev,
        ps->get_rough_trace()
    };
    if(preRendering) {
        prerender(&settings);
    } else {
        int rnd = comm->get_rank() * devNum + devId;
        render(&settings, rnd, devId, ps->get_local_chunk(), generate_rays);
    }
    printf("work thread end\n");
    statistics.end("run => work_thread");
}
