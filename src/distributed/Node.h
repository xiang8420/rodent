class Node {

protected:
    struct RayList **rayList;
    struct RayList * outList;
    
    RayStreamList * inList;  
    
    Communicator *comm;
    ProcStatus *ps;
     
    std::mutex  out_mutex, buffer_mutex, in_mutex, thread_mutex;
    
    std::condition_variable inList_not_empty; //
    std::condition_variable render_start; //

public: 
    Node(Communicator *comm, ProcStatus *ps);

    ~Node();

    ProcStatus * proc_status(){return ps;}

    int load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait);
    
    static void work_thread(void* tmp, float* processTime, int devId, int devNum, bool preprocess, bool shoot_rays);
   
    virtual void save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary) = 0;

    virtual void run(float* rProcessTime) = 0;

    int get_sent_list(); 
    
    bool rayList_empty();
    
    bool outList_empty();
    
    bool inout_list_empty();
    
};

Node::Node(Communicator *comm, ProcStatus *ps) : comm(comm), ps(ps) { }

Node::~Node() {
    printf("%d delete Node\n", comm->rank);
}

bool Node::outList_empty() {  
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++) {
        if(rayList[i]->type == "out" && !rayList[i] -> empty() ) {
            return false;
        }
    }
    return true;
}

bool Node::rayList_empty() {  
    int chunk_size = ps->get_chunk_size();
    bool res = true;
    out_mutex.lock();
    for(int i = 0;i < chunk_size; i++)
        if(!rayList[i]->empty()) { res = false; break; }
    out_mutex.unlock();
    in_mutex.lock();
    res &= inList->empty();
    in_mutex.unlock();
    return res;
}

bool Node::inout_list_empty() {
    return inList->empty() && outList->empty();
}

int Node::load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) {
    if(comm->size == 1 || ps->Exit() || ps->has_new_chunk()) {
        comm->os << "rthread exit new chunk \n";
        return -1;
    }
    comm->os<<"rthread inlist size" <<inList->size()<<"\n";
    if(primary)
        comm->os<<"rthread primary\n";
    else 
        comm->os<<"rthread secondary\n";

    
    ps->set_thread_idle(thread_id, thread_wait);
    
    std::unique_lock <std::mutex> lock(inList->mutex); 
    while (thread_wait && inList->empty() && !ps->Exit() && !ps->has_new_chunk()) {
        comm->os<<"rthread wait for incoming lock"<<ps->is_thread_idle(thread_id)<<"\n";
        inList_not_empty.wait(lock);
        comm->os<<"rthread get condition" <<thread_wait<<" "<<ps->Exit()<<"\n";
    }
    
    if(ps->Exit() || ps->has_new_chunk()) {  
        return -1;
    }
    
    comm->os<<"rthread after wait\n";
    struct Rays *queue = primary ? inList->get_primary() : inList->get_secondary();
    lock.unlock();

    if(queue != NULL) {
        ps->set_thread_idle(thread_id, false);
        int copy_size = queue->size;
        int width = queue->store_width;
        comm->os<<"rthread width "<<width <<" logic width "<<queue->logic_width<<"\n";
        memcpy(*rays, queue->get_data(), ps->get_buffer_capacity() * width * sizeof(float)); 
        queue->size = 0;

        //printf("%d queue size %d %d %d\n", comm->rank, queue->size, primary, queue->store_width);
        int* ids = (int*)(*rays);
        comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
        for(int i = 0; i < 5; i ++) {
            comm->os<<"#" <<ids[i + rays_size]<<" "<<ids[i + rays_size + ps->get_buffer_capacity() * 9];
        }
        comm->os<<"\n";
        delete queue;
        return copy_size + rays_size;
    }
    return rays_size;
}

void Node::work_thread(void* tmp, float *process_time, int devId, int devNum, bool preRendering, bool generate_rays) {
    printf("work threadstart\n");

    Node * wk = (Node*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus *ps = wk->ps;
    
    int region[4]; 
    int sppProc = ps->get_tile_info(&region[0], process_time); 
    int sppDev = sppProc / devNum;
    int seed = comm->rank * devNum + devId;
    printf("width %d height%d spp %d dev id %d local chunk %d\n", ps->width, ps->height, sppDev, devId, ps->get_local_chunk() );
    printf("work thread region %d %d %d %d\n", region[0], region[1], region[2], region[3]);
    
    ImageDecomposition * camera = ps->get_camera();
    Settings settings {
        Vec3 { camera->eye.x, camera->eye.y, camera->eye.z },
        Vec3 { camera->dir.x, camera->dir.y, camera->dir.z },
        Vec3 { camera->up.x, camera->up.y, camera->up.z },
        Vec3 { camera->right.x, camera->right.y, camera->right.z },
        camera->w, camera->h,
        Vec4_i32 { region[0], region[1], region[2], region[3]},
        sppDev
    };
    if(preRendering) {
        prerender(&settings);
    } else {
        render(&settings, seed, devId, ps->get_local_chunk(), generate_rays);
    }
}

// get chunk retun chunk id
int Node::get_sent_list() {
    int n = ps->get_chunk_size();
    int t = -1; 
    int max = 0; 
    for(int i = 0; i < n; i++) {
        comm->os<<" get sent list "<< i<<" raylist size "<<rayList[i] -> size()<<" "<< ps->get_proc(i)<<"\n";
        if(rayList[i] -> size() > max && ps->get_proc(i) != -1 && i != ps->get_local_chunk()) {
            max = rayList[i] -> size();
            t = i;
        }
    }
    comm->os<<"mthread rank " <<comm->rank<<" list "<< t<<" max "<< max<<"\n";
    if(max <= 0 || t < 0) 
        return -1;
    
    return t;
}
