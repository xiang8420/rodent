class Node {

protected:
//    struct RayList **rayList;
    std::vector <RayList> rayList;
    RayStreamList inList;  
    
    Communicator *comm;
    ProcStatus *ps;
     
    std::mutex  out_mutex, buffer_mutex, in_mutex, thread_mutex;
    
    std::condition_variable inList_not_empty; //
    std::condition_variable render_start; //

    int msg_tag;
public: 
    Node(Communicator *comm, ProcStatus *ps);

    ~Node();

    ProcStatus * get_proc_status(){return ps;}

    Communicator * get_communicator(){ return comm; }

    int load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait);
    
    static void work_thread(void* tmp, ImageDecomposition * camera, int devId, int devNum, bool preprocess, bool shoot_rays);
   
    virtual void save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary) = 0;

    virtual void run(ImageDecomposition * camera) = 0;

    int get_sent_list(); 
    
    bool rayList_empty();
    
    bool outList_empty();
    
    bool inout_list_empty();

    int get_tag();
    
};

Node::Node(Communicator *comm, ProcStatus *ps) : comm(comm), ps(ps) {
    msg_tag = 0;
    int chunk_size = ps->get_chunk_size();
    rayList.resize(chunk_size); //        = new RayList *[chunk_size];
//    rayList = new RayList *[chunk_size];
    for(int i = 0; i < chunk_size; i++) {
        rayList.emplace_back(RayList(0, "out"));
        //rayList[i].set_capacity(1048608, "out");
    }
    inList.set_capacity(1048576/*ps->get_buffer_size()*/);
  //  for(int i = 0; i < ps->get_chunk_size(); i++) { 
  //      printf("new out going ray data size %ld capacity %ld\n",rayList[i].secondary.data.size(), rayList[i].secondary.data.capacity());
  //      printf("new out going ray data size %ld capacity %ld\n",rayList[i].primary.data.size(), rayList[i].primary.data.capacity());
  //  }
}

int Node::get_tag(){
    return msg_tag++;
}

Node::~Node() {
    printf("%d delete Node\n", comm->rank);
}

bool Node::outList_empty() {  
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++) {
        if(rayList[i].type == "out" && !rayList[i].empty() ) {
            return false;
        }
    }
    return true;
}

bool Node::rayList_empty() {  
//    printf("if raylist empty\n");
    int chunk_size = ps->get_chunk_size();
    bool res = true;
    out_mutex.lock();
    for(int i = 0;i < chunk_size; i++)
        if(!rayList[i].empty()) { res = false; break; }
    out_mutex.unlock();
    in_mutex.lock();
    res &= inList.empty();
    in_mutex.unlock();
    return res;
}

bool Node::inout_list_empty() {
    return inList.empty();
}

int Node::load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) {
    if(comm->size == 1 || ps->Exit() || ps->has_new_chunk() || ps->get_chunk_size() == 1) {
        comm->os << "rthread exit new chunk \n";
        return -1;
    }
    comm->os<<"rthread inlist size" <<inList.size()<<"\n";
    
    ps->set_thread_idle(thread_id, thread_wait);
    
    std::unique_lock <std::mutex> lock(inList.mutex); 
    while (thread_wait && inList.empty() && !ps->Exit() && !ps->has_new_chunk()) {
        comm->os<<"rthread wait for incoming lock"<<ps->is_thread_idle(thread_id)<<"\n";
        inList_not_empty.wait(lock);
        comm->os<<"rthread get condition" <<thread_wait<<" "<<ps->Exit()<<"\n";
    }
    
    if(ps->Exit() || ps->has_new_chunk()) {  
        return -1;
    }
    
    comm->os<<"rthread after wait\n";
    if(primary && inList.primary_size() > 0) {
        rays_stream = inList.get_primary();
    } else if (!primary && inList.secondary_size() > 0) {
        rays_stream = inList.get_secondary();
    } else {
        return rays_size;
    }

    lock.unlock();

    ps->set_thread_idle(thread_id, false);
    int copy_size = rays_stream->size;
    int width = rays_stream->width;
    printf("copy primary size %d\n", copy_size);
  //  comm->os<<"rthread width "<<width <<" logic width "<<rays_stream->logic_width<<"\n";
    memcpy(*rays, rays_stream->get_data(), ps->get_buffer_capacity() * width * sizeof(float));
    rays_stream->size = 0;

     //printf("%d rays_stream size %d %d %d\n", comm->rank, rays_stream->size, primary, rays_stream->store_width);
    int* ids = (int*)(*rays);
    comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
    for(int i = 0; i < copy_size; i ++) {
        comm->os<<"#" <<ids[i + rays_size]<<" "<<ids[i + rays_size + ps->get_buffer_capacity() * 9];
    }
    comm->os<<"\n";
    delete rays_stream;
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
    printf("work thread region %d %d %d %d\n", region[0], region[1], region[2], region[3]);
    
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
        int rnd = comm->rank * devNum + devId;
        render(&settings, rnd, devId, ps->get_local_chunk(), generate_rays);
    }
    printf("work thread end\n");
}

// get chunk retun chunk id
int Node::get_sent_list() {
    int chunk_size = ps->get_chunk_size();
    int dst_loaded = -1, dst_unloaded = -1;
    int max_loaded = 0,  max_unloaded = 0; 
    for(int i = 0; i < chunk_size; i++) {
        if(i == ps->get_local_chunk()) continue;

        if(ps->get_proc(i) >= 0 && rayList[i].size() > max_loaded) {
            max_loaded = rayList[i].size();
            dst_loaded = i;
        }
        if(ps->get_proc(i) < 0 && rayList[i].size() > max_unloaded) {
            max_unloaded = rayList[i].size();
            dst_unloaded = i;
        }
    }
    if(comm->isMaster()) {
        return dst_loaded;
    } else {
        if(dst_loaded >= 0 ) {
            if(ps->all_thread_waiting() || max_loaded > 0.4 * ps->get_buffer_size()) 
                 return dst_loaded;
            return -1;
        } 
        else 
           return dst_unloaded; 
    } 
}
