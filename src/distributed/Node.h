class Node {

protected:
    std::vector <RayList> rayList;
    RayStreamList inList;  
    
    Communicator *comm;
    ProcStatus *ps;
     
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

    void save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary);

    int load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait);
    
    static void work_thread(void* tmp, ImageDecomposition * camera, int devId, int devNum, bool preprocess, bool shoot_rays);

    virtual void run(ImageDecomposition * camera) = 0;

    int get_sent_list(); 
    
    bool rayList_empty();
    
    bool outList_empty();
    
    bool inout_list_empty();

    int get_tag();
    
    bool all_queue_empty();

    void clear_outlist();
};

Node::Node(Communicator *comm, ProcStatus *ps) : comm(comm), ps(ps) {
    msg_tag = 0;
    sync = false;
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++) {
        rayList.emplace_back(RayList());
        //rayList[i].set_capacity(1048608, "out");
    }
    inList.set_capacity(ps->get_buffer_size());
  //  for(int i = 0; i < ps->get_chunk_size(); i++) { 
  //      printf("new out going ray data size %ld capacity %ld\n",rayList[i].secondary.data.size(), rayList[i].secondary.data.capacity());
  //      printf("new out going ray data size %ld capacity %ld\n",rayList[i].primary.data.size(), rayList[i].primary.data.capacity());
  //  }
}

int Node::get_tag(){
    return msg_tag++;
}

Node::~Node() {
    printf("%d delete Node\n", comm->get_rank());
}

bool Node::all_queue_empty(){ 
    for(int i = 0; i < ps->get_chunk_size(); i++) {
        if(!rayList[i].empty()) {return false;}
    }
    return true;
}

bool Node::outList_empty() {  
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++) {
        if( !rayList[i].empty() ) {
            return false;
        }
    }
    return true;
}

void Node::clear_outlist() {  
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++) {
        rayList[i].get_primary().clear();
        rayList[i].get_secondary().clear();
    }
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

void Node::save_outgoing_buffer(float *retired_rays, size_t size, size_t capacity, bool primary){
    out_mutex.lock(); 
    int width = primary?21:14; 
    comm->os<<"rthread save outgoing buffer "<< size <<" width "<<width<<"\n"; 
    int* ids = (int*)(retired_rays);
    for(int i = 0; i < 5; i ++) {
        comm->os<<"| "<< ids[i * width] <<" "<< ids[i * width + 9] << " ";
    }
    comm->os<<"\n";
//    for(int i = 0; i < ps->get_chunk_size(); i++) { 
//        printf("save out going ray data size %ld capacity %ld\n",rayList[i].secondary.data.size(), rayList[i].secondary.data.capacity());
//        printf("save out going ray data size %ld capacity %ld\n",rayList[i].primary.data.size(), rayList[i].primary.data.capacity());
//    }
    RayList::read_from_device_buffer(rayList.data(), retired_rays, size, capacity, primary, comm->get_rank(), ps->get_chunk_size());
    out_mutex.unlock(); 
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
    if(!inList.empty() && ps->Exit())
        error("inlist not empty but prepared to exit\n");

    std::unique_lock <std::mutex> lock(inList.mutex); 
    if(!sync) {
        ps->set_thread_idle(thread_id, thread_wait);
        
        std::cout<<"rthread "<<thread_wait<< " inlist priamry " <<inList.primary_size()<<" secondary "<<inList.secondary_size()<<"\n";
        while (thread_wait && inList.empty() && !ps->Exit() && !ps->has_new_chunk()) {
            comm->os<<"rthread wait for incoming lock"<<ps->is_thread_idle(thread_id)<<"\n";
            inList_not_empty.wait(lock);
            comm->os<<"rthread get condition" <<thread_wait<<" "<<ps->Exit()<<"\n";
        }
        if(ps->Exit() || ps->has_new_chunk()) {  
            return -1;
        }
    }

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

    ps->set_thread_idle(thread_id, false);
    int copy_size = rays_stream->size;
    int width = rays_stream->width;
    printf("copy primary size %d\n", copy_size);
 //   comm->os<<"rthread width "<<width <<" logic width "<<rays_stream->logic_width<<"\n";
    memcpy(*rays, rays_stream->get_data(), ps->get_buffer_capacity() * width * sizeof(float)); 
    rays_stream->size = 0;
 
     //printf("%d rays_stream size %d %d %d\n", comm->get_rank(), rays_stream->size, primary, rays_stream->store_width);
    int* ids = (int*)(*rays);
    comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
    for(int i = 0; i < std::min(5, copy_size); i ++) {
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
    
    
    Settings settings {
        Vec3 { camera->eye.x, camera->eye.y, camera->eye.z },
        Vec3 { camera->dir.x, camera->dir.y, camera->dir.z },
        Vec3 { camera->up.x, camera->up.y, camera->up.z },
        Vec3 { camera->right.x, camera->right.y, camera->right.z },
        camera->w, camera->h,
        Vec4_i32 { region[0], region[1], region[2], region[3]},
        sppDev,
        comm->get_size()
    };
    if(preRendering) {
        prerender(&settings);
    } else {
        int rnd = comm->get_rank() * devNum + devId;
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
