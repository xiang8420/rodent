
void swap_array(float *rays, int k, int j, int width) {
    int *tmp = new int[width];
    memcpy(tmp, &rays[width * k], width * sizeof(float));
    memcpy(&rays[width * k], &rays[width * j], width * sizeof(float));
    memcpy(&rays[width * j], tmp, width * sizeof(float));
    delete[] tmp; 
}
    
struct RayListManager {
    RayStreamList  inList;  
    RayArrayList * outArrayList;
    RayStreamList* outStreamList;
    bool stream;
    int chunk_size;
    std::mutex  out_mutex;

    RayListManager(bool stream, int chunk_size, int logic_capacity, int store_capacity)
        : stream(stream), chunk_size(chunk_size) 
    {
        if(stream) {
            outStreamList = new RayStreamList[chunk_size];
            for(int i = 0; i < chunk_size; i++)
                outStreamList[i].set_capacity(logic_capacity, store_capacity);
            outArrayList = NULL;
        } else {
            outArrayList = new RayArrayList[chunk_size];
            outStreamList = NULL;
        }
        inList.set_capacity(logic_capacity, store_capacity);
    }
    
    ~RayListManager(){
        if(outArrayList != NULL) delete[] outArrayList;
        if(outStreamList != NULL) delete[] outStreamList; 
    } 

    bool inList_empty() { return inList.empty(); }
    
    size_t inList_size() { return inList.size(); }

    bool outList_empty(ProcStatus *ps) {
        std::lock_guard <std::mutex> lock(out_mutex); 
        if(stream) {
            for(int i = 0; i < chunk_size; i++)
                if(!outStreamList[i].empty() && !ps->is_local_chunk(i)) 
                    return false;
            return true;
        } else {
            for(int i = 0; i < chunk_size; i++)
                if(!outArrayList[i].empty() && !ps->is_local_chunk(i)) 
                    return false;
            return true;
        }

    }
    
    void clear_outList() {
        if(stream) { 
            for(int i = 0; i < chunk_size; i++) {
                outStreamList[i].clear();
            }
        } else {
            for(int i = 0; i < chunk_size; i++) {
                outArrayList[i].get_primary().clear();
                outArrayList[i].get_secondary().clear();
            }
        }
    }

    bool allList_empty() {
        std::lock_guard <std::mutex> lock(out_mutex); 
        if(stream) {
            for(int i = 0; i < chunk_size; i++)
                if(!outStreamList[i].empty()) 
                    return false;
        } else {
            for(int i = 0; i < chunk_size; i++)
                if(!outArrayList[i].empty()) 
                    return false;
        }
        return inList_empty();
    }

    size_t outList_size(int i) {
        if(stream)
            return outStreamList[i].get_head_primary_size() 
                 + outStreamList[i].get_head_secondary_size();
        else
            return outArrayList[i].size();
    }

    RayArrayList * get_outArrayList() { return outArrayList; }

    RayStreamList * get_outStreamList() { return outStreamList; }

    RayStreamList& get_inList() { return inList; }

    void sort_ray_array(float* rays, int size, bool primary){
        int width = primary ? PRIMARY_WIDTH : SECONDARY_WIDTH;
        std::vector<int> end;
        std::vector<int> begin;

        for(int i = 0; i < chunk_size; i++) { 
            begin.emplace_back(0);
            end.emplace_back(0);
        }
        int * iptr = (int*) rays;
        for(int i = 0; i < size; i++) { 
            int chunk = iptr[i * width + 9];
            end[chunk]++;
        }
        
        int n = 0;
        for(int i = 0; i < chunk_size; i++) { 
            begin[i] = n;
            n += end[i];
            end[i] = n;
        }

        for(int i = 0; i < chunk_size; i++) {
            int st = begin[i];
            int ed = end[i];
            int j = st;
            while ( j < ed) {
                int cid = iptr[j * width + 9]; 
                if (cid != i) {
                    int k = begin[cid]++;
                    swap_array(rays, k, j, width);                
                } else {
                    j++;
                }
            }
        }
        
    }
    
    void save_out_list(float *rays, size_t size, bool primary, int rank){
        statistics.start("run => work_thread => send => get_mutex");
        
        printf("save rays %d \n", size);
        sort_ray_array(rays, size, primary);
        std::lock_guard <std::mutex> lock(out_mutex); 
        if(stream) {
            RayStreamList::read_from_device_buffer(outStreamList, rays, size, primary, chunk_size, rank);
            for(int i = 0; i < chunk_size; i++) 
                printf("outstreamlist size %d p %d s %d\n", i, outStreamList[i].primary_size(), outStreamList[i].secondary_size());
        } else {
            int width = primary?PRIMARY_WIDTH : SECONDARY_WIDTH; 
            RayArrayList::read_from_device_buffer(outArrayList, rays, size, primary, chunk_size);
        }
        printf("save rays %d \n", size);
        
        statistics.end("run => work_thread => send => get_mutex");
    }

    RayMsg* export_ray_msg(int cId, int rank, int dst, bool idle, int tag) {
        RayMsg * msg;
        std::lock_guard <std::mutex> lock(out_mutex); 
        if(stream)      
            msg = new RayMsg(outStreamList[cId], rank, dst, cId, idle, tag); 
        else 
            msg = new RayMsg(outArrayList[cId], rank, dst, cId, idle, tag); 
        return msg;
    }
    
    void copy_to_inlist(int current_chunk, int rank) {
        std::unique_lock <std::mutex> lock(inList.mutex); 
        if(stream && outStreamList[current_chunk].size() > 0) {
            RayStreamList::swap(inList, outStreamList[current_chunk]); 
            outStreamList[current_chunk].clear(); 
        } 
        if(!stream && outArrayList[current_chunk].size() > 0) {
            RaysArray &primary   = outArrayList[current_chunk].get_primary(); 
            RaysArray &secondary = outArrayList[current_chunk].get_secondary(); 
            inList.read_from_ptr(primary.get_data(), 0, primary.get_size(), true, rank);
            inList.read_from_ptr(secondary.get_data(), 0, secondary.get_size(), false, rank);
            outArrayList[current_chunk].clear();
            //read to inList
        //    clear_outlist();
        }
    }

};

class Node {

protected:
    Communicator *comm;
    ProcStatus *ps;
    RayListManager *rlm; 

    std::mutex  thread_mutex;
    
    std::condition_variable inList_not_empty; //
    std::condition_variable render_start; //

    int msg_tag;
public: 
    int max_used_mem, pre_mem;
    bool sync, stream;

    Node(Communicator *comm, ProcStatus *ps);

    ~Node();

    ProcStatus * get_proc_status(){return ps;}

    Communicator * get_communicator(){ return comm; }

    int load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait);
    
    static void work_thread(void* tmp, ImageDecomposition * camera, int devId, int devNum, bool preprocess, bool shoot_rays);

    virtual void run(ImageDecomposition * camera) = 0;

    void save_outgoing_buffer(float *, size_t, bool); 
    
    int get_tag();
    
    void loop_check(float i); 
};

Node::Node(Communicator *comm, ProcStatus *ps) : comm(comm), ps(ps) {
    stream = false;
    int chunk_size = ps->get_chunk_size();
    int store_capacity = ps->get_stream_store_capacity();
    int logic_capacity = ps->get_stream_logic_capacity();
    rlm = new RayListManager(stream,  chunk_size, logic_capacity, store_capacity);

    msg_tag = 0;
    pre_mem = 0;
    sync = false;
    
    comm->outStreamList = rlm->get_outStreamList();
    comm->outArrayList = rlm->get_outArrayList();
    comm->inList = &(rlm -> get_inList());

    max_used_mem = 0;
}

int Node::get_tag(){
    return msg_tag++;
}

Node::~Node() {
    delete rlm;
    printf("%d delete Node\n", comm->get_rank());
}

void Node::save_outgoing_buffer(float *rays, size_t size, bool primary) {
    rlm->save_out_list(rays, size, primary, comm->get_rank());
}

int Node::load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) {
    printf("load incoming buffer\n");
    RayStreamList& inList = rlm->get_inList();
    comm->os<<"rthread load incoming buffer inlist size "<<inList.size()<<"\n";
    if(inList.empty()) {
        if(sync)
            return -1;
        if(ps->Exit() || ps->has_new_chunk() || ps->get_chunk_size() == 1) {
            comm->os << "rthread exit new chunk exit : "<< ps->Exit() << " new chunk " << ps->has_new_chunk() <<"\n";
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
            printf("rthread wait for incoming lock\n");
            inList.cond_full.wait(lock);
            comm->os<<"rthread get condition  "<<ps->Exit()<<"\n";
            printf("rthread get condition\n");
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
        if(inList.secondary_size() > inList.primary_size())
            return rays_size;
        rays_stream = inList.get_primary();
    } else if (!primary && inList.secondary_size() > 0) {
        rays_stream = inList.get_secondary();
    } else {
        return rays_size;
    }
    inList.empty_notify(); //tell mthread inlist size changed 
    lock.unlock();

    statistics.start("run => work_thread => load_incoming_buffer-copy");

    ps->set_thread_idle(thread_id, false);
    int copy_size = rays_stream->size;
    comm->os<<"all rays chunk "<<ps->get_current_chunk()<<" size "<<copy_size<<"\n";
    int width = rays_stream->width;
    printf("copy primary size %d\n", copy_size);
 //   comm->os<<"rthread width "<<width <<" logic width "<<rays_stream->logic_width<<"\n";
    memcpy(*rays, rays_stream->get_data(), ps->get_stream_store_capacity() * width * sizeof(float)); 
 
//     //printf("%d rays_stream size %d %d %d\n", comm->get_rank(), rays_stream->size, primary, rays_stream->store_width);
    delete rays_stream;

    statistics.end("run => work_thread => load_incoming_buffer-copy");
    return copy_size + rays_size;
        
}

void Node::work_thread(void* tmp, ImageDecomposition * splitter, int devId, int devNum, bool preRendering, bool generate_rays) {
    printf("work threadstart\n");

    Node * wk = (Node*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus *ps = wk->ps;
    
    Camera *camera = splitter->camera;
    int* region = splitter->get_render_block();
    int sppProc = splitter->get_spp(); 
    int sppDev = sppProc / devNum;
   
    comm->os<<"Load new chunk start rendering\n"; 
    printf("width %d height%d spp %d dev id %d local chunk %d\n", camera->width, camera->height, sppDev, devId, ps->get_current_chunk() );
    printf("work thread rank %d region %d %d %d %d local chunk %d\n", comm->get_rank(), region[0], region[1], region[2], region[3], ps->get_current_chunk());
    comm->os<<"start render thread chunk "<<ps->get_current_chunk()<<"\n";
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
        int rnd = camera->iter * (comm->get_rank() + 1) * devNum + devId;
        render(&settings, rnd, devId, ps->get_current_chunk(), generate_rays);
    }
    printf("work thread end\n");
    statistics.end("run => work_thread");
}

void Node::loop_check(float i) {
    if(1) {    
        int mem = physical_memory_used_by_process();
        if(/*i> 1000 ||*/mem != pre_mem) {
            comm->os<<i<<" memory use "<<mem<<"\n"; 
            max_used_mem = std::max(mem, max_used_mem);
            pre_mem = mem;
        }
        //printf("%d mark %f\n", comm->get_rank(), i);
    }
}

