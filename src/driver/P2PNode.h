
int32_t * get_prerender_result();

// Asynchronous p2p model without master

class P2PNode : public Node {

protected:    
    int * statistic;
    struct RayList **List;
    struct RayList **outList;
    struct RayList * inList;
    struct RayList * buffer;
    
    std::mutex  out_mutex, buffer_mutex, in_mutex;
   
    int renderer_get_rays, renderer_save_rays, recv_rays, write_rays, sent_rays; 
public:    
    P2PNode(struct Communicator *comm, struct ProcStatus *ps);

    ~P2PNode();

    bool incoming_empty(){ return inList->empty() && buffer->empty();}

    bool all_queue_empty() { 
        bool res = true;
        out_mutex.lock();
        for(int i = 0;i < worker_size; i++)
            if(!outList[i]->empty()) { res = false; break; }

        out_mutex.unlock();
        in_mutex.lock();
        res &= inList->empty() && buffer->empty();
        in_mutex.unlock();
        return res;
    }
    
    bool check_rendering_status(); 
   
    void set_distributed_buffer(); 
    // send primary and secondary
    void save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary);
  
    void write_rays_buffer();
    
    int load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait);
   
    void mpi_thread();
    
    void count_rays();

    static void message_thread(void* tmp);
    
    static void work_thread(void* tmp, float* processTime, int devId, int devNum, bool preprocess);

    void run(float* processTime);

}; 


void P2PNode::set_distributed_buffer() {
    int chunk_size = ps->get_chunk_size();
    List = new RayList *[chunk_size];
    for(int i = 0; i < chunk_size; i++) {
        List[i] = new RayList(ps->get_buffer_size() * 4, "out", false);
    }
    for(auto c : ps->get_local_chunk()) {
        List[c]->type = "in";
    }
    outList = List;
    comm->os << "out list size"<<outList[0]->size() << "capacity" << outList[0]->get_primary()->get_capacity() << std::endl; 
    comm->os << "loaded chunk()"<<ps->get_loaded_chunk()<< std::endl; 
    //when we need to load a new chunk updata inList pointer.
    inList  = List[ps->get_loaded_chunk()];
    buffer  = new RayList(ps->get_buffer_size(), "buffer", false);

    ps->set_chunks();
}

P2PNode::P2PNode(struct Communicator *comm, struct ProcStatus *ps) 
    : Node(comm, ps) 
{
    int chunk_size = ps->get_chunk_size();
    statistic = new int[chunk_size * 4/*dep*/ * 2]; 
    std::fill(statistic, statistic + chunk_size * 4/*dep*/ * 2, 0);
    renderer_get_rays = 0;
    renderer_save_rays = 0;
    write_rays = 0; 
    sent_rays = 0;
}

P2PNode::~P2PNode() {
    delete buffer;
    for(int i = 0; i < ps->get_chunk_size(); i++){
        delete List[i];
    }
}

void P2PNode::write_rays_buffer() {
    comm->os<<"mthread write ray buffer"<<inList->size()<<"\n";
//    comm->os<<"mthread buffer size "<<buffer->size()<<"\n";
    if(inList->empty()) return;
    
    std::unique_lock <std::mutex> lock(buffer->mutex); 
    while(!buffer->empty()) {
        buffer_not_full.wait(lock);
    }
    comm->os<<"mthread get inlist lock\n";
    
    inList->write_to_device_buffer(buffer, comm->rank);
    
    buffer_not_empty.notify_one();
    lock.unlock();
    comm->os<<"mthread release inlist lock\n";
}

int P2PNode::load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) {
    if(comm->size == 1 || ps->Exit()) {
        comm->os << "rthread exit\n";
        return -1;
    }
    
    comm->os<<"rthread buffer size" <<buffer->size()<<"\n";
    if(buffer->empty()) { 
        comm->os<<"rthread notify" <<buffer->size()<<"\n";
        buffer_not_full.notify_all();
        comm->os<<"rthread notify" <<buffer->size()<<"\n";
    }

    struct Rays *queue = primary ? buffer->get_primary() : buffer->get_secondary();
    comm->os <<"rthread "<<thread_id<<"read incoming buffer"<<thread_wait<< "size "<<queue->size<<"\n";
    int width = primary ? 21 : 14;
    std::unique_lock <std::mutex> lock(buffer->mutex); 
    
    ps->set_thread_idle(thread_id, thread_wait);
    comm->os <<"rthread idle "<< ps->is_thread_idle(thread_id) 
             <<" width "<<width
             <<" queue->size"<< queue->size
             <<" exit"<<ps->Exit() 
             <<" buffersize"<<buffer->size()
             <<" thread wait"<<thread_wait
             << std::endl;
    while (thread_wait && buffer->empty() && !ps->Exit()) {
        comm->os<<"rthread wait for incoming lock"<<ps->is_thread_idle(thread_id)<<"\n";
        buffer_not_empty.wait(lock);
        comm->os<<"rthread get not empty condition" <<thread_wait<<" "<<buffer->empty()<<ps->Exit()<<"\n";
    }

    if(!queue->empty()) {
        ps->set_thread_idle(thread_id, false);
        int copy_size = queue->size; 
        memcpy(*rays, queue->get_data(), ps->get_buffer_capacity() * width * sizeof(float)); 
        queue->size = 0;

        //printf("%d queue size %d %d %d\n", comm->rank, queue->size, primary, queue->store_width);
        int* ids = (int*)(*rays);
        comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
        for(int i = 0; i < 10; i ++) {
            comm->os<<"#" <<ids[i + rays_size]<<" "<<ids[i + rays_size + ps->get_buffer_capacity() * 9];
        }
        comm->os<<"\n";
        renderer_get_rays += copy_size;
        return copy_size + rays_size;
    }
    if(buffer->empty()) { 
        comm->os<<"rthread notify" <<buffer->size()<<"\n";
        buffer_not_full.notify_all();
        comm->os<<"rthread notify" <<buffer->size()<<"\n";
        
    }
    lock.unlock(); 

    if(ps->Exit()) {  
//        printf("ray size %ld   queue primary  size %d queue secondary size %d\n", 
//                rays_size, inList->get_primary()->size, inList->get_secondary()->size);
//        printf("recv stop mpi %d thread %d %ld\n", comm->rank, thread_id, rays_size);
        return -1;
    }
    return rays_size;
}

void P2PNode::save_outgoing_buffer(float *retired_rays, size_t size, size_t capacity, bool primary){
    comm->os<<"rthread save outgoing buffer"<<size<<"\n";
    renderer_save_rays += size;
    out_mutex.lock(); 
    int width = primary?21:14; 
    int* ids = (int*)(retired_rays);
    for(int i = 0; i < 5; i ++) {
        comm->os<<"| "<< ids[i] <<" "<< ids[i + ps->get_buffer_capacity() * 9] << " ";
    }
    comm->os<<"\n";
    RayList::read_from_device_buffer(outList, retired_rays, size, capacity, primary, comm->rank);
    out_mutex.unlock(); 
}

void P2PNode::work_thread(void* tmp, float *process_time, int devId, int devNum, bool preRendering) {
  
    P2PNode * wk = (struct P2PNode*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus *ps = wk->ps;
   
    int region[4]; 
    int sppProc = ps->get_tile_info(&region[0], process_time); 
    int sppDev = sppProc / devNum;
    int seed = comm->rank * devNum + devId;
    printf("width %d height%d spp %d dev id %d\n", ps->width, ps->height, sppDev, devId);
    Settings settings {
        Vec3 { ps->eye.x, ps->eye.y, ps->eye.z },
        Vec3 { ps->dir.x, ps->dir.y, ps->dir.z },
        Vec3 { ps->up.x, ps->up.y, ps->up.z },
        Vec3 { ps->right.x, ps->right.y, ps->right.z },
        ps->w, ps->h,
        Vec4_i32 { region[0], region[1], region[2], region[3]},
        sppDev
    };
    if(preRendering) {
        prerender(&settings);
    } else {
        render(&settings, seed, devId, wk->ps->get_loaded_chunk(), true);
    }
}

int get_sent_list(RayList ** raylist, ProcStatus *ps) {
    int n = ps->get_chunk_size();
    int t = -1; 
    int max = 0; 
    for(int i = 0; i < n; i++) {
        if(raylist[i] -> size() > max && raylist[i]->type == "out") {
            max = raylist[i] -> size();
            t = i;
        }
    }
    if(max <= 0 || t < 0) 
        return -1;
    
    return t;
}

//check proc status, return if proc need to wait 
bool P2PNode::check_rendering_status() {
    if(ps->all_thread_waiting() && buffer->empty() ) {
        if(all_queue_empty()) {
            comm->os<< "mthread if all thread idle all queue empty, set itself idle\n";
            ps->set_self_idle();
            if(ps->all_proc_idle() && ps->all_rays_received()) {
                int *s = ps->get_status();
                comm->os<< "exit all rays received "<< s[0] <<" "<< s[1] <<" "<<s[2]<<" "<<s[3]<<"\n";

                //if all worker end synchronous 
                comm->os<<"mthread all proc idle send collective\n";
                QuitMsg quit_msg(comm->rank); 
                comm->send_message(&quit_msg, ps);
                ps->set_exit();
                comm->os<<"set exit\n";
                buffer_not_empty.notify_all();
            } else {
                int *s = ps->get_status();
                comm->os<< "all rays received "<< s[0] <<" "<< s[1] <<" "<<s[2]<<" "<<s[3]<<"\n";
                comm->os<<"mthread proc idle send status\n";
                StatusMsg status(comm->rank, ps->get_status(), comm->size); 
                comm->send_message(&status, ps);
            }
            return true; 
        } else {
            //maybe load new chunk
            comm->os<<"need another chunk or need send rays\n";
            ps->set_proc_busy(comm->rank);
            return false; 
        }
    } else {
//        comm->os<<"mthread proc busy\n";
        ps->set_proc_busy(comm->rank);
        return false;
    }
} 

void P2PNode::message_thread(void* tmp) {
  
    P2PNode *wk = (struct P2PNode*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus * ps = wk->ps;

    comm->os<<"message thread\n";
    int recv_count = 0;
    while (!ps->Exit()) {
       
        wk->check_rendering_status(); 
        
        if(ps->Exit()) break;
        bool recv = false;
        do {
            recv = comm->recv_message(wk->List, ps);
        } while(!recv && ps->is_proc_idle() && !ps->Exit()) ;

        wk->write_rays_buffer();
        
        if(ps->Exit()) break;

        int cId = get_sent_list(wk->outList, ps);
//            wk->out_mutex.lock();
//        comm->os<<"mthread get send list"<<cId<<"size "<<wk->outList[cId]->size()<<"\n";
//            wk->out_mutex.unlock(); 
        
        if(cId >= 0 ) {
            ps->set_proc_busy(ps->get_dst_proc(cId));
            wk->out_mutex.lock();
            comm->os<<"mthread new RayMsg"<<cId<<"size "<<wk->outList[cId]->size()<<"\n";
            wk->sent_rays += wk->outList[cId]->size();
            RayMsg *ray_msg = new RayMsg(wk->outList[cId], comm->rank, ps->get_dst_proc(cId), cId, false); 
            comm->os<<"mthread RayMsg"<<ray_msg->get_chunk()<<"\n";
            wk->out_mutex.unlock(); 
            comm->send_message(ray_msg, ps);
        }
        comm->purge_completed_mpi_buffers();
        comm->os<<"mthread single loop end \n";
        usleep(100);
	}
    comm->os <<" end message thread"<<ps->all_thread_waiting()<<"\n";
    comm->os <<" inlist "<< wk->inList->size()
             <<" buffersize"<<wk->buffer->size()
             <<" get render rays "<<wk->renderer_get_rays
             <<" save render rays" <<wk->renderer_save_rays
             <<" sent rays" << wk->sent_rays
             <<" writer_rays " << wk->write_rays
             <<" recv "<<ps->global_rays[ps->proc_rank + ps->proc_size]
             <<std::endl;
    wk->buffer_not_empty.notify_all();
    return;
} 

void P2PNode::count_rays() {
    int width = ps->width; 
    int height = ps->height;
    int chunk_size = ps->get_chunk_size();
    int *data = get_prerender_result();
    int pixel_size = width * height;
    for(int dep = 0; dep < 8; dep++) {
        for(int i = 0; i < height; i ++) {
            for(int j = 0; j < width; j++) {
                int chunk = data[j + i * width + dep * pixel_size];
                if(chunk != -1)
                    statistic[dep * chunk_size + chunk] ++;
            }
        }
    }
    printf("statistic: ");
    for(int dep = 0; dep < 8; dep++) {
        for(int i = 0; i < chunk_size; i ++) {
            printf("%d ", statistic[dep * chunk_size + i]);
        }
        printf("\n");
    }
}
void P2PNode::run(float* frame_time) {
    int deviceNum = ps->get_dev_num();
    bool PreRendering = false;
    if(PreRendering && comm->rank == 0) {
        work_thread(this, frame_time, 0, 1, true);
        count_rays(); 
    }
    // set domain and image distribution
    set_distributed_buffer();
    
    std::vector<std::thread> workThread;
    if(comm->size == 1) {
        printf("only one worker %d  ", comm->rank);
        for(int i = 0; i < deviceNum; i++) 
            workThread.emplace_back(std::thread(work_thread, this, frame_time, i, deviceNum, false));
        
        for( auto &thread: workThread) 
            thread.join();
    	
    } else {
        std::thread mthread(message_thread, this);
        //    //all proc start proc 0 send schedule? 
        while(!ps->Exit()) {
            int chunk = ps->get_new_chunk(); 
            for(int i = 0; i < deviceNum; i++) 
                workThread.emplace_back(std::thread(work_thread, this, frame_time, i, deviceNum, false));
            
            for( auto &thread: workThread) 
                thread.join();
            
        }
        printf("%d worker exit\n", comm->rank);
        mthread.join();
    }
    return;
} 
