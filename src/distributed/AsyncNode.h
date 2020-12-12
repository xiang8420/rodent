// Original Master-Worker mode 
struct AsyncNode : public Node {

    AsyncNode(struct Communicator *comm, struct ProcStatus *ps);
    
    ~AsyncNode();

    int get_sent_list(); 
    
    void send_message(); 

    void save_outgoing_buffer(float *, size_t, bool); 

    void clear_retired_rays(); 
    
    void copy_to_inlist(int current_chunk);

    void save_out_list(float *, size_t, bool, int);

    RayMsg* export_ray_msg(int, int, int, bool, int); 

    static void message_thread(void* tmp);
    
    int load_incoming_buffer(float **, size_t, bool, int, bool); 
    
    void run(Scheduler * camera);

    size_t outList_size(int i); 

    bool outList_empty(); 

    bool allList_empty();

    RayStreamList  inList;  
    RayStreamList* outStreamList;
    
    std::vector<RetiredRays*> retiredRays; 
    std::mutex  retired_mutex;
    std::condition_variable retired_cond_full; 

    int chunk_size;
};

size_t AsyncNode::outList_size(int i) {
    return outStreamList[i].get_head_primary_size() 
         + outStreamList[i].get_head_secondary_size();
}

bool AsyncNode::outList_empty() {
   // std::lock_guard <std::mutex> lock(out_mutex); 
    for(int i = 0; i < chunk_size; i++)
        if(!outStreamList[i].empty() && !ps->is_local_chunk(i)) 
            return false;
    return retiredRays.empty();
}

void sort_ray_array(float* rays, int size, int chunk_size, bool primary){
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

void AsyncNode::save_out_list(float *rays, size_t size, bool primary, int rank){
    sort_ray_array(rays, size, chunk_size, primary);

    if(OUT_BUFFER) {
        printf("save to retired rays %d\n", size); 
        RetiredRays *retired_rays = new RetiredRays(rays, size, primary);
        std::lock_guard <std::mutex> lock(retired_mutex); 
        retiredRays.emplace_back(retired_rays);
    } else {
        std::lock_guard <std::mutex> lock(outStreamList[0].mutex); 
        
        statistics.start("run => work_thread => send => get_mutex");
        RayStreamList::read_from_device_buffer(outStreamList, rays, size, primary, chunk_size, rank);
        for(int i = 0; i < chunk_size; i++) 
            printf("outstreamlist size %d p %d s %d\n", i, outStreamList[i].primary_size(), outStreamList[i].secondary_size());
        statistics.end("run => work_thread => send => get_mutex");
    }
    retired_cond_full.notify_all();
}

void AsyncNode::clear_retired_rays() {
    std::lock_guard <std::mutex> lock(out_mutex); 
    while(!retiredRays.empty()) {
        RetiredRays* rays = retiredRays.back();
        retiredRays.pop_back();
        RayStreamList::read_from_device_buffer(outStreamList, rays->data, rays->size, rays->primary, chunk_size, comm->get_rank());
        delete rays;
    }
}

    
void AsyncNode::copy_to_inlist(int current_chunk) {
    std::unique_lock <std::mutex> lock(inList.mutex); 
    if(outStreamList[current_chunk].size() > 0) {
        RayStreamList::swap(inList, outStreamList[current_chunk]); 
        outStreamList[current_chunk].clear(); 
    } 
}

RayMsg* AsyncNode::export_ray_msg(int cId, int rank, int dst, bool idle, int tag) {
    RayMsg * msg;
    if(!OUT_BUFFER)
        std::unique_lock <std::mutex> lock(outStreamList[0].mutex); 
    msg = new RayMsg(outStreamList[cId], rank, dst, cId, idle, tag); 
    return msg;
}
    
bool AsyncNode::allList_empty() {
   // std::lock_guard <std::mutex> lock(out_mutex); 
    int chunk_size = ps->get_chunk_size(); 
    for(int i = 0; i < chunk_size; i++)
        if(!outStreamList[i].empty()) 
            return false;
    return inList.empty() && retiredRays.empty();
}

AsyncNode::AsyncNode(struct Communicator *comm, struct ProcStatus *ps)
    :Node(comm, ps)
{
    printf("new AsyncNode\n");
    int store_capacity = ps->get_stream_store_capacity();
    int logic_capacity = ps->get_stream_logic_capacity();
   
    chunk_size = ps->get_chunk_size(); 
    outStreamList = new RayStreamList[chunk_size];
    for(int i = 0; i < chunk_size; i++)
        outStreamList[i].set_capacity(logic_capacity, store_capacity);
    inList.set_capacity(logic_capacity, store_capacity);
}

AsyncNode::~AsyncNode() {
    printf("delete AsyncNode\n");
    delete[] outStreamList; 
}

// get chunk retun chunk id
int AsyncNode::get_sent_list() {
    int dst_loaded = -1, dst_unloaded = -1;
    int max_loaded = 0,  max_unloaded = 0; 
   
    for(int i = 0; i < chunk_size; i++) {
        if(ps->is_local_chunk(i)) continue;
        int cur_size = outList_size(i);
        if(cur_size > max_loaded) {
            if(ps->get_proc(i) >= 0 ) {
                max_loaded = cur_size;
                dst_loaded = i;
            } else {
                max_unloaded = cur_size; 
                dst_unloaded = i;
            }
        }
    }
    if(dst_loaded >= 0 ) {
        if(ps->all_thread_waiting() || max_loaded > 0.4 * ps->get_stream_logic_capacity()) 
             return dst_loaded;
        return -1;
    } 
    else 
       return dst_unloaded; 
}

void AsyncNode::send_message() {
    //comm->os<<"mthread status inlist size "<<inList.primary_size()<<" "<<inList.secondary_size()<<" thread wait "<<ps->all_thread_waiting()<<"\n";
    //comm->os<<"all thread waitting "<<ps->all_thread_waiting() <<"inList size"<<inList.size()<<"\n";
    if (ps->all_thread_waiting() && inList.empty()) {
        if(allList_empty() ) {
            ps->set_proc_idle(comm->get_rank());
            if(ps->all_proc_idle() && ps->all_rays_received()) {
                QuitMsg quit_msg(comm->get_rank(), comm->get_tag()); 
                comm->send_message(&quit_msg, ps);
                ps->set_exit();
                inList.cond_full.notify_all();
                return;
            } else {
                if(!comm->isMaster()) {
                    statistics.start("run => message_thread => send_messagei => status");
                    StatusMsg status(comm->get_rank(), comm->get_master(), ps->get_status(), ps->get_current_chunk(), comm->get_size(), comm->get_tag()); 
                    comm->send_message(&status, ps);
                    statistics.end("run => message_thread => send_messagei => status");
                }
            }
        } else if(outList_empty()) {
            ///load new chunk
            for(int i = 0;i < chunk_size; i++)
                comm->os<<outList_size(i)<<" \n";
          //  ps->switch_current_chunk(outStreamList);
            ps->switch_current_chunk();
            int current_chunk = ps->get_current_chunk();
            comm->os<<"need another chunk or need send rays "<< current_chunk << "\n"; ;
            inList.cond_full.notify_all();
            return;
        } else {
            comm->os<<"still has rays to send\n";
        }
    } else {
        ps->set_proc_busy(comm->get_rank());
    }

    int cId = get_sent_list();
    if (cId >= 0) {
        comm->os<<"mthread outlist empty "<<cId<<" "<<outList_size(cId)<<"\n";
        if(cId >= 0) {
            int dst_proc = ps->get_proc(cId);
            comm->os<<"chunk "<<cId<<" get proc "<<dst_proc<<"\n";
            if(dst_proc >= 0) {
                ps->set_proc_busy(dst_proc);
                statistics.start("run => message_thread => send_message => new RayMsg");
                RayMsg *ray_msg = export_ray_msg(cId, comm->get_rank(), dst_proc, false, comm->get_tag()); 
                statistics.end("run => message_thread => send_message => new RayMsg");
                comm->send_message(ray_msg, ps);
                delete ray_msg;
            } else {
                error("ray msg");
            }
        }
    }
}

void AsyncNode::save_outgoing_buffer(float *rays, size_t size, bool primary) {
    save_out_list(rays, size, primary, comm->get_rank());
}

int AsyncNode::load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) {
    printf("load incoming buffer\n");
    comm->os<<"rthread load incoming buffer inlist size "<<inList.size()<<"\n";

    if(ps->has_new_chunk()) 
        return -2 - ps->get_current_chunk(); 

    if(ps->Exit() ||  ps->get_chunk_size() == 1) {
        printf(" recv exit \n");
        return -1;
    }
    statistics.start("run => work_thread => load_incoming_buffer-wait");
    statistics.start("run => work_thread => load_incoming_buffer-wait => new chunk ");
    if(!inList.empty() && ps->Exit())
        warn("inlist not empty but prepared to exit\n");

    std::unique_lock <std::mutex> lock(inList.mutex); 
    ps->set_thread_idle(thread_id, true);
    retired_cond_full.notify_all();
    
    std::cout<<"rthread  inlist priamry " <<inList.primary_size()<<" secondary "<<inList.secondary_size()<<"\n";
    while (inList.empty() && !ps->Exit() && !ps->has_new_chunk()) {
        comm->os<<"rthread wait for incoming lock"<<ps->is_thread_idle(thread_id)<<"\n";
        printf("rthread wait for incoming lock\n");
        inList.cond_full.wait(lock);
        comm->os<<"rthread get condition  "<<ps->Exit()<<"\n";
        printf("rthread get condition\n");
    }
    statistics.end("run => work_thread => load_incoming_buffer-wait");
    
    if(ps->has_new_chunk()) { 
        statistics.end("run => work_thread => load_incoming_buffer-wait => new chunk ");
        assert(inList.empty()); 
        return -2 - ps->get_current_chunk(); 
    }

    if(ps->Exit()) { 
        assert(inList.empty()); 
        return -1; 
    }
    
    statistics.end("run => work_thread => load_incoming_buffer-wait");

    struct RaysStream *rays_stream;
    if(primary && inList.primary_size() > 0) {
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
    int* ids = (int*)(*rays);
    //int psize = std::min(5, copy_size);
    int psize = copy_size;
    
    
//    if(primary) {
//        comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
//        for(int i = 0; i < psize; i ++) {
//        //    comm->os<<"#" <<ids[i + rays_size]<<" "<<ids[i + rays_size + ps->get_stream_store_capacity() * 9]
//        //                                      <<" "<<ids[i + rays_size + ps->get_stream_store_capacity() * 10];
//            if((0xFF & ids[i + rays_size + ps->get_stream_store_capacity() * 10]) > 4 )
//                comm->os<<"pri pri chk "<<(0xFF & ids[i + rays_size + ps->get_stream_store_capacity() * 10])<<" read error\n";
//        }
//        comm->os<<"\n";
//    } else {
//        comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
//        for(int i = 0; i < psize; i ++) {
//         //   comm->os<<"#" <<ids[i + rays_size]<<" "<<ids[i + rays_size + ps->get_stream_store_capacity() * 9]
//         //                                     <<" "<<ids[i + rays_size + ps->get_stream_store_capacity() * 10]
//         //                                     <<" "<<ids[i + rays_size + ps->get_stream_store_capacity() * 11];
//            if((0xFF & ids[i + rays_size + ps->get_stream_store_capacity() * 11]) > 4 )
//                comm->os<<"sec pri chk "<<(0xFF & ids[i + rays_size + ps->get_stream_store_capacity() * 11])<<" read error\n";
//        }
//        comm->os<<"\n";
//    }
    delete rays_stream;

    statistics.end("run => work_thread => load_incoming_buffer-copy");
    return copy_size + rays_size;
        
}

void AsyncNode::message_thread(void* tmp) {
  
    AsyncNode *wk = (struct AsyncNode*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus * ps = wk->ps;

    statistics.start("run => message_thread");
    printf("%d message thread\n", comm->get_rank());
    int recv_count = 0;
    while (!ps->Exit()) {
        statistics.start("run => message_thread => loop");

        statistics.start("run => message_thread => recv_message");

        if(ps->is_proc_idle() && !ps->Exit())
            comm->recv_message(ps, wk->outStreamList, wk->inList, true);
        else
            comm->recv_message(ps, wk->outStreamList, wk->inList, false);
        
        statistics.end("run => message_thread => recv_message");
         
        if(ps->Exit()) {
            break;
        }

        wk->loop_check(1);
        statistics.start("run => message_thread => inlist not empty");
        
        if(!wk->inList.empty()) 
            wk->inList.cond_full.notify_one();

        statistics.end("run => message_thread => inlist not empty");
        wk->loop_check(3);

        if(ps->Exit()) {
            break;
        }
        
        wk->loop_check(3.5);
        statistics.start("run => message_thread => send_message");
        //clear_outlist();
        
        wk->clear_retired_rays();
        
        wk->send_message();
        
        wk->loop_check(4);
        comm->purge_completed_mpi_buffers();
        statistics.end("run => message_thread => send_message");
        statistics.start("run => message_thread => new_chunk");
        {
            if(ps->has_new_chunk()) {
                wk->inList.cond_full.notify_all();
                std::unique_lock <std::mutex> lock(wk->thread_mutex); 
                while( ps->has_new_chunk()) {
                    wk->render_start.wait(lock); 
                }
            }
        }
        statistics.end("run => message_thread => new_chunk");

        statistics.start("run => message_thread => sleep");
        usleep(500);
        statistics.end("run => message_thread => sleep");
        statistics.end("run => message_thread => loop");
	}
    wk->inList.cond_full.notify_all();
    statistics.end("run => message_thread => loop");
    statistics.end("run => message_thread");
    comm->os <<" end message thread"<<ps->all_thread_waiting()<<"\n";
    comm->os <<" inlist "<< wk->inList.size()
             <<" recv "<<ps->global_rays[comm->get_rank() + comm->get_size()]
             <<std::endl;
    return;
} 

void AsyncNode::run(Scheduler * camera) {
    
    comm->os <<" start run message thread \n";

    int deviceNum = ps->get_dev_num();
    int iter = 0; 
   
    std::thread mthread(message_thread, this);
    do {
        //load new chunk;
        int current_chunk = ps->get_current_chunk();
        comm->os<<"mthread copy new chunk " << current_chunk << "\n";
        copy_to_inlist(current_chunk);
        statistics.start("run => render thread => get start ");
        launch_rodent_render(camera, deviceNum, iter==0);

        iter++;
    } while(!ps->Exit());
    mthread.join();
    return;
}

