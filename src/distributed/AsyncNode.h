// Original Master-Worker mode 
struct AsyncNode : public Node {

    AsyncNode(struct Communicator *comm, struct ProcStatus *ps);
    
    ~AsyncNode();

    int get_sent_list(); 
    
    void send_message(); 

    void save_outgoing_buffer(float *, size_t, bool); 

    void clear_retired_rays(); 
    
    void compact_retired_rays(); 
    
    void copy_to_inlist(int current_chunk);

    RayMsg* export_ray_msg(int, int, int, bool, int); 

    static void message_thread(void* tmp);
    
    int load_incoming_buffer(float **, bool, int); 
    
    void run(Scheduler * camera);

    bool outList_empty(); 

    bool allList_empty();

    RayStreamList  inlist_comm;  
    RayStreamList* outlist_comm;
    
    RayStreamList inlist_render;  
    RayStreamList* outlist_render;

    std::vector<RetiredRays*> retiredRays; 
    std::mutex  retired_mutex;
    std::condition_variable retired_cond_full; 

    std::condition_variable thread_wakeup; 
    int chunk_size;
    bool wait_respond;
};

bool AsyncNode::outList_empty() {
   // std::lock_guard <std::mutex> lock(out_mutex); 
    for(int i = 0; i < chunk_size; i++)
        if(!outlist_comm[i].empty() && !ps->is_local_chunk(i)) 
            return false;
    return true;//retiredRays.empty();
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

    
void AsyncNode::copy_to_inlist(int current_chunk) {
    std::unique_lock <std::mutex> lock(inlist_comm.mutex); 
    if(outlist_comm[current_chunk].size() > 0) {
        RayStreamList::swap(inlist_comm, outlist_comm[current_chunk]); 
        outlist_comm[current_chunk].clear(); 
    } 
}

RayMsg* AsyncNode::export_ray_msg(int cId, int rank, int dst, bool idle, int tag) {
    RayMsg * msg;
    if(!OUT_BUFFER)
        std::unique_lock <std::mutex> lock(outlist_comm[0].mutex); 
    msg = new RayMsg(outlist_comm[cId], rank, dst, cId, idle, tag); 
    return msg;
}
    
bool AsyncNode::allList_empty() {
   // std::lock_guard <std::mutex> lock(out_mutex); 
    int chunk_size = ps->get_chunk_size(); 
    for(int i = 0; i < chunk_size; i++)
        if(!outlist_comm[i].empty()) 
            return false;
    return inlist_comm.empty();// && retiredRays.empty();
}

AsyncNode::AsyncNode(struct Communicator *comm, struct ProcStatus *ps)
    :Node(comm, ps)
{
    printf("new AsyncNode\n");
    int store_capacity = ps->get_stream_store_capacity();
    int logic_capacity = ps->get_stream_logic_capacity();
   
    chunk_size = ps->get_chunk_size(); 
    outlist_comm = new RayStreamList[chunk_size];
    for(int i = 0; i < chunk_size; i++)
        outlist_comm[i].set_capacity(logic_capacity, store_capacity);
    inlist_comm.set_capacity(logic_capacity, store_capacity);
    wait_respond = false;
}

AsyncNode::~AsyncNode() {
    printf("delete AsyncNode\n");
    delete[] outlist_comm; 
}

// get chunk retun chunk id
int AsyncNode::get_sent_list() {
    int chk = -1;
    int max = 0; 
   
    for(int i = 0; i < chunk_size; i++) {
        if(ps->is_local_chunk(i)) continue;
        int cur_size = outlist_comm[i].size();
        if(cur_size > max) {
            max = cur_size;
            chk = i;
        }
    }
    if(chk >= chunk_size || chk < 0) return -1;
    
    if(ps->is_proc_idle() || max > MIN_SEND_STREAM_SIZE) 
         return chk;
    else
        return -1;
}

void AsyncNode::send_message() {
    //comm->os<<"mthread status inlist size "<<inlist_comm.primary_size()<<" "<<inlist_comm.secondary_size()<<" thread wait "<<ps->all_thread_waiting()<<"\n";
    //comm->os<<"proc idle "<<ps->is_proc_idle() <<"inlist_comm size"<<inlist_comm.size()<<"\n";
    if (ps->is_proc_idle()) {
        if(allList_empty()) {
            if(ps->all_proc_idle() && ps->all_rays_received()) {
                comm->os<<"mthread send quit\n";
                QuitMsg quit_msg(comm->get_rank(), comm->get_tag()); 
                comm->send_message(&quit_msg, ps);
                ps->set_exit();
                return;
            } else {
                if(!comm->isMaster()) {
                    statistics.start("run => message_thread => send_messagei => status");
                    StatusMsg status(comm->get_rank(), comm->get_master(), ps->get_status(), ps->get_current_chunk(), comm->get_size(), comm->get_tag()); 
                    comm->send_message(&status, ps);
                    statistics.end("run => message_thread => send_messagei => status");
                }
                comm->os<<"mthread send status\n";
                wait_respond = true;
            }
        } else if(outList_empty()) {
            ///load new chunk
            ps->switch_current_chunk(outlist_comm);
          //  ps->switch_current_chunk();
            
            int current_chunk = ps->get_current_chunk();
            copy_to_inlist(current_chunk);
            inlist_comm.cond_full.notify_all();
            inlist_render.cond_full.notify_all();
            return;
        } else {
            comm->os<<"still has rays to send\n";
        }
    } 

    int cId = get_sent_list();
    if (cId >= 0) {
        comm->os<<"mthread outlist empty "<<cId<<" "<<outlist_comm[cId].size()<<"\n";
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
    statistics.start("run => work_thread => save_outgoing");
    sort_ray_array(rays, size, chunk_size, primary);

    if(OUT_BUFFER) {
        printf("save to retired rays %d\n", size); 
        RetiredRays *retired_rays = new RetiredRays(rays, size, primary);
        std::lock_guard <std::mutex> lock(retired_mutex); 
        retiredRays.emplace_back(retired_rays);
    } else {
        std::lock_guard <std::mutex> lock(outlist_comm[0].mutex); 
        
        statistics.start("run => work_thread => send => get_mutex");
        RayStreamList::read_from_device_buffer(outlist_comm, rays, size, primary, chunk_size, comm->get_rank());
        for(int i = 0; i < chunk_size; i++) 
            printf("outstreamlist size %d p %d s %d\n", i, outlist_comm[i].primary_size(), outlist_comm[i].secondary_size());
        statistics.end("run => work_thread => send => get_mutex");
    }
    statistics.end("run => work_thread => save_outgoing");
}

void AsyncNode::clear_retired_rays() {
//    statistics.start("run => message_thread => clear_retired => get mutex");
    std::unique_lock <std::mutex> lock(retired_mutex); 
//    while(outList_empty() && !ps->is_proc_idle()) {
//        comm->os<<"mthread wait retired rays\n";
//        retired_cond_full.wait(lock);
//    }
//    statistics.end("run => message_thread => clear_retired => get mutex");
    while(!retiredRays.empty()) {
        RetiredRays* rays = retiredRays.back();
        comm->os<<"retired ray size "<<rays->size<<"\n"; 
        retiredRays.pop_back();
        RayStreamList::read_from_device_buffer(outlist_comm, rays->data, rays->size, rays->primary, chunk_size, comm->get_rank());
        delete rays;
    }
}

void AsyncNode::compact_retired_rays() {

    int primary_size = 0;    
    int secondary_size = 0;    
    for(int i = 0; i < retiredRays.size(); i++) {
        if(retiredRays[i]->primary) 
            primary_size += retiredRays[i]->size;
        else
            secondary_size += retiredRays[i]->size;
    }
    RetiredRays* retired_primary = new RetiredRays(primary_size, true);
    RetiredRays* retired_secondary = new RetiredRays(secondary_size, false);
    
    float *primary_ptr = retired_primary->data;
    float *secondary_ptr = retired_secondary->data;
    for(int i = 0; i < retiredRays.size(); i++) {
        if(retiredRays[i]->primary) { 
            int copy_size = retiredRays[i]->size * PRIMARY_WIDTH;
            memcpy(primary_ptr, retiredRays[i]->data, copy_size * sizeof(float));
            primary_ptr += copy_size;
        } else {
            int copy_size = retiredRays[i]->size * SECONDARY_WIDTH;
            memcpy(secondary_ptr, retiredRays[i]->data, copy_size * sizeof(float));
            secondary_ptr += copy_size;
        }
    }

    while(!retiredRays.empty()) {
        RetiredRays* rays = retiredRays.back();
        retiredRays.pop_back();
        delete rays;
    }

    sort_ray_array(retired_primary->data, primary_size, chunk_size, true);
    sort_ray_array(retired_secondary->data, secondary_size, chunk_size, false);

    retiredRays.emplace_back(retired_primary);
    retiredRays.emplace_back(retired_secondary);

    //clear retired rays
    //copy priamry and secondary to retired ryas

}

int AsyncNode::load_incoming_buffer(float **rays, bool primary, int thread_id) {
    //printf("load incoming buffer\n");
    comm->os<<"rthread load incoming buffer inlist size "<<inlist_comm.size()<<"\n";
  
    if(ps->has_new_chunk() || inlist_render.empty()) 
        return -2 - ps->get_current_chunk(); 

    if(ps->Exit()) {
        printf(" recv exit \n");
        return -1;
    }
    
    std::unique_lock <std::mutex> work_thread_lock(ps->thread_mutex); 
    std::unique_lock <std::mutex> lock(inlist_render.mutex); 
    
    struct RaysStream *rays_stream;
    if(primary && inlist_render.primary_size() > 0) {
        rays_stream = inlist_render.get_primary();
    } else if (!primary && inlist_render.secondary_size() > 0) {
        rays_stream = inlist_render.get_secondary();
    } else {
        return 0;
    }
    lock.unlock();
    work_thread_lock.unlock();

    statistics.start("run => work_thread => load_incoming_buffer-copy");

    int copy_size = rays_stream->size;
    comm->os<<"all rays chunk "<<ps->get_current_chunk()<<" size "<<copy_size<<"\n";
    int width = rays_stream->width;
    printf("copy primary size %d\n", copy_size);
    memcpy(*rays, rays_stream->get_data(), ps->get_stream_store_capacity() * width * sizeof(float)); 
    
    delete rays_stream;
    statistics.end("run => work_thread => load_incoming_buffer-copy");
    return copy_size;
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

        comm->recv_message(ps, wk->outlist_comm, wk->inlist_comm, wk->wait_respond);
        wk->wait_respond = false;
        
        statistics.end("run => message_thread => recv_message");
         
        if(ps->Exit()) {
            break;
        }

        wk->loop_check(1);
        statistics.start("run => message_thread => inlist not empty");
        
        if(!wk->inlist_comm.empty()) 
            wk->inlist_comm.cond_full.notify_one();

        statistics.end("run => message_thread => inlist not empty");
        wk->loop_check(3);
        usleep(500);

        if(ps->Exit()) {
            break;
        }
        
        wk->loop_check(3.5);
        //clear_outlist();
        
        statistics.start("run => message_thread => clear_retired");
        wk->clear_retired_rays();
        statistics.end("run => message_thread => clear_retired");
        
        statistics.start("run => message_thread => send_message");
        wk->loop_check(3.7);
        wk->send_message();
        wk->loop_check(4);
        comm->purge_completed_mpi_buffers();
        statistics.end("run => message_thread => send_message");
        statistics.start("run => message_thread => new_chunk");
        statistics.end("run => message_thread => new_chunk");

        statistics.start("run => message_thread => sleep");
        usleep(100);
        statistics.end("run => message_thread => sleep");
        statistics.end("run => message_thread => loop");
	}
    wk->inlist_comm.cond_full.notify_all();
    statistics.end("run => message_thread => loop");
    statistics.end("run => message_thread");
    comm->os <<" end message thread"<<ps->all_thread_waiting()<<"\n";
    comm->os <<" inlist "<< wk->inlist_comm.size()
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
        comm->os<<"rthread start  " << current_chunk << "\n";
        if(iter != 0) {        
            std::unique_lock <std::mutex> inlock(inlist_comm.mutex); 
            
            //compact_retired_rays()       
            
            //retired_cond_full.notify_one();
            
            if(inlist_comm.empty()) {
                ps->set_proc_idle();
                while (inlist_comm.empty() && !ps->Exit()) {
                    comm->os<<"rthread waiting \n";
                    inlist_comm.cond_full.wait(inlock);
                }
                if(ps->Exit()) {
                //    rodent_save_light_field();
                    break;
                }
            }
            RayStreamList::swap(inlist_comm, inlist_render);
            inlist_comm.empty_notify(); //tell mthread inlist size changed 
            inlist_render.empty_notify(); //tell mthread inlist size changed 
            ps->set_proc_busy(comm->get_rank());
        }

        statistics.start("run => render thread => get start ");
        launch_rodent_render(camera, deviceNum, iter==0);

        iter++;
    } while(!ps->Exit());
    comm->os<<" render thread times "<<iter<<"\n";
    mthread.join();
    return;
}

