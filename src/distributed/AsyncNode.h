struct RetiredRays {
    float* data;
    int size;
    bool primary;
    int width;

    RetiredRays(float* rays, int size, bool primary)
        :size(size), primary(primary) 
    {
        int width = primary ? PRIMARY_WIDTH : SECONDARY_WIDTH;
        int capacity = size * width;
        data = new float[capacity]; 
        printf("new retired ray size %d width %d \n", size, width);
        memcpy((char*)data, (char*)rays, capacity * sizeof(float));
    }
    RetiredRays(int size, bool primary)
        :size(size), primary(primary) 
    {
        int width = primary ? PRIMARY_WIDTH : SECONDARY_WIDTH;
        int capacity = size * width;
        data = new float[capacity]; 
    }

    ~RetiredRays() {
        delete[] data;
    }
};


// Original Master-Worker mode 
struct AsyncNode : public Node {

    AsyncNode(Communicator *, ProcStatus *, Scheduler *);
    
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
    
    void postprocess(int); 
    
    void run(Camera *);

    bool outList_empty(); 

    bool allList_empty();

    RayStreamList  inlist_comm;  
    RayStreamList * outlist_comm;
    
    RayStreamList inlist_render;  
    //RayStreamList * outlist_render;

    std::vector<RetiredRays*> retiredRays; 
    std::mutex retired_mutex;
    std::condition_variable retired_cond_full; 

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

void swap_array(float *rays, int k, int j, int width) {
    int *tmp = new int[width];
    memcpy(tmp, &rays[width * k], width * sizeof(float));
    memcpy(&rays[width * k], &rays[width * j], width * sizeof(float));
    memcpy(&rays[width * j], tmp, width * sizeof(float));
    delete[] tmp; 
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
    std::unique_lock <std::mutex> lock(outlist_comm[0].mutex); 
    RayMsg *msg = new RayMsg(outlist_comm[cId], rank, dst, cId, idle, tag); 
    return msg;
}

bool AsyncNode::allList_empty() {
   // std::lock_guard <std::mutex> lock(out_mutex); 
    for(int i = 0; i < chunk_size; i++)
        if(!outlist_comm[i].empty()) 
            return false;
    comm->os<<" all list empty "<<inlist_comm.empty()<<" "<<retiredRays.empty()<<"\n";
    return true;//inlist_comm.empty();// && retiredRays.empty();
}

AsyncNode::AsyncNode(Communicator *comm, ProcStatus *ps, Scheduler* scheduler)
    :Node(comm, ps, scheduler)
{
    printf("new AsyncNode\n");
    int store_capacity = ps->get_stream_store_capacity();
    int logic_capacity = ps->get_stream_logic_capacity();
   
    chunk_size = SIMPLE_TRACE ? ps->get_chunk_size() + 1 : ps->get_chunk_size(); 
    outlist_comm = new RayStreamList[chunk_size];
    for(int i = 0; i < chunk_size; i++)
        outlist_comm[i].set_capacity(logic_capacity, store_capacity);
    inlist_comm.set_capacity(logic_capacity, store_capacity);
    inlist_render.set_capacity(logic_capacity, store_capacity);
    
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
//    comm->os<<"proc idle "<<ps->is_proc_idle() <<"inlist_comm size"<<inlist_comm.size()<<"\n";
    if (ps->is_proc_idle()) {
        if(allList_empty()) {
            if(ps->all_proc_idle() && ps->all_rays_received()) {
                comm->os<<"mthread send quit\n";
                QuitMsg quit_msg(comm->get_rank(), comm->get_tag()); 
                comm->send_message(&quit_msg, ps);
                ps->set_exit();
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
        } else if (outList_empty()) {
            statistics.start("run => message_thread => send_message => switch_chunk");
            ///load new chunk
            ps->switch_current_chunk(outlist_comm);
            
            int current_chunk = ps->get_current_chunk();
            copy_to_inlist(current_chunk);
            inlist_comm.cond_full.notify_all();
            inlist_render.cond_full.notify_all();
        } else {
            //outlist need to send
        }
    }
    int cId = get_sent_list();
    do {
        if(cId >= 0) {
            int dst_proc = ps->get_proc(cId);
            if(dst_proc >= 0) {
                ps->set_proc_busy(dst_proc);
                statistics.start("run => message_thread => send_message => new RayMsg");
                RayMsg *ray_msg = export_ray_msg(cId, comm->get_rank(), dst_proc, false, comm->get_tag()); 
                statistics.end("run => message_thread => send_message => new RayMsg");
                statistics.start("run => message_thread => send_message => comm->send");
                comm->send_message(ray_msg, ps);
                statistics.end("run => message_thread => send_message => comm->send");
                delete ray_msg;
            } else {
                error("ray msg");
            }
        }
        cId = get_sent_list();
    } while(cId >= 0 && ps->is_proc_idle());
}

void AsyncNode::save_outgoing_buffer(float *rays, size_t size, bool primary) {

    RetiredRays *retired_rays = new RetiredRays(rays, size, primary);
    std::lock_guard <std::mutex> thread_lock(retired_mutex); 
    
    scheduler->pass_record->write_send(rays, size, primary, scheduler->chunk_manager->local_chunks.current); 

    retiredRays.emplace_back(retired_rays);

    if(retiredRays.size() > 10)
        clear_retired_rays();
}

void AsyncNode::clear_retired_rays() {
    statistics.start("run => message_thread => clear_retired => get mutex");
    compact_retired_rays();       
    std::unique_lock <std::mutex> lock(outlist_comm[0].mutex); 
    while(!retiredRays.empty()) {
        RetiredRays* rays = retiredRays.back();
        comm->os<<"retired ray size "<<rays->size<<"\n"; 
        retiredRays.pop_back();
        statistics.start("run => message_thread => clear_retired => read_from_device");
        RayStreamList::read_from_device_buffer(outlist_comm, rays->data, rays->size, rays->primary, chunk_size, comm->get_rank());
        statistics.end("run => message_thread => clear_retired => read_from_device");
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
//    comm->os<<"rthread load incoming buffer inlist size "<<inlist_comm.size()<<"\n";
  
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

    statistics.start("run => wthread => load_incoming_buffer-copy");

    int copy_size = rays_stream->size;
    comm->os<<"all rays chunk "<<ps->get_current_chunk()<<" size "<<copy_size<<"\n";
    int width = rays_stream->width;
    printf("copy primary size %d\n", copy_size);
    memcpy(*rays, rays_stream->get_data(), ps->get_stream_store_capacity() * width * sizeof(float)); 
  
    scheduler->pass_record->write_recv(copy_size, scheduler->chunk_manager->local_chunks.current);

    delete rays_stream;
    statistics.end("run => wthread => load_incoming_buffer-copy");
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
         
        if(ps->Exit()) break;

        wk->loop_check(1);
        statistics.start("run => message_thread => inlist not empty");
        
        if(!wk->inlist_comm.empty()) 
            wk->inlist_comm.cond_full.notify_one();

        statistics.end("run => message_thread => inlist not empty");
        wk->loop_check(3);
        usleep(500);

        if(ps->Exit()) break;
        
        wk->loop_check(3.5);
        //clear_outlist();
        
        statistics.start("run => message_thread => send_message");
        wk->loop_check(3.7);
        wk->send_message();
        wk->loop_check(4);
        comm->purge_completed_mpi_buffers();
        statistics.end("run => message_thread => send_message");

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

void AsyncNode::run(Camera *cam) {
    
    ps->reset();
    scheduler->preprocess(cam, comm->get_size(), SIMPLE_TRACE, false);
    
    comm->os <<" start run message thread \n";

    int deviceNum = ps->get_dev_num();
    int iter = 0;
  
    std::thread mthread(message_thread, this);
    do {
        //load new chunk;
        int current_chunk = ps->get_current_chunk();
        if(iter == 0) { 
            int* region = scheduler->get_render_block();
            int new_rays = scheduler->get_spp() * (region[2] - region[0]) * (region[3] - region[1]);
            scheduler->pass_record->write_recv(new_rays, scheduler->chunk_manager->local_chunks.current); 
        }

        comm->os<<"rthread start  " << current_chunk << "\n";
        if(iter != 0) {
            std::unique_lock <std::mutex> inlock(inlist_comm.mutex); 
            
            if(!retiredRays.empty()) {
                clear_retired_rays();
            }
            
            if(inlist_comm.empty()) {
                ps->set_proc_idle();
                while (inlist_comm.empty() && !ps->Exit()) {
                    comm->os<<"rthread waiting \n";
                    inlist_comm.cond_full.wait(inlock);
                }
                if(ps->Exit()) {
                    break;
                }
            }
            RayStreamList::swap(inlist_comm, inlist_render);
            inlist_comm.empty_notify(); //tell mthread inlist size changed 
            inlist_render.empty_notify(); //tell mthread inlist size changed 
            ps->set_proc_busy(comm->get_rank());
        }
        statistics.start("run => wthread => render");
        launch_rodent_render(deviceNum, iter==0);
        float t = statistics.end("run => wthread => render");
        scheduler->pass_record->write_time(t);
        iter++;
    } while(!ps->Exit());
    comm->os<<" render thread times "<<iter<<"\n";
    mthread.join();
    if(SIMPLE_TRACE /* && cam->iter > 1*/)
        postprocess(cam->iter);
    return;
}

void AsyncNode::postprocess(int iter) {
    int * chunk_hit = rodent_get_chunk_hit();
    int size = CHUNK_HIT_RES * CHUNK_HIT_RES * 6 * ps->get_chunk_size();
    int *reduce_buffer = new int[size];
    printf("start reduce %d\n", comm->get_rank());
    {
        statistics.start("run => light field => bcast ");
    //    comm->update_chunk_hit(chunk_hit, reduce_buffer, size);
        MPI_Reduce(chunk_hit, reduce_buffer, size, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
        MPI_Bcast(reduce_buffer, size, MPI_INT, 0, MPI_COMM_WORLD);
        statistics.end("run => light field => bcast ");
    }
    
    //updata render light field
    rodent_update_render_chunk_hit(reduce_buffer, size);
    //reduce ctib 存在了reduce 里
    statistics.start("run => light field => save img  ");
    if(comm->get_rank() == 0 && VIS_CHUNK_HIT) { 
        int* render_chunk_hit = rodent_get_render_chunk_hit_data();
        save_image_ctrb(render_chunk_hit, ps->get_chunk_size(), iter);
    }
    statistics.end("run => light field => save img  ");
    
    //pass_record get chunk speed
    scheduler->set_load_chunk_hit(); 
    delete[] reduce_buffer;
}
