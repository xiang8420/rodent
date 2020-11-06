// Original Master-Worker mode 
struct AsyncNode : public Node {

    AsyncNode(struct Communicator *comm, struct ProcStatus *ps);
    
    ~AsyncNode();

    bool rayList_empty();
    
    bool all_queue_empty();

    bool outList_empty();
    
    void clear_outlist();
    
    void save_outgoing_buffer(float *rays, size_t size, bool primary);

    int get_sent_list(); 
    
    void send_message(); 

    static void message_thread(void* tmp);
    
    void run(ImageDecomposition * camera);
    
    RayList *rayList;
        
    std::vector<RetiredRays*> retiredRays; 
};

AsyncNode::AsyncNode(struct Communicator *comm, struct ProcStatus *ps)
    :Node(comm, ps)
{
    printf("new AsyncNode\n");
    int chunk_size = ps->get_chunk_size();
    rayList = new RayList[chunk_size];
    inList.set_capacity(ps->get_stream_logic_capacity(),
                        ps->get_stream_store_capacity());
}

AsyncNode::~AsyncNode() {
    printf("delete AsyncNode\n");
}

bool AsyncNode::outList_empty() {  
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++) {
        if( !rayList[i].empty() && !ps->is_local_chunk(i)) {
            return false;
        }
    }
    return true;
}

bool AsyncNode::all_queue_empty() { 
    for(int i = 0; i < ps->get_chunk_size(); i++) {
        if(!rayList[i].empty()) {return false;}
    }
    return true;
}

bool AsyncNode::rayList_empty() {  
//    printf("if raylist empty\n");
    int chunk_size = ps->get_chunk_size();
    bool res = true;
    for(int i = 0;i < chunk_size; i++)
        if(!rayList[i].empty()) { res = false; break; }
    res &= inList.empty();
    return res;
}


void AsyncNode::clear_outlist() {  
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++) {
        rayList[i].get_primary().clear();
        rayList[i].get_secondary().clear();
    }
}

void AsyncNode::save_outgoing_buffer(float *rays, size_t size, bool primary){
    std::lock_guard <std::mutex> lock(out_mutex); 
    
    int width = primary ? 21 : 14; 
    
    int* ids = (int*)(rays);
    for(int i = 0; i < 5; i ++) {
        comm->os<<"| "<< ids[i * width] <<" "<< ids[i * width + 9] << " ";
    }
    comm->os<<"\n";
    
    RetiredRays *retired_rays = new RetiredRays(rays, size, width);
    retiredRays.emplace_back(retired_rays);
}

// get chunk retun chunk id
int AsyncNode::get_sent_list() {
    int chunk_size = ps->get_chunk_size();
    int dst_loaded = -1, dst_unloaded = -1;
    int max_loaded = 0,  max_unloaded = 0; 
   
    for(int i = 0; i < chunk_size; i++) {
        if(ps->is_local_chunk(i)) continue;

        if(ps->get_proc(i) >= 0 && rayList[i].size() > max_loaded ) {
            max_loaded = rayList[i].size();
            dst_loaded = i;
        }
        if(ps->get_proc(i) < 0 && rayList[i].size() > max_unloaded) {
            max_unloaded = rayList[i].size();
            dst_unloaded = i;
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
    printf("async node send message\n");
    for(int i = 0; i < ps->get_chunk_size(); i++)
        printf("w %d ",rayList[i].size());
    printf("\n");
    //comm->os<<"all thread waitting "<<ps->all_thread_waiting() <<"inList size"<<inList.size()<<"\n";
    if (ps->all_thread_waiting() && inList.empty()) {
        if(rayList_empty()) {
            ps->set_proc_idle(comm->get_rank());
            if(ps->all_proc_idle() && ps->all_rays_received()) {
                QuitMsg quit_msg(comm->get_rank(), get_tag()); 
                comm->send_message(&quit_msg, ps);
                ps->set_exit();
                inList_not_empty.notify_all();
                return;
            } else {
                StatusMsg status(comm->get_rank(), comm->get_master(), ps->get_status(), ps->get_current_chunk(), comm->get_size(), get_tag()); 
                comm->send_message(&status, ps);
            }
        } else if(outList_empty()) {
            ///load new chunk
            int chunk_size = ps->get_chunk_size();
            for(int i = 0;i < chunk_size; i++)
                comm->os<<rayList[i].primary_size()<<" "<<rayList[i].secondary_size()<<" \n";
            ps->switch_current_chunk();
            int current_chunk = ps->get_current_chunk();
            comm->os<<"need another chunk or need send rays "<< current_chunk << "\n"; ;
            inList_not_empty.notify_all();
            return;
        } else {
            comm->os<<"still has rays to send\n";
        }
    } else {
        ps->set_proc_busy(comm->get_rank());
    }
    
    if (!outList_empty()) {
        comm->os<<"mthread outlist empty"<<outList_empty()<<"\n";
        int cId = get_sent_list();
        if(cId >= 0) {
            int dst_proc = ps->get_proc(cId);
            comm->os<<"chunk "<<cId<<" get proc "<<dst_proc<<"\n";
            if(dst_proc >= 0) {
                ps->set_proc_busy(dst_proc);
                RayMsg ray_msg(rayList[cId], comm->get_rank(), dst_proc, cId, false, get_tag()); 
                comm->send_message(&ray_msg, ps);
            } else {
                RayMsg ray_msg(rayList[cId], comm->get_rank(), comm->get_master(), cId, false, get_tag()); 
                comm->send_message(&ray_msg, ps);
            }
        }
    }

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
        bool recv = false;

        statistics.start("run => message_thread => recv_message");

        if(!recv && ps->is_proc_idle() && !ps->Exit())
            recv = comm->recv_message(wk->rayList, &(wk->inList),  ps, true);
        else
            recv = comm->recv_message(wk->rayList, &(wk->inList),  ps, false);
        
        statistics.end("run => message_thread => recv_message");
         
        if(ps->Exit()) {
            statistics.end("run => message_thread => loop");
            break;
        }
        wk->loop_check(1);

        if(!wk->retiredRays.empty()) {
            std::lock_guard <std::mutex> lock(wk->out_mutex); 
            while(!wk->retiredRays.empty()) {
                RetiredRays* rays = wk->retiredRays.back();
                wk->retiredRays.pop_back();
                RayList::read_from_device_buffer(wk->rayList, rays->data, rays->size, rays->primary, comm->get_rank(), ps->get_chunk_size());
            }
        }
        statistics.start("run => message_thread => inlist not empty");
        
        if(!wk->inList.empty()) 
            wk->inList_not_empty.notify_one();

        statistics.end("run => message_thread => inlist not empty");
        wk->loop_check(3);

        if(ps->Exit()) {
            statistics.end("run => message_thread => loop");
            break;
        }
        
        wk->loop_check(3.5);
        statistics.start("run => message_thread => send_message");
        //clear_outlist();
        wk->send_message();
        
        wk->loop_check(4);
        comm->purge_completed_mpi_buffers();
        statistics.end("run => message_thread => send_message");
        {
            if(ps->has_new_chunk()) {
                std::unique_lock <std::mutex> lock(wk->thread_mutex); 
                while( ps->has_new_chunk()) {
                    wk->render_start.wait(lock); 
                }
            }
        }
        wk->loop_check(5);

        usleep(500);
        statistics.end("run => message_thread => loop");
	}
    comm->os <<" end message thread"<<ps->all_thread_waiting()<<"\n";
    comm->os <<" inlist "<< wk->inList.size()
             <<" recv "<<ps->global_rays[comm->get_rank() + comm->get_size()]
             <<std::endl;
    wk->inList_not_empty.notify_all();
    statistics.end("run => message_thread");
    return;
} 

void AsyncNode::run(ImageDecomposition * camera) {
    
    comm->os <<" start run message thread \n";

    int deviceNum = ps->get_dev_num();
    int iter = 0; 
    std::vector<std::thread> workThread;
   
    std::thread mthread(message_thread, this);
    do { 
        //load new chunk;
        int current_chunk = ps->get_current_chunk();
        comm->os<<"mthread copy new chunk " << current_chunk << "\n";
        if(rayList[current_chunk].size() > 0) {
            RaysArray &primary   = rayList[current_chunk].get_primary(); 
            RaysArray &secondary = rayList[current_chunk].get_secondary(); 
            
            int* ids = (int*)(primary.data);
            int copy_size = primary.get_size(); 
            for(int i = 0; i < std::min(5, copy_size); i ++) {
                comm->os<<" || " <<ids[i * 21]<<" "<<ids[i * 21 + 9];
            }
            comm->os<<"\n";
            ids = (int*)(secondary.data);
            copy_size = secondary.get_size(); 
            comm->os<<"mthread copy new chunk rays" <<copy_size<<"\n";
            for(int i = 0; i < std::min(5, copy_size); i ++) {
                comm->os<<" || " <<ids[i * 14]<<" "<<ids[i * 14 + 9];
            }
            comm->os<<"\n";

            inList.read_from_ptr(primary.get_data(), 0, primary.get_size(), true, comm->get_rank());
            inList.read_from_ptr(secondary.get_data(), 0, secondary.get_size(), false, comm->get_rank());
            rayList[current_chunk].clear();
            //read to inList
        //    clear_outlist();
        }
        ps->chunk_loaded();

        comm->os<<comm->get_rank()<<" before set render start"<<"\n";
        for(int i = 0; i < deviceNum; i++) 
            workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, iter == 0));
        
        render_start.notify_all(); 
        comm->os<<" render start notify \n";
        printf("%d render start notify \n", comm->get_rank());

        for(auto &thread: workThread)
            thread.join();

        workThread.clear();

        rodent_unload_bvh(0); 

        iter++;
    } while(!ps->Exit());
    mthread.join();
    return;
}

