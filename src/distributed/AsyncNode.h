// Original Master-Worker mode 
struct AsyncNode : public Node {

    AsyncNode(struct Communicator *comm, struct ProcStatus *ps);
    
    ~AsyncNode();

    int get_sent_list(); 
    
    void send_message(); 

    static void message_thread(void* tmp);
    
    void run(Scheduler * camera);
};

AsyncNode::AsyncNode(struct Communicator *comm, struct ProcStatus *ps)
    :Node(comm, ps)
{
    printf("new AsyncNode\n");
}

AsyncNode::~AsyncNode() {
    printf("delete AsyncNode\n");
}

// get chunk retun chunk id
int AsyncNode::get_sent_list() {
    int chunk_size = ps->get_chunk_size();
    int dst_loaded = -1, dst_unloaded = -1;
    int max_loaded = 0,  max_unloaded = 0; 
   
    for(int i = 0; i < chunk_size; i++) {
        if(ps->is_local_chunk(i)) continue;
        int cur_size = rlm->outList_size(i); 
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
    if (ps->all_thread_waiting() && rlm->inList_empty()) {
        if(rlm->allList_empty() ) {
            ps->set_proc_idle(comm->get_rank());
            if(ps->all_proc_idle() && ps->all_rays_received()) {
                QuitMsg quit_msg(comm->get_rank(), get_tag()); 
                comm->send_message(&quit_msg, ps);
                ps->set_exit();
                rlm->inList.cond_full.notify_all();
                return;
            } else {
                StatusMsg status(comm->get_rank(), comm->get_master(), ps->get_status(), ps->get_current_chunk(), comm->get_size(), get_tag()); 
                comm->send_message(&status, ps);
            }
        } else if(rlm->outList_empty(ps)) {
            ///load new chunk
            int chunk_size = ps->get_chunk_size();
            ps->switch_current_chunk();
            int current_chunk = ps->get_current_chunk();
            rlm->inList.cond_full.notify_all();
            return;
        } else {
            comm->os<<"still has rays to send\n";
        }
    } else {
        ps->set_proc_busy(comm->get_rank());
    }

    int cId = get_sent_list();
    if (cId >= 0) {
        comm->os<<"mthread outlist empty "<<cId<<" "<<rlm->outList_size(cId)<<"\n";
        if(cId >= 0) {
            int dst_proc = ps->get_proc(cId);
            comm->os<<"chunk "<<cId<<" get proc "<<dst_proc<<"\n";
            if(dst_proc >= 0) {
                ps->set_proc_busy(dst_proc);
                statistics.start("run => message_thread => send_message => new RayMsg");
                RayMsg *ray_msg = rlm->export_ray_msg(cId, comm->get_rank(), dst_proc, false, get_tag()); 
                statistics.end("run => message_thread => send_message => new RayMsg");
                comm->send_message(ray_msg, ps);
                delete ray_msg;
            } else {
                error("ray msg");
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

        //std::unique_lock <std::mutex> lock(wk->rlm->out_mutex); 
        if(!recv && ps->is_proc_idle() && !ps->Exit())
            recv = comm->recv_message(ps, true);
        else
            recv = comm->recv_message(ps, false);
        //lock.unlock();
        
        statistics.end("run => message_thread => recv_message");
         
        if(ps->Exit()) {
            statistics.end("run => message_thread => loop");
            break;
        }
        wk->loop_check(1);
        statistics.start("run => message_thread => inlist not empty");
        
        if(!wk->rlm->inList_empty()) 
            wk->rlm->inList.cond_full.notify_one();

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
    comm->os <<" inlist "<< wk->rlm->inList_size()
             <<" recv "<<ps->global_rays[comm->get_rank() + comm->get_size()]
             <<std::endl;
    wk->rlm->inList.cond_full.notify_all();
    statistics.end("run => message_thread");
    return;
} 

void AsyncNode::run(Scheduler * camera) {
    
    comm->os <<" start run message thread \n";

    int deviceNum = ps->get_dev_num();
    int iter = 0; 
    std::vector<std::thread> workThread;
   
    std::thread mthread(message_thread, this);
    do {
        //load new chunk;
        int current_chunk = ps->get_current_chunk();
        comm->os<<"mthread copy new chunk " << current_chunk << "\n";
        rlm->copy_to_inlist(current_chunk, comm->get_rank());
        
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

        for(int i = 0; i < deviceNum; i++) 
            rodent_unload_chunk_data(i); 

        iter++;
    } while(!ps->Exit());
    mthread.join();
    return;
}

