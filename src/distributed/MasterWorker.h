struct RetiredRays {
    float* data;
    int size;
    int width;
    bool primary;
    RetiredRays(float* rays, int size, int width)
        :size(size), width(width) 
    {
        primary = width == 21;
        int capacity = size * width;
        data = new float[capacity]; 
        memcpy((char*)data, (char*)rays, capacity * sizeof(float));
    }

    ~RetiredRays(){
        delete[] data;
    }
};

// Original Master-Worker mode 
struct MWNode : public Node {

    MWNode(struct Communicator *comm, struct ProcStatus *ps);
    
    ~MWNode();

    int get_unloaded_chunk();
    
    bool rayList_empty();
    
    bool all_queue_empty();

    bool outList_empty();
    
    void clear_outlist();
    
    void save_outgoing_buffer(float *rays, size_t size, bool primary);

    int get_sent_list(); 
    
    void master_send_message();
    
    void worker_send_message(); 

    static void message_thread(void* tmp);
    
    void run(ImageDecomposition * camera);
    
    RayList *rayList;

        
    std::vector<RetiredRays*> retiredRays; 
};

MWNode::MWNode(struct Communicator *comm, struct ProcStatus *ps)
    :Node(comm, ps)
{
    printf("new MWNode\n");
    int chunk_size = ps->get_chunk_size();
    rayList = new RayList[chunk_size];
    inList.set_capacity(ps->get_stream_logic_capacity(),
                        ps->get_stream_store_capacity());
}

MWNode::~MWNode() {
    printf("delete MWNode\n");
}


int MWNode::get_unloaded_chunk(){
    int unloaded_chunk = -1;
    int max_ray_size = 0;
    for(int i = 0; i < ps->get_chunk_size(); i++) {
        if(rayList[i].size() > max_ray_size && ps->get_proc(i) == -1) {
            unloaded_chunk = i;  
            max_ray_size = rayList[i].size();
        }
    }
    return unloaded_chunk;
}

bool MWNode::outList_empty() {  
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++) {
        if( !rayList[i].empty() ) {
            return false;
        }
    }
    return true;
}

bool MWNode::all_queue_empty() { 
    for(int i = 0; i < ps->get_chunk_size(); i++) {
        if(!rayList[i].empty()) {return false;}
    }
    return true;
}

bool MWNode::rayList_empty() {  
//    printf("if raylist empty\n");
    int chunk_size = ps->get_chunk_size();
    bool res = true;
    for(int i = 0;i < chunk_size; i++)
        if(!rayList[i].empty()) { res = false; break; }
    res &= inList.empty();
    return res;
}


void MWNode::clear_outlist() {  
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++) {
        rayList[i].get_primary().clear();
        rayList[i].get_secondary().clear();
    }
}

void MWNode::save_outgoing_buffer(float *rays, size_t size, bool primary){
    std::lock_guard <std::mutex> lock(out_mutex); 

    int width = primary?21:14; 
    RetiredRays *retired_rays = new RetiredRays(rays, size, width);
    retiredRays.emplace_back(retired_rays);
}

// get chunk retun chunk id
int MWNode::get_sent_list() {
    int chunk_size = ps->get_chunk_size();
    int dst_loaded = -1, dst_unloaded = -1;
    int max_loaded = 0,  max_unloaded = 0; 
   
 //   if(comm->get_rank() == 1) { 
 //       printf("proc %d raylist status: ", comm->get_rank());
 //       for(int i = 0; i < chunk_size; i++)
 //           printf("|chunk %d  proc %d rays %d ", i, ps->get_proc(i), rayList[i].get_primary().get_size());
 //       printf("| \n");
 //   }
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
            if(ps->all_thread_waiting() || max_loaded > 0.4 * ps->get_stream_logic_capacity()) 
                 return dst_loaded;
            return -1;
        } 
        else 
           return dst_unloaded; 
    }
}

void MWNode::master_send_message() {
    //master idle
    //printf("master all status all thread wait %d inlist %d %d %d %d\n", 
    //        ps->all_thread_waiting(), inList.size(), ps->all_proc_idle(), ps->all_rays_received(), rayList_empty());
//    for(int i = 0; i < ps->get_chunk_size(); i++)
//        printf("M %d ",rayList[i].size());
//    printf("\n");
//    comm->os<<"mthread status inlist size "<<inList.primary_size()<<" "<<inList.secondary_size()<<" thread wait "<<ps->all_thread_waiting()<<"\n";
    
    if(ps->all_thread_waiting() && inList.empty() && rayList_empty()) {
        comm->os<< "mthread if all thread idle all queue empty, set itself idle\n";
        ps->set_proc_idle(comm->get_rank());
    }
    if(ps->all_proc_idle() && ps->all_rays_received() && rayList_empty()) {
        QuitMsg quit_msg(comm->get_rank(), get_tag()); 
        comm->send_message(&quit_msg, ps);
        ps->set_exit();
        comm->os<<"set exit\n";
        inList_not_empty.notify_all();
        return;
    }

    //proc idle
    int idle_proc = ps->get_idle_proc();
    if(idle_proc >= 0 && ps->all_rays_received()) {
        //find a unloaded chunk contents unprocessed rays
        // if idle proc loaded chunk not empty ???

        if(idle_proc == comm->get_rank()) {
            comm->os<<"master idle\n"; 
        } else {
            int unloaded_chunk = get_unloaded_chunk(); 
            //find a idle proc send chunk id to it
            if(unloaded_chunk > -1) {
                comm->os<<"mthread scheduleMsg"<<"\n";
                ps->update_chunk(idle_proc, unloaded_chunk);
                RayMsg ray_msg(rayList[unloaded_chunk], comm->get_rank(), idle_proc, unloaded_chunk, false, get_tag()); 
                
                comm->os<<"mthread schedule RayMsg "<<ray_msg.get_chunk()<<"\n";
                comm->send_message(&ray_msg, ps);
                ps->set_proc_busy(idle_proc);
                
                ScheduleMsg msg(comm->get_rank(), ps->get_chunk_proc(), idle_proc, ps->get_chunk_size(), get_tag()); 
                comm->send_message(&msg, ps);
                
          //      return;//////??
            }
        }
    }
    if(!outList_empty()) {

        int cId = get_sent_list();
        if(cId >= 0 && ps->get_proc(cId) >= 0 ) {
            int dst_proc = ps->get_proc(cId);
            comm->os<<"mthread master new RayMsg "<<cId<<" size "<<rayList[cId].size()<<" rank "<<comm->get_rank()<<" dst "<<dst_proc<<"\n";
            
            ps->set_proc_busy(dst_proc);
            
            //printf("MasterWorker.h raylist chunk id %d dst %d  primary size %d width %d : ", 
            //          cId, ps->get_proc(cId), rayList[cId].get_primary().get_size(), rayList[cId].get_primary().get_store_width());
            RayMsg ray_msg(rayList[cId], comm->get_rank(), dst_proc, cId, false, get_tag()); 
            
            comm->os<<"mthread RayMsg "<<ray_msg.get_chunk()<<"\n";
            comm->send_message(&ray_msg, ps);
        }
    }
}

void MWNode::worker_send_message() {
 // printf("master after proc status\n");
//    comm->os<<"mthread status inlist size "<<inList.primary_size()<<" "<<inList.secondary_size()<<" thread wait "<<ps->all_thread_waiting()<<"\n";
//    for(int i = 0; i < ps->get_chunk_size(); i++)
//        printf("w %d ",rayList[i].size());
//    printf("\n");
    if (ps->all_thread_waiting() && inList.empty()) {
        if(rayList_empty()) {
            comm->os<< "mthread send status msg\n";
            int * tmp = ps->get_status();
            comm->os <<tmp[0]<<" "<<tmp[1] <<" "<<tmp[2]<<" "<<tmp[3]<<"\n";
            ps->set_proc_idle(comm->get_rank());
            comm->os<<"mthread local chunk"<<ps->get_local_chunk()<<"\n";
            StatusMsg status(comm->get_rank(), comm->get_master(), ps->get_status(), ps->get_local_chunk(), comm->get_size(), get_tag()); 
            comm->send_message(&status, ps);
            comm->os<< "mthread send status msg success\n";
        } else {
            int chunk_size = ps->get_chunk_size();
            for(int i = 0;i < chunk_size; i++)
                comm->os<<rayList[i].primary_size()<<" "<<rayList[i].secondary_size()<<" \n";
            
            comm->os<<"need another chunk or need send rays\n";
        }
    } else {
      //  comm->os<<"mthread set proc busy\n";
        ps->set_proc_busy(comm->get_rank());
    }
    
    if (ps->has_new_chunk()) {
        inList_not_empty.notify_all();
        for(int i = 0; i < ps->get_chunk_size(); i++) {
            comm->os<<"new chunk raylist "<<rayList[i].size()<<" |";
        }
        comm->os<<"\n";
    }
    if (!outList_empty()) {
        //comm->os<<"mthread outlist empty"<<outList_empty()<<"\n";
        int cId = get_sent_list();
        if(cId >= 0) {
            if(ps->get_proc(cId) >= 0) {
                int dst_proc = ps->get_proc(cId);
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

void MWNode::message_thread(void* tmp) {
  
    MWNode *wk = (struct MWNode*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus * ps = wk->ps;

    statistics.start("run => message_thread");
    printf("%d message thread\n", comm->get_rank());
//    comm->os<<"message thread\n";
    int recv_count = 0;
    while (!ps->Exit()) {
       
        statistics.start("run => message_thread => loop");
        bool recv = false;
       // if(comm->isMaster()){ 
       //     printf("master before recv thread wait %d inlist %d %d %d %d\n", 
       //        ps->all_thread_waiting(), wk->inList.size(), ps->all_proc_idle(), ps->all_rays_received(), wk->rayList_empty());
       // }
        
        //if send list empty also stay here 
        statistics.start("run => message_thread => recv_message");
        do {
      //     if(comm->isMaster()) 
      //         comm->os<<"master wait recv\n";
            recv = comm->recv_message(wk->rayList, &(wk->inList),  ps);
        } while(!recv && ps->is_proc_idle() && !ps->Exit()) ;
        statistics.end("run => message_thread => recv_message");
         
        if(ps->Exit()) {
            statistics.end("run => message_thread => loop");
            break;
        }

        if(!wk->retiredRays.empty()) {
            std::lock_guard <std::mutex> lock(wk->out_mutex); 
            while(!wk->retiredRays.empty()) {
                RetiredRays* rays = wk->retiredRays.back();
                wk->retiredRays.pop_back();
                RayList::read_from_device_buffer(wk->rayList, rays->data, rays->size, rays->primary, comm->get_rank(), ps->get_chunk_size());
            }
        }
        
        statistics.start("run => message_thread => inlist not empty");
        //if new chunk havent loaded, render thread havent start, mthread wait
        
        {
            std::unique_lock <std::mutex> lock(wk->thread_mutex); 
            if(!wk->inList.empty()) {
                if(ps->has_new_chunk() ) {
    //                comm->os<<"wait new chunk\n";
                    while( ps->has_new_chunk()) {
    //                    comm->os<<"wait rthread launch\n";
                        wk->render_start.wait(lock); 
                    }
                } 
                wk->inList_not_empty.notify_one();
            }
        }
        statistics.end("run => message_thread => inlist not empty");

        if(ps->Exit()) {
            statistics.end("run => message_thread => loop");
            break;
        }
        
        statistics.start("run => message_thread => send_message");
        //clear_outlist();
        if(comm->isMaster())
            wk->master_send_message(); 
        else
            wk->worker_send_message();
        
        comm->purge_completed_mpi_buffers();
        statistics.end("run => message_thread => send_message");
        
        statistics.start("run => message_thread => sleep");
        usleep(500);
        statistics.end("run => message_thread => sleep");
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

void MWNode::run(ImageDecomposition * camera) {
    
    comm->os <<" start run message thread \n";


    int deviceNum = ps->get_dev_num();
    int iter = 0; 
    std::vector<std::thread> workThread;
    
    std::thread mthread(message_thread, this);
    //    all proc start proc 0 send schedule? 
    while(1) {
        comm->os<<"new chunk raylist local chunk "<<ps->get_local_chunk()<<"\n";
        for(int i = 0; i < ps->get_chunk_size(); i++)
            comm->os<< rayList[i].size()<<"  ";
        comm->os<<"\n";

        comm->os<<"                                   new chunk raylist left chunk " << ps->get_local_chunk()<<" ";
        for(int i = 0; i < ps->get_chunk_size(); i++)
            comm->os<<" "<<rayList[i].size();
        comm->os<<"\n";

        if(ps->has_new_chunk()) {
            ps->chunk_loaded();

            comm->os<<comm->get_rank()<<" before set render start"<<"\n";
            for(int i = 0; i < deviceNum; i++) 
                workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, iter == 0));
             
            //set chunk loaded
            comm->os<<" set render start "<<"\n";
            comm->os<<" set render start \n";
            
            render_start.notify_all(); 
            comm->os<<" render start notify \n";
            printf("%d render start notify \n", comm->get_rank());

            for(auto &thread: workThread)
                thread.join();

            workThread.clear();
        } else {
            ps->set_exit();
            break;
        }
        iter++;
    }
    printf("%d worker exit\n", comm->get_rank());
    mthread.join();
    return;
    printf("run MWNode\n");
}

