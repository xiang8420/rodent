
// Original Master-Worker mode 
struct MWNode : public Node{
    double recv_t, wait_t, send_t, total_t, read_t, write_t, st, ed;

    MWNode(struct Communicator *comm, struct ProcStatus *ps);
    
    ~MWNode();
    
    int get_dst_worker_id();

    int get_max_rays_chunk(bool unloaded); 

    bool all_queue_empty(){ 
        for(int i = 0; i < ps->get_chunk_size(); i++) {
            if(!rayList[i]->empty()) {return false;}
        }
        return true;
    }

    int get_unloaded_chunk();
    
    void check_proc_status();
   
    void send_message();

    static void message_thread(void* tmp);
    
    void run(ImageDecomposition * camera);
    
    void save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary);

    void clear_outlist() {  
        int chunk_size = ps->get_chunk_size();
        for(int i = 0; i < chunk_size; i++) {
            if(rayList[i]->type == "out" ) {
                rayList[i]->primary->clear();
                rayList[i]->secondary->clear();
            }
        }
    }

};

MWNode::MWNode(struct Communicator *comm, struct ProcStatus *ps)
    :Node(comm, ps)
{
    printf("new MWNode\n");
    int chunk_size = ps->get_chunk_size();
    comm->os<< "master chunk size"<<chunk_size<<"\n";
  
    size_t memory[4];
    memory[0] = physical_memory_used_by_process();

    comm->os<<"master set up\n";
    rayList       = new RayList *[chunk_size];
    for(int i = 0; i < chunk_size; i++) {
        rayList[i] = new RayList(10, "out");
    }
    memory[1] = physical_memory_used_by_process();

    inList = new RayStreamList(ps->get_buffer_size());

    memory[2] = physical_memory_used_by_process();
    printf("Memory %ld kb out list %ld in list %ld\n", memory[0], memory[1] - memory[0], memory[2] - memory[1]);
    
    st = clock();
    recv_t = 0; wait_t = 0;
    send_t = 0; total_t = 0;
    write_t = 0; read_t = 0;
}

MWNode::~MWNode() {
    printf("delete MWNode\n");
    for(int i = 0; i < ps->get_chunk_size(); i++){
        delete rayList[i];
    }
    delete inList;
}

void MWNode::save_outgoing_buffer(float *retired_rays, size_t size, size_t capacity, bool primary){
    out_mutex.lock(); 
    int width = primary?21:14; 
    comm->os<<"rthread save outgoing buffer "<< size <<" width "<<width<<"\n"; 
    int* ids = (int*)(retired_rays);
    for(int i = 0; i < 5; i ++) {
        comm->os<<"| "<< ids[i * width] <<" "<< ids[i * width + 9] << " ";
    }
    comm->os<<"\n";
    RayList::read_from_device_buffer(rayList, retired_rays, size, capacity, primary, comm->rank);
    out_mutex.unlock(); 
}

int MWNode::get_unloaded_chunk(){
    int unloaded_chunk = -1;
    for(int i = 0; i < ps->get_chunk_size(); i++) {
        if(rayList[i]->size() > 0 && ps->get_proc(i) == -1) {
            unloaded_chunk = i; break; 
        }
    }
    return unloaded_chunk;
}

void MWNode::send_message() {
    
    //clear_outlist();
    if(comm->isMaster()) {
        //master idle
     //   printf("master all status all thread wait %d inlist %d %d %d %d\n", 
     //           ps->all_thread_waiting(), inList->size(), ps->all_proc_idle(), ps->all_rays_received(), rayList_empty());
     //   for(int i = 0; i < ps->get_chunk_size(); i++)
     //       printf(" %d ",rayList[i]->size());
     //   printf("\n");
        
        if(ps->all_thread_waiting() && inList->empty() && rayList_empty()) {
            comm->os<< "mthread if all thread idle all queue empty, set itself idle\n";
            ps->set_proc_idle(comm->rank);
        }

        if(ps->all_proc_idle() && ps->all_rays_received() && rayList_empty()) {
            QuitMsg quit_msg(comm->rank); 
            comm->send_message(&quit_msg, ps);
            ps->set_exit();
            comm->os<<"set exit\n";
            inList_not_empty.notify_all();
            return;
        }

        //proc idle
        int idle_proc = ps->get_idle_proc();
        out_mutex.lock();
        if(idle_proc >= 0 && ps->all_rays_received()) {
            //find a unloaded chunk contents unprocessed rays
            // if idle proc loaded chunk not empty ???

            if(idle_proc == comm->rank) {
                comm->os<<"master idle\n"; 
            } else {
                int unloaded_chunk = get_unloaded_chunk(); 
                //find a idle proc send chunk id to it
                if(unloaded_chunk > -1) {
                    comm->os<<"mthread scheduleMsg"<<"\n";
                    ps->update_chunk(idle_proc, unloaded_chunk);
                    RayMsg *ray_msg = new RayMsg(rayList[unloaded_chunk], comm->rank, idle_proc, unloaded_chunk, false); 
                    
                    comm->os<<"mthread schedule RayMsg "<<ray_msg->get_chunk()<<"\n";
                    comm->send_message(ray_msg, ps);
                    ps->set_proc_busy(idle_proc);
                    
                    ScheduleMsg *msg = new ScheduleMsg(comm->rank, ps->get_chunk_map(), idle_proc, ps->get_chunk_size()); 
                    comm->send_message(msg, ps);
                    
              //      return;//////??
                }
            }
        }
        if(!outList_empty()) {

            int cId = get_sent_list();
            if(cId >= 0 ) {
                int dst_proc = ps->get_proc(cId);
                comm->os<<"mthread master new RayMsg "<<cId<<" size "<<rayList[cId]->size()<<" rank "<<comm->rank<<" dst "<<dst_proc<<"\n";
                
                ps->set_proc_busy(dst_proc);
                
                int* ids = (int*)(rayList[cId]->get_primary()->get_data());
                for(int i = 0; i < 5; i ++) {
                    printf("%d %d |", ids[i * 16], ids[i * 16 + 9]);
                    comm->os<<"| "<< ids[i * 16] <<" "<< ids[i * 16 + 9] << " ";
                }
                printf("\n");
                
                RayMsg *ray_msg = new RayMsg(rayList[cId], comm->rank, dst_proc, cId, false); 
                
                comm->os<<"mthread RayMsg "<<ray_msg->get_chunk()<<"\n";
                comm->send_message(ray_msg, ps);
            }
        }
        out_mutex.unlock(); 
       // printf("master after proc status\n");
    } else {
    //    comm->os<<"mthread status inlist size "<<inList->size()<<" thread wait "<<ps->all_thread_waiting()<<"\n";
        if (ps->all_thread_waiting() && inList->empty()) {
            if(rayList_empty()) {
                comm->os<< "mthread send status msg\n";
                int * tmp = ps->get_status();
                comm->os <<tmp[0]<<" "<<tmp[1] <<" "<<tmp[2]<<" "<<tmp[3]<<"\n";
                
                ps->set_proc_idle(comm->rank);
                //send to master I am idle
                comm->os<<"mthread local chunk"<<ps->get_local_chunk()<<"\n";
                StatusMsg status(comm->rank, comm->master, ps->get_status(), ps->get_local_chunk(), comm->size); 
                comm->send_message(&status, ps);
                comm->os<< "mthread send status msg success\n";
            } else {
                int chunk_size = ps->get_chunk_size();
                for(int i = 0;i < chunk_size; i++)
                    comm->os<<rayList[i]->primary_size()<<" "<<rayList[i]->secondary_size()<<" "<<rayList[i]->type<<" \n";
                
                comm->os<<"need another chunk or need send rays\n";
            }
        } else {
          //  comm->os<<"mthread set proc busy\n";
            ps->set_proc_busy(comm->rank);
        }
        
        if (ps->has_new_chunk()) {
            inList_not_empty.notify_all();
        }
        if (!outList_empty()) {
            comm->os<<"mthread outlist empty"<<outList_empty()<<"\n";
            out_mutex.lock();
            int cId = get_sent_list();
            RayMsg* ray_msg;
            if(cId >= 0 && rayList[cId]->size() > 0.4 * ps->get_buffer_size()) {
                int dst_proc = ps->get_proc(cId);
                
                comm->os<<"mthread new RayMsg "<<cId<<" dst proc "<<dst_proc<<" size "<<rayList[cId]->size()<<" local chunk "<<ps->get_local_chunk()<<"\n";
                ps->set_proc_busy(dst_proc);
                ray_msg = new RayMsg(rayList[cId], comm->rank, dst_proc, cId, false); 
                comm->send_message(ray_msg, ps);
            } else {
                for(int i = 0; i < ps->get_chunk_size(); i++) {
                    if(rayList[i] -> size() > 0 && rayList[i]->type == "out") {
                        ray_msg = new RayMsg(rayList[i], comm->rank, comm->master, i, false); 
                        comm->send_message(ray_msg, ps);
                        break;
                    }
                }
            }
            out_mutex.unlock(); 
        }
    }

    comm->purge_completed_mpi_buffers();
}

void MWNode::message_thread(void* tmp) {
  
    MWNode *wk = (struct MWNode*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus * ps = wk->ps;

    printf("%d message thread\n", comm->rank);
    comm->os<<"message thread\n";
    int recv_count = 0;
    while (!ps->Exit()) {
       
        bool recv = false;
        if(comm->isMaster()){ 
            printf("master before recv thread wait %d inlist %d %d %d %d\n", 
               ps->all_thread_waiting(), wk->inList->size(), ps->all_proc_idle(), ps->all_rays_received(), wk->rayList_empty());
        }
        do {
      //      if(comm->isMaster()) 
      //          comm->os<<"master wait recv\n";
            recv = comm->recv_message(wk->rayList, wk->inList,  ps);
        } while(!recv && ps->is_proc_idle() && !ps->Exit()) ;
        
        if(ps->Exit()) break;
        
        //if new chunk havent loaded, render thread havent start, mthread wait
        std::unique_lock <std::mutex> lock(wk->thread_mutex); 
        if(!wk->inList->empty()) {
            //原有的光线 先放自己这 剩余的 发给master
            if(ps->has_new_chunk() ) {
                comm->os<<"wait new chunk\n";
                while( ps->has_new_chunk()) {
                    comm->os<<"wait rthread launch\n";

                    wk->render_start.wait(lock); 
                }
            } 
            wk->inList_not_empty.notify_one();
        }
        lock.unlock(); 

        if(ps->Exit()) break;
        
        wk->send_message(); 
        
        usleep(100);
	}
    comm->os <<" end message thread"<<ps->all_thread_waiting()<<"\n";
    comm->os <<" inlist "<< wk->inList->size()
             <<" recv "<<ps->global_rays[comm->rank + comm->size]
             <<std::endl;
    wk->inList_not_empty.notify_all();
    return;
} 

void MWNode::run(ImageDecomposition * camera) {
    
    comm->os <<" start run message thread \n";
    
    camera->decomposition(ps->get_chunk_map(), true, comm->rank, comm->size); 
    ps->updata_local_chunk();

    int deviceNum = ps->get_dev_num();
    int iter = 0; 
    std::vector<std::thread> workThread;
    if(comm->size == 1) {
        printf("only one worker %d  ", comm->rank);
        for(int i = 0; i < deviceNum; i++) 
            workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, true));
        
        for( auto &thread: workThread) 
            thread.join();
    	
    } else {
        comm->os <<"mthread  \n";
        printf("%d start mthread run\n ", comm->rank);
        std::thread mthread(message_thread, this);
        //    all proc start proc 0 send schedule? 
        while(1) {
            printf("new chunk raylist local chunk %d ", ps->get_local_chunk());
            for(int i = 0; i < ps->get_chunk_size(); i++)
                printf(" %d ", rayList[i]->size());
            printf("\n");

            comm->os<<"                                   new chunk raylist left chunk " << ps->get_local_chunk()<<" ";
            for(int i = 0; i < ps->get_chunk_size(); i++)
                comm->os<<" "<<rayList[i]->size();
            comm->os<<"\n";

            if(ps->has_new_chunk()) {
                ps->chunk_loaded();

                printf("%d before set render start \n", comm->rank);
                for(int i = 0; i < deviceNum; i++) 
                    workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, iter == 0));
                 
                //set chunk loaded
                printf("%d set render start \n", comm->rank);
                comm->os<<" set render start \n";
                
                render_start.notify_all(); 
                comm->os<<" render start notify \n";
                printf("%d render start notify \n", comm->rank);

                for(auto &thread: workThread)
                    thread.join();

                workThread.clear();
            } else {
                ps->set_exit();
                break;
            }
            iter++;
        }
        printf("%d worker exit\n", comm->rank);
        mthread.join();
    }
    return;
    printf("run MWNode\n");
}

