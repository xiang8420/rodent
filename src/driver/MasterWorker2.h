
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
    
    void run(float* rProcessTime); 
    
    void save_outgoing_buffer(float *rays, size_t size, size_t capacity, bool primary);

    void clear_outlist() {  
        int chunk_size = ps->get_chunk_size();
        for(int i = 0; i < chunk_size; i++) {
            if(rayList[i]->type == "out" ) {
                rayList[i]->primary->size = 0;
                rayList[i]->secondary->size = 0;
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
   
    comm->os<<"master set up\n";
    rayList       = new RayList *[chunk_size];
    for(int i = 0; i < chunk_size; i++) {
        rayList[i] = new RayList(4 * ps->get_buffer_size(), "out", false);
    }
    inList  = rayList[ps->get_local_chunk()];
    inList->type = "in";
    buffer = new RayList(ps->get_buffer_size(), "buffer", false); 
    printf("wmnode rank %d\n", comm->rank); 
    
    st = clock();
    recv_t = 0; wait_t = 0;
    send_t = 0; total_t = 0;
    write_t = 0; read_t = 0;
}

MWNode::~MWNode(){
    for(int i = 0; i < ps->get_chunk_size(); i++){
        delete rayList[i];
    }
    delete buffer;
    delete ps;
    delete comm;
}

void MWNode::save_outgoing_buffer(float *retired_rays, size_t size, size_t capacity, bool primary){
    comm->os<<"rthread save outgoing buffer "<< size <<"\n"; 
    out_mutex.lock(); 
    renderer_save_rays += size;
    int width = primary?21:14; 
    int* ids = (int*)(retired_rays);
    for(int i = 0; i < 5; i ++) {
        comm->os<<"| "<< ids[i] <<" "<< ids[i + ps->get_buffer_capacity() * 9] << " ";
    }
    comm->os<<"\n";
    RayList::read_from_device_buffer(rayList, retired_rays, size, capacity, primary, comm->rank);
    out_mutex.unlock(); 
}

int MWNode::get_unloaded_chunk(){
    int unloaded_chunk = -1;
    for(int i = 0; i < ps->get_chunk_size(); i++) {
        if(rayList[i]->size() > 0 && ps->get_target_proc(i) == -1) {
            unloaded_chunk = i; break; 
        }
    }
    return unloaded_chunk;
}

void MWNode::check_proc_status() {
//    if(comm->isMaster()) {
//        //master idle
//        comm->os<< "mthread master check proc status\n";
//        printf( "mthread master check proc status\n");
//        if(ps->all_thread_waiting() && buffer->empty() && rayList_empty()) {
//            comm->os<< "mthread if all thread idle all queue empty, set itself idle\n";
//            ps->set_proc_idle(comm->rank);
//            if(ps->all_proc_idle() && ps->all_rays_received()) {
//                QuitMsg quit_msg(comm->rank); 
//                comm->send_message(&quit_msg, ps);
//                ps->set_exit();
//                comm->os<<"set exit\n";
//                buffer_not_empty.notify_all();
//            }
//        }
//        //proc idle
//        int idle_proc = ps->get_idle_proc();
//        if(idle_proc >= 0) {
//            //find a unloaded chunk contents unprocessed rays
//            // if idle proc loaded chunk not empty ???
//            //
//            int unloaded_chunk = get_unloaded_chunk(); 
//            //find a idle proc send chunk id to it
//            if(unloaded_chunk > -1) {
//                comm->os<<"mthread scheduleMsg"<<"\n";
//                ps->update_chunk(idle_proc, unloaded_chunk);
//                int * a = ps->get_chunk_map();
//                comm->os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";
//
//                ScheduleMsg *msg = new ScheduleMsg(comm->rank, ps->get_chunk_map(), idle_proc, ps->get_chunk_size()); 
//                comm->send_message(msg, ps);
//                ps->set_proc_busy(idle_proc);
//            }  
//        }
//       // printf("master after proc status\n");
//    } else {
//      //  comm->os<<"mthread status buffer size "<<buffer->size()<<" thread wait "<<ps->all_thread_waiting()<<"\n";
//        if(ps->all_thread_waiting() && buffer->empty() ) {
//            if(rayList_empty()) {
//                comm->os<< "mthread send status msg\n";
//                ps->set_proc_idle(comm->rank);
//                //send to master I am idle
//                StatusMsg status(comm->rank, comm->master, ps->get_status(), ps->get_local_chunk(), comm->size); 
//                comm->send_message(&status, ps);
//            } else {
//                int chunk_size = ps->get_chunk_size();
//                for(int i = 0;i < chunk_size; i++)
//                    comm->os<<rayList[i]->primary_size()<<" "<<rayList[i]->secondary_size()<<" "<<rayList[i]->type<<" \n";
//                
//                int * a = ps->get_chunk_map();
//                comm->os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";
//                
//                comm->os<<"need another chunk or need send rays\n";
//            }
//        } else {
//            ps->set_proc_busy(comm->rank);
//        }
//        
//        if(ps->has_new_chunk()) {
//            buffer_not_empty.notify_all();
//        }
//    }
} 

void MWNode::send_message() {
    if(comm->isMaster()) {
        //master idle
        printf("master all status all thread wait %d inlist %d %d %d %d\n", 
                ps->all_thread_waiting(), inList->size(), ps->all_proc_idle(), ps->all_rays_received(), rayList_empty());
        int * a = ps->get_chunk_map();
        printf("%d %d %d %d \n", a[0], a[1], a[2], a[3]);
        printf("proc idle %d %d\n", ps->proc_idle[0], ps->proc_idle[1] );
        printf("master ray list %d %d %d %d \n", rayList[0]->size(), rayList[1]->size(), rayList[2]->size(), rayList[3]->size());

        if(ps->all_thread_waiting() && inList_empty() && rayList_empty()) {
            comm->os<< "mthread if all thread idle all queue empty, set itself idle\n";
            ps->set_proc_idle(comm->rank);
        }

        if(ps->all_proc_idle() && ps->all_rays_received() && rayList_empty()) {
            QuitMsg quit_msg(comm->rank); 
            comm->send_message(&quit_msg, ps);
            ps->set_exit();
            comm->os<<"set exit\n";
            buffer_not_empty.notify_all();
            return;
        }

        //proc idle
        int idle_proc = ps->get_idle_proc();
        comm->os<<"mthread status idle proc ";
        for(int i = 0; i< comm->size;i++){
            comm->os<<ps->proc_idle[0] <<" "<<ps->global_rays[i] <<" "<<ps->global_rays[i + comm->size]<<" ||"; 
        }
        comm->os<<"\n";
            
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
                    int * a = ps->get_chunk_map();
                    comm->os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";

                    RayMsg *ray_msg = new RayMsg(rayList[unloaded_chunk], comm->rank, idle_proc, unloaded_chunk, false); 
                    
                    comm->os<<"mthread schedule RayMsg "<<ray_msg->get_chunk()<<"\n";
                    comm->send_message(ray_msg, ps);
                    ps->set_proc_busy(idle_proc);
                    
                    ScheduleMsg *msg = new ScheduleMsg(comm->rank, ps->get_chunk_map(), idle_proc, ps->get_chunk_size()); 
                    comm->send_message(msg, ps);
                    
              //      out_mutex.unlock(); 
              //      return;//////??
                }
            }
        }
        if(!outList_empty()) {
            int cId = get_sent_list(rayList, ps);
            if(cId >= 0 ) {
                int dst_proc = ps->get_target_proc(cId);
                comm->os<<"mthread master new RayMsg "<<cId<<" size "<<rayList[cId]->size()<<" rank "<<comm->rank<<" dst "<<dst_proc<<"\n";
                comm->os<<"mthread ray list  "<<rayList[0]->size()<<" "<<rayList[1]->size()<<" "<<rayList[2]->size()<<" "<<rayList[3]->size()<<"\n";
                printf("master mthread ray listt %d %d %d %d %d\n", cId, rayList[0]->size(), rayList[1]->size(), rayList[2]->size(), rayList[3]->size());
                int * a = ps->get_chunk_map();
                comm->os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";
                
                ps->set_proc_busy(dst_proc);
              //  out_mutex.lock();
                RayMsg *ray_msg = new RayMsg(rayList[cId], comm->rank, dst_proc, cId, false); 
              //  out_mutex.unlock(); 
                
                comm->os<<"mthread RayMsg "<<ray_msg->get_chunk()<<"\n";
                comm->send_message(ray_msg, ps);
            }
        }
        out_mutex.unlock(); 
       // printf("master after proc status\n");
    } else {
        comm->os<<"mthread status buffer size "<<buffer->size()<<" thread wait "<<ps->all_thread_waiting()<<"\n";
        if (ps->all_thread_waiting() && inList_empty()) {
            comm->os<<"mthread before empty raylist left "<<rayList[0]->size()<<" "<<rayList[1]->size()<<" "<<rayList[2]->size()<<" "<<rayList[3]->size()<<"\n";
            if(rayList_empty()) {
                comm->os<< "mthread send status msg\n";
                ps->set_proc_idle(comm->rank);
                //send to master I am idle
                comm->os<<"mthread local chunk"<<ps->get_local_chunk()<<"\n";
                comm->os<<"mthread before status raylist left "<<rayList[0]->size()<<" "<<rayList[1]->size()<<" "<<rayList[2]->size()<<" "<<rayList[3]->size()<<"\n";
                StatusMsg status(comm->rank, comm->master, ps->get_status(), ps->get_local_chunk(), comm->size); 
                comm->send_message(&status, ps);
                comm->os<< "mthread send status msg success\n";
            } else {
                int chunk_size = ps->get_chunk_size();
                for(int i = 0;i < chunk_size; i++)
                    comm->os<<rayList[i]->primary_size()<<" "<<rayList[i]->secondary_size()<<" "<<rayList[i]->type<<" \n";
                int * a = ps->get_chunk_map();
                comm->os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";
                
                comm->os<<"need another chunk or need send rays\n";
            }
        } else {
          //  comm->os<<"mthread set proc busy\n";
            ps->set_proc_busy(comm->rank);
        }
        
        if (ps->has_new_chunk()) {
            buffer_not_empty.notify_all();
        }
        if (!outList_empty()) {
            comm->os<<"mthread ray list  "<<rayList[0]->size()<<" "<<rayList[1]->size()<<" "<<rayList[2]->size()<<" "<<rayList[3]->size()<<"\n";
            out_mutex.lock();
            int cId = get_sent_list(rayList, ps);
            RayMsg* ray_msg;
            if(cId >= 0 ) {
                int dst_proc = ps->get_target_proc(cId);
                comm->os<<"mthread new RayMsg "<<cId<<" dst proc "<<dst_proc<<" size "<<rayList[cId]->size()<<" local chunk "<<ps->get_local_chunk()<<"\n";
                int * a = ps->get_chunk_map();
                comm->os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";
                comm->os<<"mthread "<<rayList[0]->type<<" "<<rayList[1]->type<<" "<<rayList[2]->type<<" "<<rayList[3]->type<<"\n";
                ps->set_proc_busy(dst_proc);
                ray_msg = new RayMsg(rayList[cId], comm->rank, dst_proc, cId, false); 
                comm->send_message(ray_msg, ps);
            } else {
                comm->os<<"mthread "<<cId<<"\n";
                for(int i = 0; i < ps->get_chunk_size(); i++) {
                    if(rayList[i] -> size() > 0 && rayList[i]->type == "out") {
                        int * a = ps->get_chunk_map();
                        comm->os<<"mthread new RayMsg unloaded "<< i <<"size "<<rayList[i]->size()<<"\n";
                        comm->os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";
                        comm->os<<"mthread "<<rayList[0]->type<<" "<<rayList[1]->type<<" "<<rayList[2]->type<<" "<<rayList[3]->type<<"\n";
                        comm->os<<"mthread "<<rayList[0]->size()<<"  "<<rayList[1]->size()<<"  "<<rayList[2]->size()<<" "<<rayList[3]->size()<<"\n";
                        ray_msg = new RayMsg(rayList[i], comm->rank, comm->master, i, false); 
                        comm->send_message(ray_msg, ps);
                        break;
                    }
                }
              //  comm->os<<"clear "<<rayList[0]->size()<<" "<<rayList[1]->size()<<" "<<rayList[2]->size()<<" "<<rayList[3]->size()<<"\n";
              //  comm->os<<"after clear  "<<rayList[0]->size()<<" "<<rayList[1]->size()<<" "<<rayList[2]->size()<<" "<<rayList[3]->size()<<"\n";
              //  comm->os<<"mthread send to master\n";
              //  ray_msg = new RayMsg(rayList, comm->rank, comm->master, ps->get_chunk_size(), false); 
              //  comm->send_message(ray_msg, ps);
            }
          //  clear_outlist(); 
            out_mutex.unlock(); 
        }
    }

    comm->purge_completed_mpi_buffers();
}

void MWNode::message_thread(void* tmp) {
  
    MWNode *wk = (struct MWNode*)tmp;
    Communicator * comm = wk->comm; 
    ProcStatus * ps = wk->ps;

    comm->os<<"message thread\n";
    int recv_count = 0;
    while (!ps->Exit()) {
       
        bool recv = false;
        if(comm->isMaster()){ 
            printf("master before recv thread wait %d inlist %d %d %d %d\n", 
                    ps->all_thread_waiting(), wk->inList->size(), ps->all_proc_idle(), ps->all_rays_received(), wk->rayList_empty());
        }
        do {
            if(comm->isMaster()) {
                printf("master wait  %d %d %d %d \n", wk->rayList[0]->size(), wk->rayList[1]->size(), wk->rayList[2]->size(), wk->rayList[3]->size());
        //        comm->os<<"master wait recv\n";
            }
            recv = comm->recv_message(wk->rayList, ps);
        } while(!recv && ps->is_proc_idle() && !ps->Exit()) ;
        
        if(ps->Exit()) break;
        
        //if new chunk havent loaded, render thread havent start, mthread wait
        std::unique_lock <std::mutex> lock(wk->thread_mutex); 
        if(ps->has_new_chunk() && !wk->buffer->empty()) {
            printf("ps has unloaded chunk");
            int * a = ps->get_chunk_map();
            comm->os<<a[0]<<" "<<a[1]<<" "<<a[2]<<" "<<a[3]<<"\n";
            //原有的光线 先放自己这 剩余的 发给master
            comm->os<<"wait new chunk\n";
            wk->buffer_not_empty.notify_all();
            while( ps->has_new_chunk()) {
                comm->os<<"wait rthread launch\n";
                wk->render_start.wait(lock); 
            }
        }
        lock.unlock(); 
//        comm->os <<"mthread after recv rays \n";
        
//        if(comm->isMaster()) 
//            comm->os<<"before write rya buffer "<< wk->rayList[0]->size()<<"  "<<wk->rayList[1]->size()<<"  "<<wk->rayList[2]->size()<<" "<<wk->rayList[3]->size()<<"\n";
        wk->write_rays_buffer();
        
        if(ps->Exit()) break;

        wk->send_message(); 
        
        if(comm->isMaster()) 
            printf("master after send\n");

        usleep(100);
	}
    comm->os <<" end message thread"<<ps->all_thread_waiting()<<"\n";
    comm->os <<" inlist "<< wk->inList->size()
             <<" buffersize"<<wk->buffer->size()
             <<" get render rays "<<wk->renderer_get_rays
             <<" save render rays" <<wk->renderer_save_rays
             <<" recv "<<ps->global_rays[ps->proc_rank + ps->proc_size]
             <<std::endl;
    wk->buffer_not_empty.notify_all();
    return;
} 

void MWNode::run(float * frame_time) {
    comm->os <<" start run message thread \n";
    int deviceNum = ps->get_dev_num();
    int iter = 0; 
    std::vector<std::thread> workThread;
    if(comm->size == 1) {
        std::cerr << "please use P2PNode for one node\n";
    } else {
        comm->os <<"mthread  \n";
        std::thread mthread(message_thread, this);
        //    all proc start proc 0 send schedule? 
        while(1) {
            comm->os<<"                                   new chunk raylist left chunk " << ps->get_local_chunk()
                 <<" |"<<rayList[0]->size()<<" "<<rayList[1]->size()<<" "<<rayList[2]->size()<<" "<<rayList[3]->size()<<"\n";

            if(ps->has_new_chunk()) {
                int chunk = ps->get_local_chunk(); 
                inList = rayList[chunk];
                for(int i = 0; i < ps->get_chunk_size(); i++) {
                    rayList[i]->type = "out";
                }
                inList->type = "in";
                ps->chunk_loaded();

                printf("%d before set render start \n", comm->rank);
                for(int i = 0; i < deviceNum; i++) 
                    workThread.emplace_back(std::thread(work_thread, this, frame_time, i, deviceNum, false, iter == 0));
                 
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
