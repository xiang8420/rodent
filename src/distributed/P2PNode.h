
int32_t * get_prerender_result();

// Asynchronous p2p model without master

class P2PNode : public Node {

protected:    
    int * statistic;
    
   
public:    
    P2PNode(struct Communicator *comm, struct ProcStatus *ps);

    ~P2PNode();

    bool check_rendering_status(); 
   
    // send primary and secondary
    void count_rays(int, int);

    static void message_thread(void* tmp);
    
    void run(ImageDecomposition * camera);

}; 

P2PNode::P2PNode(struct Communicator *comm, struct ProcStatus *ps) 
    : Node(comm, ps) 
{

    int chunk_size = ps->get_chunk_size();
    assert(chunk_size == comm->get_size()); 
    
    statistic = new int[chunk_size * 4/*dep*/ * 2]; 

    std::fill(statistic, statistic + chunk_size * 4/*dep*/ * 2, 0);
    
    comm->os << "out list size"<<rayList[0].size() << "capacity" << rayList[0].get_primary().get_capacity() << std::endl; 
}

P2PNode::~P2PNode() {
    delete[] statistic;
}

//check proc status, return if proc need to wait 
bool P2PNode::check_rendering_status() {
    comm->os<<"mthread status inlist size "<<inList.size()<<" thread wait "<<ps->all_thread_waiting()<<"\n";
    for(int i = 0; i < ps->get_chunk_size(); i++)
        comm->os<<" "<<rayList[i].size();
    comm->os<<"\n";

    if(ps->all_thread_waiting() && inList.empty() ) {
        if(rayList_empty()) {
            comm->os<< "mthread if all thread idle all queue empty, set itself idle\n";
            ps->set_proc_idle(comm->get_rank());
            if(ps->all_proc_idle() && ps->all_rays_received()) {
                int *s = ps->get_status();

                //if all worker end synchronous 
                comm->os<<"mthread all proc idle send collective\n";
                QuitMsg quit_msg(comm->get_rank(), get_tag()); 
                comm->send_message(&quit_msg, ps);
                ps->set_exit();
                comm->os<<"set exit\n";
                inList_not_empty.notify_all();
            } else {
                comm->os<<"mthread proc idle send status\n";
                StatusMsg status(comm->get_rank(), -1, ps->get_status(), ps->get_local_chunk(), comm->get_size(), get_tag()); 
                comm->send_message(&status, ps);
            }
            return true; 
        } else {
            //maybe load new chunk
            comm->os<<"need another chunk or need send rays\n";
            ps->set_proc_busy(comm->get_rank());
            return false; 
        }
    } else {
        comm->os<<"mthread proc busy\n";
        ps->set_proc_busy(comm->get_rank());
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
            recv = comm->recv_message(wk->rayList.data(), &(wk->inList), ps);
        } while(!recv && ps->is_proc_idle() && !ps->Exit()) ;

        if(!wk->inList.empty()) {
            wk->inList_not_empty.notify_one();
        }

        if(ps->Exit()) break;

        int cId = wk->get_sent_list();
        
        if(cId >= 0 ) {
            ps->set_proc_busy(ps->get_proc(cId));
            wk->out_mutex.lock();
            comm->os<<"mthread new RayMsg"<<cId<<"size "<<wk->rayList[cId].size()<<"\n";
            RayMsg ray_msg(wk->rayList[cId], comm->get_rank(), ps->get_proc(cId), cId, false, wk->get_tag()); 
            comm->os<<"mthread RayMsg"<<ray_msg.get_chunk()<<"\n";
            wk->out_mutex.unlock(); 
            comm->send_message(&ray_msg, ps);
        }
        comm->purge_completed_mpi_buffers();
        usleep(100);
	}
    comm->os <<" end message thread"<<ps->all_thread_waiting()<<"\n";
    comm->os <<" inlist "<< wk->inList.size()
             <<" recv "<<ps->global_rays[comm->get_rank() + comm->get_size()]
             <<std::endl;
    wk->inList_not_empty.notify_all();
    return;
} 

void P2PNode::count_rays(int width, int height) 
{

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

void P2PNode::run(ImageDecomposition * camera) {
    

    int deviceNum = ps->get_dev_num();
//    if(comm->get_rank() == 0) {
//        work_thread(this, camera, 0, 1, true, true);
//        count_rays(camera->width, camera->height); 
//    }
    
    std::vector<std::thread> workThread;
    std::thread mthread(message_thread, this);
    comm->os << "loaded chunk()"<<ps->get_local_chunk()<< std::endl; 
    ps->chunk_loaded();
    //    //all proc start proc 0 send schedule? 
    while(!ps->Exit()) {
        int chunk = ps->get_local_chunk(); 
        for(int i = 0; i < deviceNum; i++) 
            workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, true));
        
        for( auto &thread: workThread) 
            thread.join();
        
    }
    printf("%d worker exit\n", comm->get_rank());
    mthread.join();
    return;
} 
