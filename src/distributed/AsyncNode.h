#include <cmath>

// Original Master-Worker mode 
struct AsyncNode : public Node {

    AsyncNode(Communicator *, ProcStatus *, Scheduler *);
    
    ~AsyncNode();

    int get_sent_list(); 
    void send_message(); 
    void save_outgoing_buffer(float *, size_t, bool); 
    void copy_to_inlist(int current_chunk);
    RayMsg* export_ray_msg(int, int, int, bool, int); 
    static void message_thread(void* tmp);
    int load_incoming_buffer(float **, bool, int); 
    void process_chunk_hit(int*, int); 
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
            //return;
        } else if (outList_empty()) {
            statistics.start("run => message_thread => send_message => switch_chunk");
            ///load new chunk
            comm->os<<"before switch: ";
            for(int i = 0;i < chunk_size; i++) {
                comm->os<<outlist_comm[i].size()<<" ";
            }
            comm->os<<"\n";
            ps->switch_current_chunk(outlist_comm);
            
            int current_chunk = ps->get_current_chunk();
            comm->os<<"mthread copy new chunk " << current_chunk << "\n";
            copy_to_inlist(current_chunk);
                
            inlist_comm.cond_full.notify_all();
            inlist_render.cond_full.notify_all();
            
            statistics.end("run => message_thread => send_message => switch_chunk");
           // return;
        } else {
            //outlist need to send
        }
    }
    
    int cId = get_sent_list();
    do {
        if(cId >= 0) {
            comm->os<<"mthread outlist "<<cId<<" size  "<<outlist_comm[cId].size()<<"\n";
            int dst_proc = ps->get_proc(cId);
            comm->os<<"chunk "<<cId<<" get proc "<<dst_proc<<"\n";
            if(dst_proc >= 0) {
                ps->set_proc_busy(dst_proc);
                statistics.start("run => message_thread => send_message => new RayMsg");
                comm->os<<"mthread construct msg\n";
                RayMsg *ray_msg = export_ray_msg(cId, comm->get_rank(), dst_proc, false, comm->get_tag()); 
                comm->os<<"mthread get new msg\n";
                statistics.end("run => message_thread => send_message => new RayMsg");
                statistics.start("run => message_thread => send_message => comm->send");
                comm->send_message(ray_msg, ps);
                statistics.end("run => message_thread => send_message => comm->send");
                comm->os<<"send over\n";
                delete ray_msg;
            } else {
                error("ray msg");
            }
        }
        cId = get_sent_list();
    } while(cId >= 0 && ps->is_proc_idle());
}

void AsyncNode::save_outgoing_buffer(float *rays, size_t size, bool primary) {

    statistics.start("run => wthread => thread_save_outgoing");
    RetiredRays *retired_rays = new RetiredRays(rays, size, primary);
    statistics.end("run => wthread => thread_save_outgoing");

    statistics.start("run => wthread => lock_write_outgoing");
    std::lock_guard <std::mutex> thread_lock(retired_mutex); 
    
    scheduler->pass_record->write_send(rays, size, primary, scheduler->chunk_manager->local_chunks.current); 

    retiredRays.emplace_back(retired_rays);

    if(retiredRays.size() > 10)
        RetiredRays::clear_retired_rays(retiredRays, outlist_comm, chunk_size, comm->get_rank());

    statistics.end("run => wthread => lock_write_outgoing");
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
    //printf("copy primary size %d\n", copy_size);
    memcpy(*rays, rays_stream->get_data(), ps->get_stream_store_capacity() * width * sizeof(float)); 
  
    //scheduler->chunk_manager->local_chunks.recv_rays(copy_size);

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
        statistics.start("run => message_thread => sleep");
        usleep(600);
        statistics.end("run => message_thread => sleep");

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
    scheduler->preprocess(cam, comm->get_size(), false);
    
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
            statistics.start("run => wthread => wait ");
            std::unique_lock <std::mutex> inlock(inlist_comm.mutex); 
            
            if(!retiredRays.empty()) 
                RetiredRays::clear_retired_rays(retiredRays, outlist_comm, chunk_size, comm->get_rank());
            
            if(inlist_comm.empty()) {
                ps->set_proc_idle();
                while (inlist_comm.empty() && !ps->Exit()) {
                    comm->os<<"rthread waiting \n";
                    inlist_comm.cond_full.wait(inlock);
                }
                if(ps->Exit()) {
                    statistics.end("run => wthread => wait ");
                    break;
                }
            }
            RayStreamList::swap(inlist_comm, inlist_render);
            inlist_comm.empty_notify(); //tell mthread inlist size changed 
            inlist_render.empty_notify(); //tell mthread inlist size changed 
            ps->set_proc_busy(comm->get_rank());
            statistics.end("run => wthread => wait ");
        }
        statistics.start("run => wthread => render");
        launch_rodent_render(deviceNum, iter==0);
        float t = statistics.end("run => wthread => render");
        scheduler->pass_record->write_time(t);
        iter++;
    } while(!ps->Exit());
    comm->os<<" render thread times "<<iter<<"\n";
    mthread.join();        
    if(SIMPLE_TRACE  && cam->iter == 0)
        postprocess(cam->iter);
    return;
}

int median(int a, int b, int c){
    if(a > b) {
        if(b > c) return b;
        else if(a > c)  return c;
        else return a;
    } else {
        if(b < c) return b;
        else if(a > c)  return c;
        else return a;
    }
}

int average(int a, int b, int c) {
    return (a + b + c) / 3; 
}

void AsyncNode::process_chunk_hit(int* chunk_hit, int size) {
    int res = CHUNK_HIT_RES;
    int all_face_size = res * res * 6;
    int face_size = res * res;

    for(int i = 0; i < size; i++) {
        int lu = chunk_hit[i];
        //int lu = reduce_buffer[all_face_size * i + face_size * j + (res - 1 - v) * res + u];
        
        int gradient = 0;
        int res = 0;
        int lu_pre = (lu >> 28) & 0xF;
        int lu_cur = lu & 0xF;
        int lu_next;
        for(int k = 0; k < 8; k++) {
            int bit_next = (k == 7 ? 1 : k + 1) * 4;
            lu_next = (lu >> bit_next) & 0xF;
            
           // if(lu_cur >= 13) {
           //     res += lu_cur << (k * 4);
           // } else {
           //     
           //     gradient = std::abs(lu_next - lu_cur) + std::abs(lu_pre - lu_cur);
           //     //if(lu_next != 0 && lu_cur != 0)
           //     //    gradient += std::abs(lu_next - lu_cur);
           //     if(gradient > 30)   
           //         res += lu_cur << (k * 4);
           // }
           // res += ((lu_next + lu_cur + lu_pre) / 3) << (k*4);
            res += average(lu_next, lu_cur, lu_pre) << (k*4);
            lu_pre = lu_cur;
            lu_cur = lu_next;
        }
        // if(gradient < 100)
        //     chunk_hit[all_face_size * i + face_size * j + u * res + v] = 0;
        chunk_hit[i] = res;
    }
}

void AsyncNode::postprocess(int iter) {
    int * chunk_hit = rodent_get_chunk_hit();
    int size = CHUNK_HIT_RES * CHUNK_HIT_RES * 6 * ps->get_chunk_size();
    printf("start reduce %d\n", comm->get_rank());
    {
        statistics.start("run => light field => bcast ");
        comm->update_chunk_hit(chunk_hit, size);
        statistics.end("run => light field => bcast ");
    }
    
    //updata render light field
    process_chunk_hit(chunk_hit, size); 
    rodent_update_render_chunk_hit(chunk_hit, size);
    //reduce ctib 存在了reduce 里
    statistics.start("run => light field => save img  ");
    if(comm->get_rank() == 0 && VIS_CHUNK_HIT) { 
        int* render_chunk_hit = rodent_get_render_chunk_hit_data();
        save_image_ctrb(render_chunk_hit, ps->get_chunk_size(), iter);
    }
    statistics.end("run => light field => save img  ");
    
    //pass_record get chunk speed
    //scheduler->set_load_chunk_hit(); 
}
