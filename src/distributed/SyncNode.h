// Synchronize Dynamic  
struct SyncNode : public Node{

    double recv_t, wait_t, send_t, total_t, read_t, write_t, st, ed;

    SyncNode(struct Communicator *comm, struct ProcStatus *ps);
    
    ~SyncNode();
    
    int get_unloaded_chunk(int *);
    
    void send_message();

    void gather_rays(int, int, bool);

    void schedule(int * list_size); 
    
    void synchronize (); 

    void worker_synchronize (); 
    
    static void message_thread(void* tmp);
    
    void run(ImageDecomposition * camera);
    
    int load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait);
    
};
//每一轮光线全部给master然后 master发送调度消息给各个节点 各个节点读取新场景， 接受光线
SyncNode::SyncNode(struct Communicator *comm, struct ProcStatus *ps)
    :Node(comm, ps)
{
    sync = true;
    printf("new SyncNode\n");
    int chunk_size = ps->get_chunk_size();
//    comm->os<< "master chunk size"<<chunk_size<<"\n";
  
    size_t memory[4];
    memory[0] = physical_memory_used_by_process();
    comm->os<<"master set up\n";
    
    memory[1] = physical_memory_used_by_process();
    printf("Memory %ld kb in list %ld\n", memory[0], memory[1] - memory[0]);
    
    st = clock();
    recv_t = 0; wait_t = 0;
    send_t = 0; total_t = 0;
    write_t = 0; read_t = 0;
}

SyncNode::~SyncNode() {
    printf("delete SyncNode\n");
}

void get_raylist_size(std::vector <RayList> &rayList, std::vector <int> &list_size) { 
    for(int i = 0; i < rayList.size(); i ++)
        list_size[i] = rayList[i].size();
}

void SyncNode::gather_rays(int chunk, int owner, bool primary) {

    std::vector<int> gather_rays_size(comm->get_size());
    struct RaysArray &rays = primary ? rayList[chunk].get_primary() : rayList[chunk].get_secondary();
    int width = rays.get_store_width();
    int list_size = rays.get_size();
    int proc_size = comm->get_size();
    MPI_Gather(&list_size, 1, MPI_INT, gather_rays_size.data(), 1, MPI_INT, owner, MPI_COMM_WORLD);

    std::vector<int> rcounts(proc_size);
    std::vector<int> displs(proc_size);
    std::vector<float> recv_buffer;
    int all_rays_size = 0;
    if(comm->get_rank() == owner) {
        comm->os<<"owner "<<owner<<" list size: ";
        for(int k = 0; k < comm->get_size(); k ++)
            comm->os<<gather_rays_size[k]<<" ";
        comm->os<<"\n";

        int in_size = width;
        for(int i = 0; i < proc_size; i++) {
            all_rays_size += gather_rays_size[i];
            rcounts[i] = gather_rays_size[i] * in_size; 
            displs[i] = i == 0 ? 0 : (displs[i - 1] + rcounts[i - 1]);
            comm->os<<"gather rays size "<<gather_rays_size[i]<<" rcounts " <<rcounts[i]<<" displs "<<displs[i]<<"\n";
        }
        recv_buffer.resize(all_rays_size * in_size);
    }
    {
        comm->os<<"before gatherv list size "<<list_size<<"\n";
        int* iptr = (int*) rays.get_data();
        for(int i = 0; i < 5; i ++) {
            comm->os<< iptr[i*width]<<" "<<iptr[i*width + 9]<<" | ";
        }
        comm->os<<"\n";
    }
    MPI_Gatherv((float*)rays.get_data(), list_size * width , MPI_FLOAT, (float*)recv_buffer.data(), rcounts.data(), displs.data(), MPI_FLOAT, owner, MPI_COMM_WORLD);
    if(comm->get_rank() == owner)  {
        comm->os<<"after gatherv owner ray size"<<all_rays_size<<"\n";
        int* iptr = (int*) recv_buffer.data();
        for(int i = 0; i < std::min(5, all_rays_size) ; i ++) {
            comm->os<< iptr[i*width]<<" "<<iptr[i*width + 9]<<" | ";
        }
        comm->os<<"\n";

        inList.read_from_ptr((char*)recv_buffer.data(), all_rays_size, primary, owner);
    }
     
    comm->os<<" raylist "<<chunk<<" size "<< rays.data.size()<<" capacity "<<rays.data.capacity()<<"\n";

}

int SyncNode::get_unloaded_chunk(int* list_size) {
    int unloaded_chunk = -1;
    int max_ray_size = 0;
    for(int i = 0; i < ps->get_chunk_size(); i++) {
        if(list_size[i] > max_ray_size && ps->get_proc(i) == -1) {
            unloaded_chunk = i;  
            max_ray_size = list_size[i];
        }
    }
    return unloaded_chunk;
}

void SyncNode::schedule(int * list_size) {
    int chunk_size = ps->get_chunk_size();
    int proc_size  = comm->get_size();

    if(proc_size >= chunk_size) 
        return;

    for(int i = 0; i < chunk_size; i++) {
         printf(" %d list size %d %d\n", comm->get_rank(), i, list_size[i]);
         //if proc list empty, assign a new chunk
         int chunk_proc = ps->get_proc(i);
         if(list_size[i] == 0 && chunk_proc >= 0) {
             printf("proc %d need new schedule\n ", comm->get_rank());
             int unloaded_chunk = get_unloaded_chunk(list_size);
             if(unloaded_chunk >= 0)
                 ps->update_chunk(chunk_proc, unloaded_chunk);
             //proc_[i] load new chunk 
         } 
    } 
}

void SyncNode::synchronize () {
  
    printf("master sync raylist size %d\n", rayList.size());
    int chunk_size =  ps->get_chunk_size();
    std::vector<int> list_size(chunk_size); 

    get_raylist_size(rayList, list_size);
    MPI_Allreduce(list_size.data(), list_size.data(), chunk_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
    
    //exit
    int all_rays = 0;
    for(int i = 0; i < chunk_size; i++) {
        all_rays += list_size[i];
    }
    if(all_rays == 0) {
        ps->set_exit(); 
        return;
    }
    
    schedule(list_size.data()); //静态调度时无用
    printf("proc %d chunk %d\n", comm->get_rank(), ps->get_local_chunk());
    //如果proc空了 则给一个chunk并进行gather
    //implement get dst chunk, add proc_chunk_map
    for(int chunk = 0; chunk < chunk_size; chunk++) {
        int owner = ps->get_proc(chunk);
        if(owner < 0) {
            //or send to master? 
            continue; 
        } 
        //gether all chunk rays 
        gather_rays(chunk, owner, true); 
        gather_rays(chunk, owner, false); 
        rayList[chunk].clear();
    }  
    printf("rthread after gather inlist size %d\n", inList.size());
    

//    for()
//
}

//void SyncNode::worker_synchronize () {
//    
//    int chunk_size =  ps->get_chunk_size();
//    std::vector<int> list_tmp_size(chunk_size); 
//    printf("worker sync\n");
//    for(int i = 0; i < rayList.size(); i ++)
//        printf("list tmp size %d\n", list_tmp_size[i]);
//
//    std::vector<int> list_size(chunk_size); 
//    comm->reduce((float*)list_tmp_size.data(), (float*)list_size.data(), chunk_size); 
////    for(dd) 
////    
////    if(recv exit)
////        ps->set_working();
//    
//}

int SyncNode::load_incoming_buffer(float **rays, size_t rays_size, bool primary, int thread_id, bool thread_wait) {
    printf("syn load incoming buffer\n");
    return -1;
}

void SyncNode::run(ImageDecomposition * camera) {
    
    comm->os <<" start run message thread \n";

    int deviceNum = ps->get_dev_num();
    int iter = 0; 

    std::vector<std::thread> workThread;
    while(!ps->Exit()) {
        comm->os<<"run while\n";
        ps->chunk_loaded();

        for(int i = 0; i < deviceNum; i++) 
            workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, iter == 0));
        
        render_start.notify_all(); 
        for(auto &thread: workThread)
            thread.join();
    
        comm->os<<"rthread end thread\n";

        workThread.clear();
        
        synchronize();
    
        iter++;

    }
    return;
}

