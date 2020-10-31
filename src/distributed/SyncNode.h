// Synchronize Dynamic  
struct SyncNode : public Node{

    SyncNode(struct Communicator *comm, struct ProcStatus *ps);
    
    ~SyncNode();
    
    int get_unloaded_chunk(int *);
    
    bool all_queue_empty();

    void save_outgoing_buffer(float *rays, size_t size, bool primary);

    void send_message();

    void gather_rays(int, int, bool);

    void schedule(int * list_size); 
    
    void synchronize (); 

    void worker_synchronize (); 
    
    static void message_thread(void* tmp);
    
    void run(ImageDecomposition * camera);
   
    std::vector <RayList> rayList;
    
};
//每一轮光线全部给master然后 master发送调度消息给各个节点 各个节点读取新场景， 接受光线
SyncNode::SyncNode(struct Communicator *comm, struct ProcStatus *ps)
    :Node(comm, ps)
{
    sync = true;
    printf("new SyncNode\n");
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i++) {
        rayList.emplace_back(RayList());
    }
    inList.set_capacity(ps->get_stream_logic_capacity(), 
                        ps->get_stream_store_capacity());
  
}

bool SyncNode::all_queue_empty() { 
    for(int i = 0; i < ps->get_chunk_size(); i++) {
        if(!rayList[i].empty()) {return false;}
    }
    return true;
}

SyncNode::~SyncNode() {
    printf("delete SyncNode\n");
}

void get_raylist_size(std::vector <RayList> &rayList, std::vector <int> &list_size) { 
    for(int i = 0; i < rayList.size(); i ++)
        list_size[i] = rayList[i].size();
}

void SyncNode::save_outgoing_buffer(float *retired_rays, size_t size, bool primary){
    out_mutex.lock(); 
    int width = primary?21:14; 
    comm->os<<"rthread save outgoing buffer "<< size <<" width "<<width<<"\n"; 
    int* ids = (int*)(retired_rays);
    for(int i = 0; i < 5; i ++) {
        comm->os<<"| "<< ids[i * width] <<" "<< ids[i * width + 9] << " ";
    }
    comm->os<<"\n";
//    for(int i = 0; i < ps->get_chunk_size(); i++) { 
//        printf("save out going ray data size %ld capacity %ld\n",rayList[i].secondary.data.size(), rayList[i].secondary.data.capacity());
//        printf("save out going ray data size %ld capacity %ld\n",rayList[i].primary.data.size(), rayList[i].primary.data.capacity());
//    }
    //RayStreamList::read_from_device_buffer(rayList.data(), retired_rays, size, primary, comm->get_rank());
    RayList::read_from_device_buffer(rayList.data(), retired_rays, size, primary, comm->get_rank(), ps->get_chunk_size());
    out_mutex.unlock(); 
//    for(int i = 0; i < ps->get_chunk_size(); i++)
//        printf("Node.h raylist primary size %d width %d : ", rayList[i].get_primary().get_size(), rayList[i].get_primary().get_store_width());
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

        inList.read_from_ptr(recv_buffer.data(), 0, all_rays_size, primary, owner);
    }

}

int SyncNode::get_unloaded_chunk(int* list_size) {
    int unloaded_chunk = -1;
    int max_ray_size = 0;
    int chunk_size = ps->get_chunk_size();
  
    int simple_mesh = chunk_size - 1; 
    if(ps->get_rough_trace() && list_size[simple_mesh] > max_ray_size ) { 
        unloaded_chunk = simple_mesh;
        max_ray_size   = list_size[simple_mesh] / 2;
    }

    for(int i = 0; i < chunk_size - 1; i++) {
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
    printf("proc_size %d > chunk size %d\n", proc_size, chunk_size);
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
    printf("schedule : ");
    for(int i = 0; i < chunk_size; i++) {
        printf("|c %d p %d ", i, ps->get_proc(i));
    }printf("\n");
}

void SyncNode::synchronize () {
  
    statistics.start("run => synchronize => MPI_Allreduce");
    printf("master sync raylist size %d\n", rayList.size());
    int chunk_size =  ps->get_chunk_size();
    std::vector<int> list_size(chunk_size); 

    get_raylist_size(rayList, list_size);
    MPI_Allreduce(list_size.data(), list_size.data(), chunk_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
    
    statistics.end("run => synchronize => MPI_Allreduce");

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
    comm->os<<"proc raylist : ";
    for(int chunk = 0; chunk < chunk_size; chunk++) {
        comm->os<<"( "<<chunk<<" "<<rayList[chunk].size()<<" ) | ";
    }
    comm->os<<"\n";

    statistics.start("run => synchronize => gather_rays");
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
    statistics.end("run => synchronize => gather_rays");
    printf("rthread proc %d after gather inlist size %d\n", comm->get_rank(), inList.size());
}

void SyncNode::run(ImageDecomposition * camera) {
    
    comm->os <<" start run message thread \n";

    int deviceNum = ps->get_dev_num();
    int iter = 0; 

    std::vector<std::thread> workThread;
    while(!ps->Exit()) {
        comm->os<<"run while\n";
        ps->chunk_loaded();

        std::cout<<"rthread rank " << comm->get_rank()<<" generate rays "<<ps->generate_rays()<<"\n";
        for(int i = 0; i < deviceNum; i++)
            workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, iter == 0 && ps->generate_rays()));
        
        render_start.notify_all(); 
        for(auto &thread: workThread)
            thread.join();
    
        comm->os<<"rthread end thread\n";
        workThread.clear();
        statistics.start("run => synchronize");
        synchronize();
        statistics.end("run => synchronize");
        iter++;
    }
    //

    int chunk_size = ps->get_chunk_size();
    printf("proc %d left rays :", comm->get_rank());
    for(int i = 0; i < chunk_size; i++) 
        printf(" %d || ", rayList[i].size()); 
    printf("\n");

    std::cout<<" chunk size "<<chunk_size<<" ray left "<<rayList[chunk_size - 1].size()<<"\n";
//    if(rayList[chunk_size - 1].empty()) {
//        //load rayList[chunk_size - 1] to liststream
//        ps->update_chunk(comm->get_rank(), chunk_size - 1);
//        RaysArray &primary   = rayList[chunk_size - 1].get_primary(); 
//        RaysArray &secondary = rayList[chunk_size - 1].get_secondary(); 
//
//        inList.read_from_ptr((char*)primary.get_data(), primary.get_size(), true, comm->get_rank());
//        inList.read_from_ptr((char*)secondary.get_data(), secondary.get_size(), false, comm->get_rank());
//        rayList[chunk_size - 1].clear();
//        
//        for(int i = 0; i < deviceNum; i++)
//            workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, false/*no new rays*/));
//        
//        for(auto &thread: workThread)
//            thread.join();
//         
//    }
//    assert(chunk_size % 2 == 1 || );
    if(chunk_size % 2 == 1 && rayList[chunk_size].size() > 0 ) {
        comm->os<<"proc "<<comm->get_rank()<<" bounce ray left "<<rayList[chunk_size - 1].size()<<"\n";      
    }
    
    return;
}

