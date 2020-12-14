// Synchronize Dynamic  
struct SyncNode : public Node{

    SyncNode(struct Communicator *comm, struct ProcStatus *ps);
    
    ~SyncNode();
    
    int get_unloaded_chunk(int *);
    
    int load_incoming_buffer(float **, bool, int); 
    
    void save_outgoing_buffer(float *rays, size_t size, bool primary);
    
    void send_message();

    void gather_rays(int, int, bool);
    
    void get_raylist_size(RayArrayList *, std::vector <int> &); 

    void schedule(int * list_size); 
    
    void synchronize (); 

    void worker_synchronize (); 
    
    static void message_thread(void* tmp);
    
    void run(Scheduler * camera);
   
    RayArrayList * outArrayList;
    RayStreamList inList;  
};
//每一轮光线全部给master然后 master发送调度消息给各个节点 各个节点读取新场景， 接受光线
SyncNode::SyncNode(struct Communicator *comm, struct ProcStatus *ps)
    :Node(comm, ps)
{
    int chunk_size = ps->get_chunk_size();
    outArrayList = new RayArrayList[chunk_size];
}

SyncNode::~SyncNode() {
    delete[] outArrayList;
    printf("delete SyncNode\n");
}

void SyncNode::save_outgoing_buffer(float *rays, size_t size, bool primary) {
    int width = primary ? PRIMARY_WIDTH : SECONDARY_WIDTH; 
    RayArrayList::read_from_device_buffer(outArrayList, rays, size, primary, ps->get_chunk_size());
}

void SyncNode::get_raylist_size(RayArrayList *outArrayList, std::vector <int> &list_size) { 
    int chunk_size = ps->get_chunk_size();
    for(int i = 0; i < chunk_size; i ++)
        list_size[i] = outArrayList[i].size();
}

void SyncNode::gather_rays(int chunk, int owner, bool primary) {

    std::vector<int> gather_rays_size(comm->get_size());
    struct RaysArray &rays = primary ? outArrayList[chunk].get_primary() : outArrayList[chunk].get_secondary();
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
    chunk_size = ps->get_rough_trace() ? chunk_size - 1 : chunk_size;
 
    if(ps->get_rough_trace()) { 
        int simple_mesh = chunk_size - 1; 
        if(ps->get_rough_trace() && list_size[simple_mesh] / 2 > max_ray_size ) { 
            unloaded_chunk = simple_mesh;
            max_ray_size   = list_size[simple_mesh] / 2;
        }
    }
    for(int i = 0; i < chunk_size; i++) {
        printf("get_unloaded_chunk %d %d\n", i, list_size[i]);
        if(list_size[i] > max_ray_size && ps->get_proc(i) == -1) {
            unloaded_chunk = i;  
            max_ray_size = list_size[i];
        }
    }

    return unloaded_chunk;
}

int SyncNode::load_incoming_buffer(float **rays, bool primary, int thread_id) {
    printf("load incoming buffer\n");
    comm->os<<"rthread load incoming buffer inlist size "<<inList.size()<<"\n";
    if(inList.empty()) { return -1; } 
    if(!inList.empty() && ps->Exit())
        warn("inlist not empty but prepared to exit\n");


    struct RaysStream *rays_stream;
    if(primary && inList.primary_size() > 0) {
        if(inList.secondary_size() > inList.primary_size())
            return 0;
        rays_stream = inList.get_primary();
    } else if (!primary && inList.secondary_size() > 0) {
        rays_stream = inList.get_secondary();
    } else {
        return 0;
    }

    statistics.start("run => work_thread => load_incoming_buffer-copy");

    ps->set_thread_idle(thread_id, false);
    int copy_size = rays_stream->size;
    comm->os<<"all rays chunk "<<ps->get_current_chunk()<<" size "<<copy_size<<"\n";
    int width = rays_stream->width;
    printf("copy primary size %d\n", copy_size);
 //   comm->os<<"rthread width "<<width <<" logic width "<<rays_stream->logic_width<<"\n";
    memcpy(*rays, rays_stream->get_data(), ps->get_stream_store_capacity() * width * sizeof(float)); 
 
//     //printf("%d rays_stream size %d %d %d\n", comm->get_rank(), rays_stream->size, primary, rays_stream->store_width);
    int* ids = (int*)(*rays);
    //int psize = std::min(5, copy_size);
    int psize = copy_size;
    
    
//    if(primary) {
//        comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
//        for(int i = 0; i < psize; i ++) {
//        //    comm->os<<"#" <<ids[i + rays_size]<<" "<<ids[i + rays_size + ps->get_stream_store_capacity() * 9]
//        //                                      <<" "<<ids[i + rays_size + ps->get_stream_store_capacity() * 10];
//            if((0xFF & ids[i + rays_size + ps->get_stream_store_capacity() * 10]) > 4 )
//                comm->os<<"pri pri chk "<<(0xFF & ids[i + rays_size + ps->get_stream_store_capacity() * 10])<<" read error\n";
//        }
//        comm->os<<"\n";
//    } else {
//        comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
//        for(int i = 0; i < psize; i ++) {
//         //   comm->os<<"#" <<ids[i + rays_size]<<" "<<ids[i + rays_size + ps->get_stream_store_capacity() * 9]
//         //                                     <<" "<<ids[i + rays_size + ps->get_stream_store_capacity() * 10]
//         //                                     <<" "<<ids[i + rays_size + ps->get_stream_store_capacity() * 11];
//            if((0xFF & ids[i + rays_size + ps->get_stream_store_capacity() * 11]) > 4 )
//                comm->os<<"sec pri chk "<<(0xFF & ids[i + rays_size + ps->get_stream_store_capacity() * 11])<<" read error\n";
//        }
//        comm->os<<"\n";
//    }
    delete rays_stream;

    statistics.end("run => work_thread => load_incoming_buffer-copy");
    return copy_size;
        
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
    int* chunk_proc =  ps->get_chunk_proc();
    for(int i = 0; i < chunk_size; i++)
        printf("%d %d |", chunk_proc[i*2], chunk_proc[i*2+1]);
    printf("\n");
}

void SyncNode::synchronize () {
  
    statistics.start("run => synchronize => MPI_Allreduce");
    int chunk_size =  ps->get_chunk_size();
    std::vector<int> list_size(chunk_size); 

    get_raylist_size(outArrayList, list_size);
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
    comm->os<<"proc "<<comm->get_rank()<<" chunk "<<ps->get_current_chunk()<<"\n";
    //如果proc空了 则给一个chunk并进行gather
    //implement get dst chunk, add proc_chunk_map
    comm->os<<"proc raylist : ";
    for(int chunk = 0; chunk < chunk_size; chunk++) {
        comm->os<<"( "<<chunk<<" "<<outArrayList[chunk].size()<<" ) | ";
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
        outArrayList[chunk].clear();
    }  
    statistics.end("run => synchronize => gather_rays");
    printf("rthread proc %d after gather inlist size %d\n", comm->get_rank(), inList.size());
}

void SyncNode::run(Scheduler * camera) {
    
    comm->os <<" start run message thread \n";

    int deviceNum = ps->get_dev_num();
    int iter = 0; 
    std::vector<std::thread> workThread;
    while(!ps->Exit()) {
        comm->os<<"run while\n";
        ps->chunk_loaded();

        int current_chunk = ps->get_current_chunk();
        if(outArrayList[current_chunk].size() > 0) {
            RaysArray &primary   = outArrayList[current_chunk].get_primary(); 
            RaysArray &secondary = outArrayList[current_chunk].get_secondary(); 
            
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
            outArrayList[current_chunk].clear();
            //read to inList
        //    clear_outlist();
        }

        launch_rodent_render(camera, deviceNum, iter==0);
        std::cout<<"rthread rank " << comm->get_rank()<<" generate rays "<<ps->generate_rays()<<"\n";
        statistics.start("run => synchronize");
        synchronize();
        statistics.end("run => synchronize");
        iter++;
    }
    //

    int chunk_size = ps->get_chunk_size();
    printf("proc %d left rays :", comm->get_rank());
    for(int i = 0; i < chunk_size; i++) 
        printf(" %d || ", outArrayList[i].size()); 
    printf("\n");

    std::cout<<" chunk size "<<chunk_size<<" ray left "<<outArrayList[chunk_size - 1].size()<<"\n";
//    if(outArrayList[chunk_size - 1].empty()) {
//        //load outArrayList[chunk_size - 1] to liststream
//        ps->update_chunk(comm->get_rank(), chunk_size - 1);
//        RaysArray &primary   = outArrayList[chunk_size - 1].get_primary(); 
//        RaysArray &secondary = outArrayList[chunk_size - 1].get_secondary(); 
//
//        inList.read_from_ptr((char*)primary.get_data(), primary.get_size(), true, comm->get_rank());
//        inList.read_from_ptr((char*)secondary.get_data(), secondary.get_size(), false, comm->get_rank());
//        outArrayList[chunk_size - 1].clear();
//        
//        for(int i = 0; i < deviceNum; i++)
//            workThread.emplace_back(std::thread(work_thread, this, camera, i, deviceNum, false, false/*no new rays*/));
//        
//        for(auto &thread: workThread)
//            thread.join();
//         
//    }
//    assert(chunk_size % 2 == 1 || );
    if(chunk_size % 2 == 1 && outArrayList[chunk_size].size() > 0 ) {
        comm->os<<"proc "<<comm->get_rank()<<" bounce ray left "<<outArrayList[chunk_size - 1].size()<<"\n";      
    }
    
    return;
}

