// Synchronize Dynamic  
struct SyncNode : public Node{

    SyncNode(Communicator *comm, ProcStatus *ps, Scheduler*);
    ~SyncNode();
    int get_unloaded_chunk(int *);
    int load_incoming_buffer(float **, bool, int); 
    void save_outgoing_buffer(float *rays, size_t size, bool primary);
    void send_message();
    void gather_rays(int, int);
    void schedule(int * list_size); 
    void synchronize (); 
    void worker_synchronize (); 
    static void message_thread(void* tmp);
    void run(Camera *);
   
    RayStreamList * outlist;
    RayStreamList inlist;  
    
    std::vector<RetiredRays*> retiredRays; 
    std::mutex retired_mutex;
    int chunk_size;
};
//每一轮光线全部给master然后 master发送调度消息给各个节点 各个节点读取新场景， 接受光线
SyncNode::SyncNode(struct Communicator *comm, struct ProcStatus *ps, Scheduler *scheduler)
    :Node(comm, ps, scheduler)
{
    chunk_size = ps->get_chunk_size();
    int store_capacity = ps->get_stream_store_capacity();
    int logic_capacity = ps->get_stream_logic_capacity();
    
    outlist = new RayStreamList[chunk_size];
    for(int i = 0; i < chunk_size; i++)
        outlist[i].set_capacity(logic_capacity, store_capacity);
    inlist.set_capacity(logic_capacity, store_capacity);
}

SyncNode::~SyncNode() {
    delete[] outlist;
    printf("delete SyncNode\n");
}

void SyncNode::save_outgoing_buffer(float *rays, size_t size, bool primary) {

    RetiredRays *retired_rays = new RetiredRays(rays, size, primary);
    
    printf("rank %d before save outgoing\n", comm->get_rank());
    std::lock_guard <std::mutex> lock(retired_mutex); 
    retiredRays.emplace_back(retired_rays);
   // printf("rank %d after save outgoing\n", comm->get_rank());
}

void SyncNode::gather_rays(int chunk, int owner) {
    RayMsg msg(outlist[chunk], comm->get_rank(), owner, chunk, false, 0); 
    int ray_size[3];
    ray_size[0] = msg.get_content_size();
    ray_size[1] = msg.primary_size();
    ray_size[2] = msg.secondary_size();
    
    printf("before gather rays %d %d %d chunk %d chunk size %d\n", ray_size[0], ray_size[1], ray_size[2], chunk, outlist[chunk].size());
//    printf("\n owner %d chunk %d | rank %d ray_size %d %d %d|\n", owner, chunk, comm->get_rank(), ray_size[0], ray_size[1], ray_size[2]);
    std::vector<int> gather_ray_size(comm->get_size() * 3);
    MPI_Gather(&ray_size[0], 3, MPI_INT, gather_ray_size.data(), 3, MPI_INT, owner, MPI_COMM_WORLD);
     
    int proc_size = comm->get_size();
    int proc_rank = comm->get_rank();
    std::vector<int> rcounts(proc_size);
    std::vector<int> displs(proc_size);
    std::vector<char> recv_buffer;
    int all_rays_size = 0;
    if(comm->get_rank() == owner) {
        for(int k = 0; k < comm->get_size(); k ++) {
            printf("owner %d chunk %d id %d raysize %d %d %d \n"
                    , owner, chunk, k, gather_ray_size[k*3]
                    , gather_ray_size[k*3+1], gather_ray_size[k*3+2]);
        }

        printf("owner : ");
        for(int i = 0; i < proc_size; i++) {
            all_rays_size += gather_ray_size[i * 3];
            rcounts[i] = gather_ray_size[i * 3];
            displs[i] = i == 0 ? 0 : (displs[i - 1] + rcounts[i - 1]);
            comm->os<<"gather rays size "<<gather_ray_size[i * 3]<<" rcounts " <<rcounts[i]<<" displs "<<displs[i]<<"\n";
            printf(" count %d  displs %d: ", rcounts[i], displs[i]);
        }
        recv_buffer.resize(all_rays_size);
    }
//    printf("before gather %d %d\n", proc_rank, ray_size[proc_rank * 3] / 4); 
    
    MPI_Gatherv(msg.get_content(), msg.get_content_size(), MPI_CHAR, recv_buffer.data(), rcounts.data(), displs.data(), MPI_CHAR, owner, MPI_COMM_WORLD);
    if(comm->get_rank() == owner)  {
        char* ptr = (char*)recv_buffer.data();
        for(int i = 0; i < proc_size; i++) { 
            printf("before read message chunk %d owner %d %d %d %d\n", chunk, owner
                    , gather_ray_size[3 * i], gather_ray_size[3 * i + 1], gather_ray_size[3 * i + 2]);
            inlist.read_from_stream_message(ptr, gather_ray_size[3 * i + 1], gather_ray_size[3 * i + 2], proc_rank); 
            ptr += gather_ray_size[3 * i];
        }
    }
}

int SyncNode::get_unloaded_chunk(int* list_size) {
    int unloaded_chunk = -1;
    int max_ray_size = 0;
    int chunk_size = ps->get_chunk_size();
    bool simple_mesh = ps->get_simple_trace();
    chunk_size = simple_mesh ? chunk_size - 1 : chunk_size;
 
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
    printf("rank %d before load incoming\n", comm->get_rank()); 
    comm->os<<"rthread load incoming buffer inlist size "<<inlist.size()<<"\n";
    std::unique_lock <std::mutex> lock(inlist.mutex); 
    if(inlist.empty()) { return -1; } 
    if(!inlist.empty() && ps->Exit())
        warn("inlist not empty but prepared to exit\n");

    struct RaysStream *rays_stream;
    if(primary && inlist.primary_size() > 0) {
        rays_stream = inlist.get_primary();
    } else if (!primary && inlist.secondary_size() > 0) {
        rays_stream = inlist.get_secondary();
    } else {
        return 0;
    }
    statistics.start("run => work_thread => load_incoming_buffer-copy");

    ps->set_thread_idle(thread_id, false);
    lock.unlock();
    int copy_size = rays_stream->size;
    comm->os<<"all rays chunk "<<ps->get_current_chunk()<<" size "<<copy_size<<"\n";
    int width = rays_stream->width;
    printf("copy primary size %d\n", copy_size);
 //   comm->os<<"rthread width "<<width <<" logic width "<<rays_stream->logic_width<<"\n";
    memcpy(*rays, rays_stream->get_data(), ps->get_stream_store_capacity() * width * sizeof(float)); 
 
    //printf("%d rays_stream size %d %d %d\n", comm->get_rank(), rays_stream->size, primary, rays_stream->store_width);
    int* ids = (int*)(*rays);
    int psize = std::min(5, copy_size);
    
    if(primary) {
        comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
        for(int i = 0; i < psize; i ++) {
            comm->os<<"#" <<ids[i]<<" "<<ids[i]
                                       <<" "<<ids[i + ps->get_stream_store_capacity() * 9];
         //   if((0xFF & ids[i + ps->get_stream_store_capacity() * 10]) > 4 )
         //       comm->os<<"pri pri chk "<<(0xFF & ids[i + ps->get_stream_store_capacity() * 10])<<" read error\n";
        }
        comm->os<<"\n";
    } else {
        comm->os<<"rthread recv ray render size" <<copy_size<<"\n";
        for(int i = 0; i < psize; i ++) {
            comm->os<<"#" <<ids[i]<<" "<<ids[i]
                                       <<" "<<ids[i + ps->get_stream_store_capacity() * 9]
                                       <<" "<<ids[i + ps->get_stream_store_capacity() * 11];
         //   if((0xFF & ids[i + ps->get_stream_store_capacity() * 11]) > 4 )
         //       comm->os<<"sec pri chk "<<(0xFF & ids[i + ps->get_stream_store_capacity() * 11])<<" read error\n";
        }
        comm->os<<"\n";
    }
    delete rays_stream;

    statistics.end("run => work_thread => load_incoming_buffer-copy");
    return copy_size;
}

bool ray_size_cmp(const std::pair<int, int> a, const std::pair<int, int> b) {
    return a.second > b.second;
}

void SyncNode::schedule(int * list_size) {
    scheduler->chunk_manager->reset_chunk_list();
    //only set local chunks empty, but chunk data keep in memory
    scheduler->chunk_manager->local_chunks.chunks.clear();
    comm->os<<"before schedule chunk "<<ps->get_current_chunk()<<" "<<scheduler->chunk_manager->local_chunks.chunks[0].id<<"\n";
    
    int chunk_size = ps->get_chunk_size();
    int proc_size  = comm->get_size();
    if(proc_size >= chunk_size) 
        return;
    //chunk 按照光线数量进行排序
    std::vector<std::pair<int, int>> chunk_ray_list;
    for(int i = 0; i < chunk_size; i++) 
        chunk_ray_list.emplace_back(std::make_pair(i, list_size[i]));
    //sort
    sort(chunk_ray_list.begin(), chunk_ray_list.end(), ray_size_cmp);

    //找含有最多光线的proc_size个chunk并分配下区
    int max_loaded = chunk_size / proc_size;
    for(int i = 0; i < proc_size; i++) { 
        for(int j = 0; j < chunk_size; j++) {
            if(chunk_ray_list[j].first / max_loaded == i) {
                ps->update_chunk(chunk_ray_list[j].first, i);
                break;
            }
        }
    }
    comm->os<<"new schedule: \n";
    for(int i = 0; i < chunk_size; i++)
        comm->os<<chunk_ray_list[i].first<<" "<<chunk_ray_list[i].second<<" | ";

    scheduler->chunk_manager->update_local_chunks(comm->get_rank());
    scheduler->chunk_manager->local_chunks.init(true/*sync*/);
    comm->os<<"after schedule chunk "<<ps->get_current_chunk()<<" "<<scheduler->chunk_manager->local_chunks.chunks[0].id<<"\n";
}

void SyncNode::synchronize () {
    //MPI_Barrier(MPI_COMM_WORLD);
    statistics.start("run => synchronize => MPI_Allreduce");
    int chunk_size =  ps->get_chunk_size();
    std::vector<int> list_size(chunk_size); 
    
    for(int i = 0; i < chunk_size; i ++){ 
        list_size[i] = outlist[i].size();
    }
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
    
    //如果proc空了 则给一个chunk并进行gather
    //implement get dst chunk, add proc_chunk_map

    int proc_size = comm->get_size();
    statistics.start("run => synchronize => gather_rays");
    for(int chunk = 0; chunk < chunk_size; chunk++) {
        int owner = ps->get_proc(chunk);
        if(owner >= 0) {
            //gether all chunk rays 
            gather_rays(chunk, owner); 
            outlist[chunk].clear();
        }
    }
    statistics.end("run => synchronize => gather_rays");
    printf("rthread proc %d after gather inlist size %d\n", comm->get_rank(), inlist.size());
}

void SyncNode::run(Camera *cam) {
    
    ps->reset();
    scheduler->preprocess(cam, comm->get_size(), true);
    comm->os <<" start run message thread \n";

    int deviceNum = ps->get_dev_num();
    int iter = 0; 
    std::vector<std::thread> workThread;
    while(!ps->Exit()) {
        comm->os<<"run while\n";
        ps->chunk_loaded();
        //开始渲染
        launch_rodent_render(deviceNum, iter==0);
        printf("rank %d after render\n", comm->get_rank()); 
        statistics.start("run => synchronize");
       
        //整理光线 
        RetiredRays::clear_retired_rays(retiredRays, outlist, chunk_size, comm->get_rank());
        //渲染结束，同步调度
        synchronize();
        printf("rank %d after sync\n", comm->get_rank()); 
        printf("after sync\n");
        statistics.end("run => synchronize");
        iter++;
    }

    int chunk_size = ps->get_chunk_size();
    printf("proc %d left rays :", comm->get_rank());
    for(int i = 0; i < chunk_size; i++) 
        printf(" %d || ", outlist[i].size()); 
    printf("\n");

    std::cout<<" chunk size "<<chunk_size<<" ray left "<<outlist[chunk_size - 1].size()<<"\n";
    scheduler->chunk_manager->local_chunks.FIFO_clear();
    return;
}

