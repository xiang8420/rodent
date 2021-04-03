
class Node {

protected:
    Communicator *comm;
    ProcStatus *ps;
    Scheduler *scheduler;

    std::mutex  thread_mutex;
    
    std::condition_variable inList_not_empty; //
    std::condition_variable render_start; //

public: 
    int max_used_mem, pre_mem;
    
    Node(Communicator *comm, ProcStatus *ps, Scheduler *scheduler);

    ~Node();

    ProcStatus * get_proc_status(){return ps;}

    Communicator * get_communicator(){ return comm; }

    void launch_rodent_render(int, bool); 

    virtual void run(Camera *) = 0;

    virtual int load_incoming_buffer(float **rays, bool primary, int thread_id) = 0;
    
    virtual void save_outgoing_buffer(float *, size_t, bool) = 0; 

    void loop_check(float i); 
};

Node::Node(Communicator *comm, ProcStatus *ps, Scheduler *scheduler) 
    : comm(comm), ps(ps), scheduler(scheduler)
{
    pre_mem = 0;
    max_used_mem = 0;
}

Node::~Node() {}

void Node::launch_rodent_render(int devNum, bool generate_rays) {

    Camera *camera = scheduler->camera;
    int sppProc = scheduler->get_spp(); 
    int sppDev = sppProc / devNum;

    statistics.start("run => work_thread => render => run");
    // only used for async dynamic schedule
    ps->chunk_loaded(); 

    std::vector<std::thread> workThread;
    for(int i = 0; i < devNum; i++) {
        int* region = scheduler->get_render_block();
        Settings settings {
            Vec3 { camera->eye.x, camera->eye.y, camera->eye.z },
            Vec3 { camera->dir.x, camera->dir.y, camera->dir.z },
            Vec3 { camera->up.x, camera->up.y, camera->up.z },
            Vec3 { camera->right.x, camera->right.y, camera->right.z },
            camera->w, camera->h,
            Vec4_i32 { region[0], region[1], region[2], region[3]},
            sppDev,
            ps->get_simple_trace()
        };
        int rnd = camera->iter * (comm->get_rank() + 1) * devNum + i;

        //printf("region %d %d %d %d\n", region[0], region[1], region[2], region[3]);
        //printf("launch render chunk %d %d\n", ps->get_current_chunk(), ps->get_next_chunk());
        workThread.emplace_back(render, &settings, rnd, i, ps->get_current_chunk(), ps->get_next_chunk(), generate_rays);
    }
    render_start.notify_all();

    for(auto &thread: workThread)
        thread.join();
    
    statistics.end("run => work_thread => render => run");
    comm->os<<comm->get_rank()<<" before set render start"<<"\n";
}


void Node::loop_check(float i) {
    if(1) {    
//        int mem = physical_memory_used_by_process();
//        if(/*i> 1000 ||*/mem != pre_mem) {
// //           comm->os<<i<<" memory use "<<mem<<"\n"; 
//            max_used_mem = std::max(mem, max_used_mem);
//            pre_mem = mem;
//        }
// //       printf("%d mark %f\n", comm->get_rank(), i);
// //       comm->os<<comm->get_rank()<<" mark "<< i <<"\n";
    }
}

void swap_array(float *rays, int k, int j, int width) {
    int *tmp = new int[width];
    memcpy(tmp, &rays[width * k], width * sizeof(float));
    memcpy(&rays[width * k], &rays[width * j], width * sizeof(float));
    memcpy(&rays[width * j], tmp, width * sizeof(float));
    delete[] tmp; 
}

void sort_ray_array(float* rays, int size, int chunk_size, bool primary){
    int width = primary ? PRIMARY_WIDTH : SECONDARY_WIDTH;
    std::vector<int> end;
    std::vector<int> begin;

    for(int i = 0; i < chunk_size; i++) { 
        begin.emplace_back(0);
        end.emplace_back(0);
    }
    int * iptr = (int*) rays;
    for(int i = 0; i < size; i++) { 
        int chunk = iptr[i * width + 9];
        end[chunk]++;
    }
    
    int n = 0;
    for(int i = 0; i < chunk_size; i++) { 
        begin[i] = n;
        n += end[i];
        end[i] = n;
    }

    for(int i = 0; i < chunk_size; i++) {
        int st = begin[i];
        int ed = end[i];
        int j = st;
        while ( j < ed) {
            int cid = iptr[j * width + 9]; 
            if (cid != i) {
                int k = begin[cid]++;
                swap_array(rays, k, j, width);                
            } else {
                j++;
            }
        }
    }
}

struct RetiredRays {
    float* data;
    int size;
    bool primary;
    int width;

    RetiredRays(float* rays, int size, bool primary)
        :size(size), primary(primary) 
    {
        int width = primary ? PRIMARY_WIDTH : SECONDARY_WIDTH;
        int capacity = size * width;
        data = new float[capacity]; 
        //printf("new retired ray size %d width %d \n", size, width);
        memcpy((char*)data, (char*)rays, capacity * sizeof(float));
    }
    
    RetiredRays(int size, bool primary)
        :size(size), primary(primary) 
    {
        int width = primary ? PRIMARY_WIDTH : SECONDARY_WIDTH;
        int capacity = size * width;
        data = new float[capacity]; 
    }

    ~RetiredRays() {
        delete[] data;
    }
    
    static void compact_retired_rays(std::vector<RetiredRays*> &retiredRays, int chunk_size) {
        int primary_size = 0;    
        int secondary_size = 0;    
        for(int i = 0; i < retiredRays.size(); i++) {
            if(retiredRays[i]->primary) 
                primary_size += retiredRays[i]->size;
            else
                secondary_size += retiredRays[i]->size;
        }
        RetiredRays* retired_primary = new RetiredRays(primary_size, true);
        RetiredRays* retired_secondary = new RetiredRays(secondary_size, false);
        
        float *primary_ptr = retired_primary->data;
        float *secondary_ptr = retired_secondary->data;
        for(int i = 0; i < retiredRays.size(); i++) {
            if(retiredRays[i]->primary) { 
                int copy_size = retiredRays[i]->size * PRIMARY_WIDTH;
                memcpy(primary_ptr, retiredRays[i]->data, copy_size * sizeof(float));
                primary_ptr += copy_size;
            } else {
                int copy_size = retiredRays[i]->size * SECONDARY_WIDTH;
                memcpy(secondary_ptr, retiredRays[i]->data, copy_size * sizeof(float));
                secondary_ptr += copy_size;
            }
        }

        while(!retiredRays.empty()) {
            RetiredRays* rays = retiredRays.back();
            retiredRays.pop_back();
            delete rays;
        }

        sort_ray_array(retired_primary->data, primary_size, chunk_size, true);
        sort_ray_array(retired_secondary->data, secondary_size, chunk_size, false);

        retiredRays.emplace_back(retired_primary);
        retiredRays.emplace_back(retired_secondary);

        //clear retired rays
        //copy priamry and secondary to retired ryas

    }
    
    static void clear_retired_rays(std::vector<RetiredRays*> &retiredRays, RayStreamList * outlist, int chunk_size, int proc_rank) {
        //printf("clear retired rays %d \n", retiredRays.size());
       // for(auto& rays: retiredRays) 
       //     printf("%d \n", rays->size);
        
        statistics.start("run => message_thread => clear_retired => get mutex");
        compact_retired_rays(retiredRays, chunk_size);       
        std::unique_lock <std::mutex> lock(outlist[0].mutex); 
        while(!retiredRays.empty()) {
            RetiredRays* rays = retiredRays.back();
            retiredRays.pop_back();
            statistics.start("run => message_thread => clear_retired => read_from_device");
            RayStreamList::read_from_device_buffer(outlist, rays->data, rays->size, rays->primary, chunk_size, proc_rank);
            statistics.end("run => message_thread => clear_retired => read_from_device");
            delete rays;
        }
       // for(int i = 0; i < chunk_size; i ++) {
       //     if(outlist[i].ray_size() > 0)
       //         comm->os<<"clear retired rays "<<i<<" "<<outlist[i].ray_size()<<"\n";
       // }
    }

};
