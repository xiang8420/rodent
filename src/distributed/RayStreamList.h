#ifndef DEF_RAYSTREAMLIST
#define DEF_RAYSTREAMLIST

#include<queue>

#define MAX_STREAM_LIST 32 

struct RaysStream {
    
    float * data;
    //size capacity are ray size,
    int size, width;
    int logic_capacity, store_capacity; 
    
    RaysStream(int lcap, int scap, int width)
        :logic_capacity(lcap), store_capacity(scap), width(width) 
    {
        printf("construct empty RaysStream\n");
        data = new float[store_capacity * width];
        size = 0; 
    }
    //read from float Array ptr, get device stream rays  capacity = 1048608
    RaysStream(char * ptr, int lcap, int scap, int width, int copy_size, bool stream) 
        :logic_capacity(lcap), store_capacity(scap), width(width) 
    {
        if(copy_size > logic_capacity)
            error("RaysStream copy size  > logic_capacity \n", copy_size, " ", logic_capacity);

        data = new float[store_capacity * width];
        if(stream) {
            printf("new raysstream :");
            int* iptr = (int*)ptr;
            if(copy_size > 7)
                //for(int i = copy_size-6; i < copy_size; i++) {
                for(int i = 0; i < 5; i++) {
                    printf("%d %d &&", iptr[i], iptr[i + 9 * store_capacity]);
                }
                printf("\n");
            memcpy((char*)data, ptr, store_capacity * width * sizeof(float));
            iptr = (int*)data;
            if(copy_size > 7)
                //for(int i = copy_size-6; i < copy_size; i++) {
                for(int i = 0; i < 5; i++) {
                    printf("%d %d &&", iptr[i], iptr[i + 9 * store_capacity]);
                }
        } else {
            // mask used as ptr mask not this 
            int  ptr_logic_width = width; 
            float* fptr = (float*)ptr;
            for(int i = 0; i < ptr_logic_width; i++) {
                for(int j = 0; j < copy_size; j++) 
                    data[j + i * store_capacity] = fptr[j * width + i];
            }
        }
        //mask change to ray itself
        size = copy_size; 
        //statistics.end("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays => copy");
    }

    void fill(float * ptr, int st, int copy_size, int rank) {
        int* iptr = (int*)ptr;
        int ed = size + copy_size;
        //printf("fill size %d copy_size %d logic_capacity%d\n", size, copy_size, logic_capacity);
        assert(ed <= logic_capacity);
        for(int i = 0; i < width; i++) {
            for(int j = 0; j < copy_size; j++) {
                data[i * store_capacity + j + size] = ptr[(j + st) * width + i];
            }
        }
        size = ed; 
    }

    ~RaysStream(){ delete[] data;}
    bool full() {return size == logic_capacity;}
    float* get_data() {return data; }
    int get_size() {return size; }
};

class RayStreamList {
public:
    std::mutex mutex;
    std::condition_variable cond_empty; 
    std::condition_variable cond_full; 

    RayStreamList(){}
    
    RayStreamList (int, int); 
    
    void set_capacity (int, int);
    
    //别忘了delete rays
    RaysStream* get_primary();

    RaysStream* get_secondary(); 
     
    void read_from_array_message(char*, int, int, int); 
    void read_from_stream_message(char*, int, int, int); 
    int primary_size(){ return primary.size();}    
    int secondary_size(){ return secondary.size();}    
    
    int get_head_primary_size() { return primary.size() > 0 ? primary.front()->get_size() : 0; }
    int get_head_secondary_size() { return secondary.size() > 0 ? secondary.front()->get_size() : 0; }
    int size() { return primary.size() + secondary.size(); }
    int get_logic_capacity() {return logic_capacity; }
    int get_store_capacity() {return store_capacity; }
    
    bool empty() { return primary.empty() && secondary.empty(); } 
    bool full() { return primary.size() + secondary.size() > MAX_STREAM_LIST; }
    void lock() { mutex.lock(); }
    void unlock() { mutex.unlock(); }
    void full_notify() { cond_full.notify_all(); }
    void empty_notify() { cond_empty.notify_all(); }

    void clear();
    void read_from_ptr(float*, int, int, bool, int); 
    RaysStream* get_free_stream(bool);
    static void read_from_device_buffer(RayStreamList * , float *, size_t , bool, int, int);
    static void swap(RayStreamList &, RayStreamList &); 


private:
    std::queue<RaysStream *> primary;
    std::queue<RaysStream *> secondary;
    
    int logic_capacity, store_capacity;  //1048576 1048608 
};

void RayStreamList::set_capacity (int logic, int store) {
    logic_capacity = logic;
    store_capacity = store;
}

RayStreamList::RayStreamList (int logic, int store) {
    logic_capacity = logic;
    store_capacity = store;
}

RaysStream* RayStreamList::get_primary() {
    if(!primary.empty()) {
        RaysStream * rays = primary.front();
        primary.pop();
        return rays;
    } else {
        printf("return get primary null \n");
        return NULL; 
    }
}

RaysStream* RayStreamList::get_secondary() {
    if(!secondary.empty()) {
        RaysStream * rays = secondary.front();
        secondary.pop();
        return rays;
    } else {
        printf("return get secondary null \n");
        return NULL;
    }
}

// msg_primary_size is rays num 
void RayStreamList::read_from_array_message(char* src_ptr, int msg_primary_size, int msg_secondary_size, int rank) {
    std::ofstream os; 
    os.open("out/proc_buffer_" + std::to_string(rank), std::ios::out | std::ios::app ); 
    printf("read from message\n");
    os<<"read from message\n"; 

    int copy_size = std::min(logic_capacity, msg_primary_size);
    os<<"copy size "<<copy_size<<" logic capacity "<<logic_capacity<<" msg primary "<<msg_primary_size<<"\n";
    float *ptr = (float*)src_ptr;
     
    while(copy_size > 0) {
        
        //statistics.start("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays");
        RaysStream * rays   = new struct RaysStream((char*)ptr, logic_capacity, store_capacity, PRIMARY_WIDTH, copy_size, false);
        //statistics.end("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays");
        
        int width = PRIMARY_WIDTH; 




        msg_primary_size -= copy_size; 
        ptr += copy_size * width; 
        primary.push(rays);
        
        copy_size = std::min(logic_capacity, msg_primary_size);
    } 
    
    copy_size = std::min(logic_capacity, msg_secondary_size);
    while(copy_size > 0) {
        //statistics.start("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays secondary");
        RaysStream * rays   = new struct RaysStream((char*)ptr, logic_capacity, store_capacity, SECONDARY_WIDTH, copy_size, false);
        //statistics.end("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays secondary");

        int * ids = (int*)ptr;
        os<<"secondary copy size "<<copy_size<<"\n"<<" ptr ray";
        for(int i = 0; i < std::min(copy_size, 5); i ++) {
            os<<"| "<< ids[i * SECONDARY_WIDTH] <<" "<< ids[i * SECONDARY_WIDTH + 9] << " ";
        } os<<"\n ";

        ids = (int*)(rays->get_data());
        os<<"rays ray ";
        for(int i = 0; i < std::min(copy_size, 5); i ++) {
            os<<"| "<< ids[i] <<" "<< ids[i + 9 * store_capacity] << " ";
        } os<<"\n";
 

        msg_secondary_size -= copy_size; 
        int width = SECONDARY_WIDTH; 
        ptr += copy_size * width; 
        secondary.push(rays);
        
        copy_size = std::min(logic_capacity, msg_secondary_size);
    } 
    os<<"RayStreamList end read from ptr\n"; 
}
// msg_primary_size is stream num
void RayStreamList::read_from_stream_message(char* src_ptr, int msg_primary_size, int msg_secondary_size, int rank) {
    printf("m %d  recv RayMsg primary %d secondary %d\n", rank, msg_primary_size, msg_secondary_size) ;
    for(int i = 0; i < msg_primary_size; i++) {
        int num = ((int*)src_ptr)[0];
        printf("read from stream msg primary num %d\n", num);
        src_ptr += sizeof(int);
        if(num > 0) {
            RaysStream * rays   = new struct RaysStream(src_ptr, logic_capacity, store_capacity, PRIMARY_WIDTH, num, true);
            printf("recv get new ray  %d : ", num);
           // int* iptr = (int*)rays->data;
           // if(num > 7)
           //     //for(int i = num-6; i < num; i++) 
           //     for(int i = 0; i < 5; i++) 
           //         printf("%d %d msg ", iptr[i], iptr[i + 9 * store_capacity]);
           // printf("\n");
            src_ptr += store_capacity * PRIMARY_WIDTH * sizeof(float); 
            primary.push(rays);
        }
    }
    for(int i = 0; i < msg_secondary_size; i++) {
        int num = ((int*)src_ptr)[0];
        printf("read from stream msg secondary num %d\n", num);
        src_ptr += sizeof(int);
        if(num > 0) {
            RaysStream * rays   = new struct RaysStream(src_ptr, logic_capacity, store_capacity, SECONDARY_WIDTH, num, true);
            int* iptr = (int*)rays->data;
          //  if(num > 7)
          //      //for(int i = num-6; i < num; i++) 
          //      for(int i = 0; i < 5; i++) 
          //          printf("%d %d msg ", iptr[i], iptr[i + 9 * store_capacity]);
          //  printf("\n");
            src_ptr += store_capacity * SECONDARY_WIDTH * sizeof(float); 
            secondary.push(rays);
        }
    }
}
void RayStreamList::read_from_ptr(float* rays_ptr, int st, int rays_size, bool isPrimary, int rank) {
//    std::ofstream os; 
//    os.open("out/proc_buffer_" + std::to_string(rank), std::ios::out | std::ios::app ); 
//    printf("read from message\n");
//    os<<"read from message\n"; 

    int copy_size = std::min(logic_capacity, rays_size);
    int width = isPrimary ? PRIMARY_WIDTH : SECONDARY_WIDTH; 
//    os<<"copy size "<<copy_size<<" logic capacity "<<logic_capacity<<" priamry ? "<<isPrimary<<" msg ray "<<rays_size<<"\n";
    float *ptr = (float*)rays_ptr;
    while(copy_size > 0) {
        RaysStream * rays   = new struct RaysStream((char*)ptr, logic_capacity, store_capacity, width, copy_size, false);

        int * ids = (int*)ptr;
//        os<<"rays copy size "<<copy_size<<"\n"<<" ptr ray";
//        for(int i = 0; i < std::min(copy_size, 5); i ++) {
//            os<<"| "<< ids[i * width] <<" "<< ids[i * width + 9] << " ";
//        } os<<"\n ";

//        ids = (int*)(rays->get_data());
//        os<<"rays ray ";
//        for(int i = 0; i < std::min(copy_size, 5); i ++) {
//            os<<"| "<< ids[i] <<" "<< ids[i + 9 * store_capacity] << " ";
//        } os<<"\n";


        rays_size -= copy_size; 
        ptr += copy_size * width; 
        if(isPrimary)
            primary.push(rays);
        else 
            secondary.push(rays);
        
        copy_size = std::min(logic_capacity, rays_size);
    } 
//    os<<"RayStreamList end read from ptr\n"; 
}

void RayStreamList::clear() {
    while (!primary.empty()) { 
        RaysStream * rays = primary.front();
        primary.pop();
        delete rays;
    }
    while (!secondary.empty()) { 
        RaysStream * rays = secondary.front();
        secondary.pop();
        delete rays;
    }
}

RaysStream* RayStreamList::get_free_stream(bool Primary) {
    std::queue<RaysStream *> &queue = Primary ? primary : secondary;
    if(!queue.empty() && !queue.back()->full()) {
        return queue.back();
    } else {
        RaysStream* rays = new RaysStream(logic_capacity, store_capacity, Primary ? PRIMARY_WIDTH : SECONDARY_WIDTH);
        queue.push(rays);
        return rays;
    }
}

void RayStreamList::read_from_device_buffer(RayStreamList * outList, float *out_buffer,  size_t buffer_size, bool primary, int chunk_size, int rank) {
//    auto os = std::ofstream("out/proc_read_from_buffer_" + std::to_string(rank));

    RaysStream **writeList = new RaysStream*[chunk_size];
    for(int i = 0; i < chunk_size; i ++) 
        writeList[i] = NULL; //outList[i].get_free_stream(primary);
    
    printf("read from divice buffer %d \n", buffer_size); 
    
    int width   = primary ? PRIMARY_WIDTH : SECONDARY_WIDTH;
    int* iptr   = (int*)(out_buffer);
    int st = 0, num = 0; 
    int tmp = iptr[9];
    for(int i = 0; i < buffer_size; i++) {
        int chunk = iptr[i * width + 9];
        if(chunk == tmp && i != buffer_size - 1) {
            num ++;
        } else {
//            os<<"before list read st "<< st <<" num "<< num <<" tmp "<< tmp<<"\n";
            //printf("before fill left %d num %d\n", left, num);
            while(num > 0) {
                outList[tmp].mutex.lock(); 
                if(writeList[tmp] == NULL || (writeList[tmp]->full() && num > 0))
                    writeList[tmp] = outList[tmp].get_free_stream(primary);
                int left = std::min(num, writeList[tmp]->logic_capacity - writeList[tmp]->size);
                writeList[tmp]->fill(out_buffer, st, left, rank);
                outList[tmp].mutex.unlock(); 
                
                st += left;
                num -= left;
            }
            //printf("after list read st %d num %d capacity %d dst chunk %d rank %d ray size %d ray width %d\n", st, num, capacity, tmp, rank, ray.get_size(), ray.get_store_width());
            tmp = chunk;
            st = i; 
            num = 1; 
        }
    }
//    printf("after save rays outlist size chunk size%d: ", chunk_size);
//    for(int i = 0; i < chunk_size; i ++) {
//        printf("chunk %d:\n", i);
//        int rays_size = outList[i].primary_size();
////        printf("ray_size = %d", rays_size);
//        RaysStream * rays; 
//        for(int k = 0; k < rays_size; k++) {
//            rays = outList[i].get_primary();
//            if(rays->get_size() > 0)
//                outList[i].primary.push(rays);
//            else 
//                delete rays;
//        }
//        
//        rays_size = outList[i].secondary_size();
////        printf("ray_size = %d\n", rays_size);
//        for(int k = 0; k < rays_size; k++) {
//            rays = outList[i].get_secondary();
//            if(rays->get_size() > 0)
//                outList[i].secondary.push(rays);
//            else 
//                delete rays;
//        }
//    }
}

void RayStreamList::swap(RayStreamList &list1, RayStreamList &list2) {
    std::swap(list1.primary,   list2.primary);
    std::swap(list1.secondary, list2.secondary); 
}


#endif
