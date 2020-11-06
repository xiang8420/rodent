#ifndef DEF_RAYSTREAMLIST
#define DEF_RAYSTREAMLIST

#include<queue>


struct RaysStream {
    
    float * data;
    //size capacity are ray size,
    int size, width, capacity; 
    
    RaysStream(float *, int, int, int); 

    ~RaysStream(){ delete[] data;}

    float* get_data() {return data; }
};

class RayStreamList {
public:
    std::mutex mutex;

    RayStreamList(){}

    void set_capacity (int, int);
    
    //别忘了delete rays
    RaysStream* get_primary();

    RaysStream* get_secondary(); 
   
     
    void read_from_message(char*, int, int, int); 

    int primary_size(){ return primary.size();}    
    int secondary_size(){ return secondary.size();}    
    bool empty() { return primary.empty() && secondary.empty(); } 
    int size() { return primary.size() + secondary.size(); }
    void lock() { mutex.lock();}
    void unlock() { mutex.unlock();}

    void read_from_ptr(float*, int, int, bool, int); 
    static void read_from_device_buffer(RayStreamList * , float *, size_t , bool , int);
    static void read_from_message(RayList *, char*, int, int);

private:
    std::queue<RaysStream *> primary;
    std::queue<RaysStream *> secondary;
    
    int logic_capacity, store_capacity;  //1048576 1048608 

};

//read from float Array ptr, get device stream rays  capacity = 1048608
RaysStream::RaysStream(float * ptr, int capacity, int width, int copy_size) 
     :capacity(capacity), width(width) 
{
    if(copy_size > capacity)
        error("RaysStream copy size > capacity\n");

    size_t memory = physical_memory_used_by_process();
    //statistics.start("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays => resize");
    data = new float[capacity * width];

   // std::vector<float> tmp(capacity * width);
   // data = tmp;

    //data.resize(capacity * width);
    //statistics.end("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays => resize");
    printf("2new rays resize cap %d %ld kb %ld mb\n", capacity, memory, memory / 1024);
    
    //statistics.start("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays => copy");
    // mask used as ptr mask not this 
    bool *ptr_mask = new bool[width];
    int  ptr_logic_width = width; 
    int  ptr_store_width = set_ray_mask(ptr_mask, ptr_logic_width, RAY_COMPACT);
    
    printf("new RaysStream, ");
    int*a = (int*)data;
    int*b = (int*)ptr;
    for(int i = 0, k = 0; i < ptr_logic_width; i++) {
        if(ptr_mask[i]) {
            for(int j = 0; j < copy_size; j++) {

                data[j + i * capacity] = ptr[j * ptr_store_width + k];
               // if(j < 5  && (i == 0 || i ==9) ) 
               //     printf("| %d %d ", a[j + i * capacity], b[j * ptr_store_width + k]);
            }
            k++;
        } 
    }
    //mask change to ray itself
    size = copy_size; 
    //statistics.end("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays => copy");
}

void RayStreamList::set_capacity (int logic, int store) {
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

void RayStreamList::read_from_message(char* src_ptr, int msg_primary_size, int msg_secondary_size, int rank) {
    std::ofstream os; 
    os.open("out/proc_buffer_" + std::to_string(rank), std::ios::out | std::ios::app ); 
    printf("read from message\n");
    os<<"read from message\n"; 

    int copy_size = std::min(logic_capacity, msg_primary_size);
    os<<"copy size "<<copy_size<<" logic capacity "<<logic_capacity<<" msg primary "<<msg_primary_size<<"\n";
    float *ptr = (float*)src_ptr;
     
    while(copy_size > 0) {
        
        //statistics.start("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays");
        RaysStream * rays   = new struct RaysStream(ptr, store_capacity, 21, copy_size);
        //statistics.end("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays");
        
        int width = RAY_COMPACT ? 16 : 21; 

        int * ids = (int*)ptr;
        os<<"primary copy size "<<copy_size<<"\n"<<" ptr ray";
        for(int i = 0; i < std::min(copy_size, 5); i ++) {
            os<<"| "<< ids[i * width] <<" "<< ids[i * width + 9] << " ";
        } os<<"\n ";

        ids = (int*)(rays->get_data());
        os<<"rays ray ";
        for(int i = 0; i < std::min(copy_size, 5); i ++) {
            os<<"| "<< ids[i] <<" "<< ids[i + 9 * store_capacity] << " ";
        } os<<"\n";


        msg_primary_size -= copy_size; 
        ptr += copy_size * width; 
        primary.push(rays);
        
        copy_size = std::min(logic_capacity, msg_primary_size);
    } 
    
    copy_size = std::min(logic_capacity, msg_secondary_size);
    while(copy_size > 0) {
        //statistics.start("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays secondary");
        RaysStream * rays   = new struct RaysStream(ptr, store_capacity, 14, copy_size);
        //statistics.end("run => message_thread => recv_message => RecvMsg => read_from_message => new Rays secondary");

        int * ids = (int*)ptr;
        os<<"secondary copy size "<<copy_size<<"\n"<<" ptr ray";
        for(int i = 0; i < std::min(copy_size, 5); i ++) {
            os<<"| "<< ids[i * 14] <<" "<< ids[i * 14 + 9] << " ";
        } os<<"\n ";

        ids = (int*)(rays->get_data());
        os<<"rays ray ";
        for(int i = 0; i < std::min(copy_size, 5); i ++) {
            os<<"| "<< ids[i] <<" "<< ids[i + 9 * store_capacity] << " ";
        } os<<"\n";
 

        msg_secondary_size -= copy_size; 
        int width = 14; 
        ptr += copy_size * width; 
        secondary.push(rays);
        
        copy_size = std::min(logic_capacity, msg_secondary_size);
    } 
    os<<"RayStreamList end read from ptr\n"; 
}

void RayStreamList::read_from_ptr(float* rays_ptr, int st, int rays_size, bool isPrimary, int rank) {
    std::ofstream os; 
    os.open("out/proc_buffer_" + std::to_string(rank), std::ios::out | std::ios::app ); 
    printf("read from message\n");
    os<<"read from message\n"; 

    int copy_size = std::min(logic_capacity, rays_size);
    int width = isPrimary ? (RAY_COMPACT ? 16 : 21) : 14; 
    os<<"copy size "<<copy_size<<" logic capacity "<<logic_capacity<<" priamry ? "<<isPrimary<<" msg ray "<<rays_size<<"\n";
    float *ptr = (float*)rays_ptr;
    while(copy_size > 0) {
        RaysStream * rays   = new struct RaysStream(ptr, store_capacity, width, copy_size);

        int * ids = (int*)ptr;
        os<<"rays copy size "<<copy_size<<"\n"<<" ptr ray";
        for(int i = 0; i < std::min(copy_size, 5); i ++) {
            os<<"| "<< ids[i * width] <<" "<< ids[i * width + 9] << " ";
        } os<<"\n ";

        ids = (int*)(rays->get_data());
        os<<"rays ray ";
        for(int i = 0; i < std::min(copy_size, 5); i ++) {
            os<<"| "<< ids[i] <<" "<< ids[i + 9 * store_capacity] << " ";
        } os<<"\n";


        rays_size -= copy_size; 
        ptr += copy_size * width; 
        if(isPrimary)
            primary.push(rays);
        else 
            secondary.push(rays);
        
        copy_size = std::min(logic_capacity, rays_size);
    } 
    os<<"RayStreamList end read from ptr\n"; 
}

void RayStreamList::read_from_device_buffer(RayStreamList * raylist, float *out_buffer,  size_t size, bool primary, int rank) {
//    int width       = primary ? 21 : 14;
//    int* iptr   = (int*)(out_buffer);
//
//    int st = 0, num = 0; 
//    int *chunk_list = (int *)(&iptr[9 * store_capacity]);
//
//    int tmp = chunk_list[0];
//    for(int i = 0; i < size; i++) {
//        int chunk = chunk_list[i] >> 12;
//        if(chunk == tmp) {
//            num ++;
//        } else {
//            //printf("before list read st %d num %d capacity %d dst chunk %d rank %d\n", st, num, capacity, tmp, rank);
//            raylist[tmp].read_from_ptr(out_buffer, st, num, primary, rank);
//            //printf("after list read st %d num %d capacity %d dst chunk %d rank %d ray size %d ray width %d\n", st, num, capacity, tmp, rank, ray.get_size(), ray.get_store_width());
//            tmp = chunk;
//            st = i; 
//            num = 1; 
//        }
//    }
//    raylist[tmp].read_from_ptr(out_buffer, st, num, primary, rank);
}

#endif
