#ifndef DEF_RAYARRAYLIST
#define DEF_RAYARRAYLIST

#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <mutex>

size_t physical_memory_used_by_process(); 

//array of struct ray, used for storage primary or secondary depends on width 
class RaysArray {
public:
    float* data;

    RaysArray(){}; 

    void clear(); 
    
    ~RaysArray(){clear();}

    void resize(int capacity, int width); 

    void read_device_buffer(float *rays, int src, int num);
    
    void read_from_ptr(char *src_ptr, int copy_size);
    
    int check_capacity(int);

    float& operator[](int id) { return data[id];}
    bool full() { return size == capacity; }
    bool empty() {return size <= 0;}
    float* get_data() {return data; }
    int get_size() { return size;}
    void set_size(int n) { size = n;}
    int get_capacity() {return capacity;}
    int get_store_width(){return store_width;}
    int get_logic_width(){return logic_width;}
private:
    
    //size capacity are ray size,
    int size, capacity; 
    int logic_width, store_width;

};

class RayArrayList {

public:
    std::mutex mutex;

    RayArrayList(const RayArrayList&, int); 

    RayArrayList();

    void set_capacity(); 
    
    void clear();

    //void read_from_message(char* , int, int);
    
    static void read_from_device_buffer(RayArrayList * , float *, size_t , bool , int);
    
    static void read_from_message(RayArrayList *, char*, int, int);

    RaysArray& get_primary() { return primary;}
    RaysArray& get_secondary() { return secondary;}
    int primary_size(){ return primary.get_size();}    
    int secondary_size(){ return secondary.get_size();}    
    int primary_store_width() { return primary.get_store_width(); }
    int secondary_store_width() { return secondary.get_store_width(); }
    int size() { return primary.get_size() + secondary.get_size(); }
    bool empty() { return primary.empty() && secondary.empty(); } 
    void lock() {mutex.lock();}
    void unlock() {mutex.unlock();}

private:
    struct RaysArray primary;
    struct RaysArray secondary;
    int rank;
};


void RaysArray::resize(int cap, int width) 
{
    capacity = cap;
    
    logic_width = width;
    store_width = width; 
    
    data = new float[capacity * store_width];

    size_t memory = physical_memory_used_by_process();
    printf("1new rays resize cap %d %ld kb %ld mb\n", capacity, memory, memory / 1024);

    size = 0;
}

void RaysArray::clear()
{
    printf("clear rayArray\n");
    //delete[] data;
    size = 0;
}

int RaysArray::check_capacity(int num) {
    printf("check size %d num %d capacity %d \n", size, num, capacity);
    //statistics.start("run => message_thread => RayArrayList => check_capacity");
    if(num + size > capacity) {
        size_t memory = physical_memory_used_by_process();
        //printf("rays resize %ld kb %ld mb\n", memory, memory / 1024);
        capacity += std::max(num, 1024);
        
        //statistics.start("run => message_thread => check_capacity");
        
        float *new_data = new float[capacity * store_width]; 
        if(size > 0) {
            printf("before memcpy check size %d num %d  width %d capacity %d \n", size, num, store_width, capacity);
            memcpy((char*)new_data, (char*)data, size * store_width * sizeof(float));
            printf("after memcpy check size %d num %d  width %d capacity %d \n", size, num, store_width, capacity);
       
            printf("check capacity: \n"); 
            int* ids = (int*)(data);
            for(int i = 0; i < 5; i ++) {
                printf(" %d %d$", ids[i*store_width], ids[i*store_width + 9]);
            }
            printf("\n");
            ids = (int*)(new_data);
            for(int i = 0; i < 5; i ++) {
                printf(" %d %d$", ids[i*store_width], ids[i*store_width + 9]);
            }
            printf("\n");
        }
        delete[] data;
        data = new_data;

        //statistics.end("run => message_thread => check_capacity");
        //printf("store width %d capacity %d\n", capacity, store_width);
        //printf("after rays resize %ld kb %ld mb\n", memory, memory / 1024);
    }
    printf("after check size %d num %d capacity %d \n", size, num, capacity);
    //statistics.end("run => message_thread => RayArrayList => check_capacity");
    return capacity;
}

void RaysArray::read_device_buffer(float *rays, int src, int num) {
    int* iptr = (int *) data;
    for(int i = 0; i < num; i++) {
        int k = 0;
        for(int j = 0; j < logic_width; j++) {
        //    if(j < 5)
        //        printf("| %d %d %d", (size + i) * store_width + k, (src + i) * logic_width + j, data.size(), num);
            data[(size + i) * store_width + k]  = rays[(src + i) * logic_width + j];
            k++; 
        }
    }
    size += num;
}

void RaysArray::read_from_ptr(char *src_ptr, int copy_size) {
    if(copy_size == 0 || check_capacity(copy_size)==0) return;

    char* dst_ptr = (char*)get_data() + size * store_width * sizeof(float);
    int copy_length = copy_size * store_width * sizeof(float);
    printf("RaysArray read from copy size %d length %d \n",  copy_size, copy_length);
    memcpy(dst_ptr, src_ptr, copy_length);
    size += copy_size; 
}

inline void swap(struct RaysArray* &a, struct RaysArray* &b) {
    struct RaysArray* tmp;
    tmp = a;
    a = b;
    b = tmp;
}

inline void swap(struct Ray* &a, struct Ray* &b) {
    struct Ray* tmp;
    tmp = a;
    a = b;
    b = tmp;
}
     
inline void swap(float **a, float **b) {
    float** tmp;
    *tmp = *a;
    *a = *b;
    *b = *tmp;
}

//存储问题 可能从这里解决
void RayArrayList::clear() {
    //printf("\nraylist clear\n\n");
    primary.clear() ; 
    secondary.clear(); 
}

RayArrayList::RayArrayList(const RayArrayList& a, int r) {
    //printf("cant copy construct\n");
    primary.resize(0, PRIMARY_WIDTH);
    secondary.resize(0, SECONDARY_WIDTH);
    rank = r;
}

RayArrayList::RayArrayList() {
    primary.resize(0, PRIMARY_WIDTH);
    secondary.resize(0, SECONDARY_WIDTH);
}

void RayArrayList::set_capacity() {
    primary.resize(0, PRIMARY_WIDTH);
    secondary.resize(0, SECONDARY_WIDTH);
}

void RayArrayList::read_from_device_buffer(RayArrayList * raylist, float *out_buffer,  size_t size, bool primary, int chunk_size) {
    int width       = primary ? PRIMARY_WIDTH : SECONDARY_WIDTH;
    int* chunkIds   = (int*)(out_buffer);

    std::lock_guard <std::mutex> lock(raylist[0].mutex); 
    int st = 0, num = 0; 
    int tmp = chunkIds[9];
    printf("read from divice buffer %ld\n", size);
    for(int i = 0; i < chunk_size; i++) { 
        if(primary && raylist[i].primary.check_capacity(size) == 0) return;
        if(!primary && raylist[i].secondary.check_capacity(size) == 0) return;
    }


    for(int i = 0; i < size; i++) {
        int chunk = chunkIds[i * width + 9];
        if(chunk == tmp) {
            num ++;
        } else {
            RaysArray &ray = primary ? raylist[tmp].primary : raylist[tmp].secondary;
            //printf("before list read st %d num %d capacity %d dst chunk %d rank %d\n", st, num, capacity, tmp, rank);
            ray.read_device_buffer(out_buffer, st, num);
            //printf("after list read st %d num %d capacity %d dst chunk %d rank %d ray size %d ray width %d\n", st, num, capacity, tmp, rank, ray.get_size(), ray.get_store_width());
            tmp = chunk;
            st = i; 
            num = 1; 
        }
    }
    RaysArray &ray = primary ? raylist[tmp].primary : raylist[tmp].secondary;
    ray.read_device_buffer(out_buffer, st, num); 
}

void RayArrayList::read_from_message(RayArrayList* rayList, char* ptr, int msg_primary_size, int msg_secondary_size) {
    std::ofstream os; 
    os.open("out/proc_buffer_master", std::ios::out | std::ios::app ); 
    os<<"master read from message: all primary size "<<msg_primary_size<<" secondary "<<msg_secondary_size<<"\n"; 
    if(msg_primary_size > 0) {
        int* id_ptr = (int*) ptr;
        int width = PRIMARY_WIDTH;
        int st = 0, ed = 0; 
        int tmp = id_ptr[9];
        while(ed < msg_primary_size - 1) {
            ed ++;
            int chunk = id_ptr[ed * width + 9];
            if( chunk != tmp || ed == msg_primary_size - 1) {
                os<<"copy st "<<st<<" ed "<<ed<<" ed*width+9 "<<id_ptr[ed * width + 9]<<" chunk "<<chunk<<" tmp "<<tmp<<"\n"; 
                rayList[tmp].get_primary().read_from_ptr((char*)ptr + st * width * sizeof(float), ed - st);
                tmp = chunk;
                st = ed; 
            }
        } 
        ptr = ptr + msg_primary_size * width * sizeof(float); 
    }
    os<<"secondary \n"; 
    if(msg_secondary_size > 0) {
        int width = SECONDARY_WIDTH;// rayList[0].get_secondary().get_store_width();
        int* id_ptr = (int*) ptr;
        int st = 0, ed = 0; 
        int tmp = id_ptr[9];
        os<<"msg_secondary_size "<< msg_secondary_size<<" :";
        while(ed < msg_secondary_size) {
            ed ++;
            int chunk = id_ptr[ed * width + 9];
            os <<id_ptr[ed * width + 9]<<" "<<chunk<<" ";
            if( chunk != tmp || ed == msg_secondary_size - 1) {
         //       os<<"copy st "<<st<<" ed "<<ed<<"chunk "<<chunk<<" tmp "<<tmp<<"\n"; 
                rayList[tmp].get_secondary().read_from_ptr((char*)ptr + st * width * sizeof(float), ed - st);
                tmp = chunk;
                st = ed; 
            } 
        } 
        os<<"\n";
    }
    os<<"end read from message\n\n"; 
}

#endif
