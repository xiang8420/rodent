#include "raylist.h"
#include <cstring>
#include <fstream>
#include <iostream>

inline int set_ray_mask(bool *mask, int width, bool compact) {
    
    for(int i = 0; i < width; i++ ) {
        mask[i] = true;
    }
    if (compact) {
        if(width == 21) {
            for(int i = 10; i < 15; i++ ) {
                mask[i] = false; 
            }
            return 16; 
        } else {
            return 14;   
        }
    } else {
        return width; 
    } 
}

Rays::Rays(int capacity, int width, bool compact) 
    : capacity(capacity), logic_width(width), compact(compact) {
    mask = new bool[width];
    store_width = set_ray_mask(mask, width, compact);
    
//    queue.reserve(2 * capacity * store_width);
    queue.resize(capacity * store_width);
    size = 0;
}

int Rays::check_capacity(int num) {
    printf("check size %d num %d capacity %d queue size %d queue capacity%ld\n", size, num, capacity, queue.size(), queue.capacity());
    if(num + size >= capacity) {
        printf("put resize %ld capacity%ld %ld\n", queue.size(), queue.capacity(), queue.data());
        capacity += std::max(num, 1048608);
        queue.resize(capacity * store_width);
        printf("end put resize %ld capacity%ld %ld\n", queue.size(), queue.capacity(), queue.data());
    }
}

void Rays::read_from_device(float *rays, int src, int num, int rays_capacity, int rank) {
    std::ofstream os; 
    os.open("out/proc_buffer_" + std::to_string(rank), std::ios::out | std::ios::app ); 
    os<< "check size "<< size <<" num " << num <<" capacity " << capacity <<" queue size "<< queue.size() << " q capacity "<<queue.capacity()<<"\n";
    if(size + num >= capacity) {
        os<<"put resize size + num > capacity"<< size << " "<< num << " "<<capacity<<"\n" ;
    }
    
    check_capacity(num);
    for(int i = 0, k = 0; i < logic_width; i++) {
        if(mask[i]) {
            for(int j = 0; j < num; j++) {
                queue[(size + j) * store_width + k]  = rays[src + j + i * rays_capacity];
            }
            k++;
        }
    }

        os<<"read dev rays \n";
        os<<"rays \n";
        int *ids = (int*)rays;
        for(int i = 0; i < 10; i ++) {
            os<<"||"<<ids[i]<<" "<< ids[i + rays_capacity * 9];
        }
        os<<"\n";
        os<<"list \n";
        ids = (int*)queue.data();
        for(int i = 0; i < 10; i ++) {
            os<<"||"<<ids[i * store_width]<<" "<< ids[i * store_width + 9];
        }
        os<<"\n";
        ids = (int*)queue.data();
        os<<"list 0 " << num<<" "<<size<<"\n";
        for(int i = 0; i < 10; i ++) {
            os<<"#"<<ids[(i + size) * store_width]<<" "<<ids[(i + size) * store_width + 9];
        }
        os<<"\n\n";
    size += num;
}

int Rays::write_to_device_buffer(Rays * buffer, int buffer_size, int rays_capacity, int rank, bool primary) {

    float *rays = buffer->get_data();

    int copy_size = std::min(size, buffer_size); 
    size -= copy_size;
    int src = size; 
    for(int i = 0, k = 0; i < logic_width; i++) {
        if(mask[i]) {
            for(int j = 0; j < copy_size; j++) {
                rays[j + i * rays_capacity] = queue[(src + j) * store_width + k];
            }
            k++;
        } else {
            for(int j = 0; j < copy_size; j++) {
                rays[j + i * rays_capacity] = 0 ; //data[(src + j) * store_width + k];
            }
        }
    }
    buffer->size = copy_size;
    return copy_size;
}

void Rays::read_from_rays(struct Rays* a, int st){
    printf("write to rays ptr\n");
    check_capacity(1);
    if(a->store_width != store_width) printf("width not match");
    for(int i = 0; i < store_width; i++){
        queue[size * store_width + i] = a->queue[st * store_width + i];
    }
    size ++;
}

void Rays::read_from_ptr(char *src_ptr, int copy_size) {
    printf("read frome ptr\n");
    check_capacity(copy_size);
    char* dst_ptr = (char*)get_data() + size * store_width * sizeof(float);
    int copy_length = copy_size * store_width * sizeof(float);
    memcpy(dst_ptr, src_ptr, copy_length);
    size += copy_size; 
    printf("end read frome ptr\n");
}

inline void swap(struct Rays* &a, struct Rays* &b) {
    struct Rays* tmp;
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

void RayList::read_from_device_buffer(RayList ** raylist, float *raybuffer,  size_t size, size_t capacity, bool primary, int rank) {
    int width       = primary ? 21 : 14;
    int* chunkIds   = (int*)(raybuffer + 9 * capacity);
    std::ofstream os; 
    os.open("out/proc_buffer_" + std::to_string(rank), std::ios::out | std::ios::app ); 
    os<<"classification\n"; 
    printf("classification\n");
    int chunk = chunkIds[0] >> 12;
    
    Rays *rays = primary ? raylist[chunk]->primary : raylist[chunk]->secondary;
    rays->read_from_device(raybuffer, 0, size, capacity, rank); 

    os<<"end classification\n"; 
    printf("end classification\n");
//    for(int i = 0; i < size; i++) {
//        int chunk = chunkIds[i] >> 12;
////        printf("chunk %d \n", chunkIds[i]);
//        Rays *list = primary ? raylist[chunk]->primary : raylist[chunk]->secondary;
//        list->read_from_device(raybuffer, i, 1, capacity, rank); 
//        if(list->full()){
//            printf("error queue is is full %d %d \n", chunk, raylist[chunk]->primary->get_size());   
//        }
//    }
}

void RayList::write_to_device_buffer(RayList *buffer, int rank) {

    std::ofstream os; 
    printf("write to device buffer\n");
    os.open("out/proc_buffer_" + std::to_string(rank), std::ios::out | std::ios::app ); 
    os<<"write to device buffer\n"; 
    if(buffer->primary->empty() && !get_primary()->empty()) {
        primary->write_to_device_buffer(buffer->primary, buffer->logic_capacity, buffer->store_capacity, rank, true);
    }
    if(buffer->secondary->empty() && !get_secondary()->empty()) {
        secondary->write_to_device_buffer(buffer->secondary, buffer->logic_capacity, buffer->store_capacity, rank, false);
    }
    os<<"end write to device buffer\n"; 
}

void RayList::read_from_message(char* src_ptr, int msg_primary_size, int msg_secondary_size, int rank) {
    std::ofstream os; 
    os.open("out/proc_buffer_" + std::to_string(rank), std::ios::out | std::ios::app ); 
    printf("read from message\n");
    os<<"read from message\n"; 
    if(msg_primary_size > 0 ) {
        if(primary->size + msg_primary_size > primary->capacity){
            os<<" read from message resize primary size"<<primary->size<<" "<<msg_primary_size<<" "<<primary->capacity<<"\n"; 
        }
        
        primary->read_from_ptr(src_ptr, msg_primary_size); 
        int   primary_length = msg_primary_size * primary->store_width * sizeof(float);
        src_ptr += msg_primary_size * primary->store_width * sizeof(float); 
    } 
    if(msg_secondary_size > 0) {
        if(secondary->size + msg_secondary_size > secondary->capacity){
            os<<" read from message resize secondary size"<<secondary->size<<" "<<msg_secondary_size<<" "<<secondary->capacity<<"\n"; 
        }
        secondary->read_from_ptr(src_ptr, msg_secondary_size); 
//        secondary->check_capacity(msg_secondary_size); 
//        char* secondary_copy_ptr  = (char*)secondary->get_data() + secondary->size * secondary->store_width * sizeof(float);
//        int secondary_length = msg_secondary_size * secondary->store_width * sizeof(float);
//        memcpy(secondary_copy_ptr, src_ptr, secondary_length);
//        secondary->size += msg_secondary_size; 
    }
    os<<"end read from message\n"; 
}

//void RayList::write_to_message(Message *msg){
//
//    if(outList->empty()) return 0;
//  
//    Rays* primary = outList->get_primary(); 
//    int width = primary->store_width; 
//    int* ids = (int*)(primary->get_data());
//    for(int i = 0; i < 5; i ++) {
//        printf(" %d %d$", ids[i*width], ids[i*width + 9]);
//    }
//    printf("\n");
//
//
//    header->primary = outList->get_primary()->size; 
//    header->secondary = outList->get_secondary()->size; 
//
//    printf("serialize size %d %d:", header->primary, header->secondary);
//    
//    int primary_length = header->primary * outList->get_primary()->store_width * sizeof(float);
//    int secondary_length = header->secondary * outList->get_secondary()->store_width * sizeof(float);
//    content.resize(primary_length + secondary_length);
//    printf("Semd Message size %d\n", content.size());  
//
//    char* ptr = content.data();
//    memcpy(ptr, outList->get_primary()->get_data(), primary_length);
//    memcpy(ptr + primary_length, outList->get_secondary()->get_data(), secondary_length);
//    
//    return primary_length + secondary_length;
//}

void RayList::classification(RayList ** raylist, Rays *raybuffer) {
    int *bufferdata   = (int*) raybuffer->get_data();
    int recv_size     = raybuffer->size;
    int width   = raybuffer->store_width;
    if( raybuffer->logic_width == 21) {
        for(int i = 0; i < recv_size; i++) {
            int chunk = bufferdata[i * width + 9]  >> 12;
            raylist[chunk]->primary->read_from_rays(raybuffer, i);
            if(raylist[chunk]->primary->full()){
                printf("error queue is is full %d %d \n", chunk, raylist[chunk]->primary->get_size());   
            }
        }
    } else {
        for(int i = 0; i < recv_size; i++){
            int chunk = bufferdata[i * width + 9]  >> 12;
            raylist[chunk]->secondary->read_from_rays(raybuffer, i);
            if(raylist[chunk]->secondary->full()){
                printf("error queue is is full %d %d \n", chunk, raylist[chunk]->secondary->get_size());   
            }
        }
    }
    raybuffer->clear();
}

