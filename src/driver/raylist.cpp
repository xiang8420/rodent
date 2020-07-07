#include "raylist.h"
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
    data = queue.data();
    size = 0;
}
void Rays::read_dev_rays(float *rays, int src, int num, int rays_capacity, int rank) {
    std::ofstream os; 
    os.open("out/proc_buffer_" + std::to_string(rank), std::ios::out | std::ios::app ); 
    
    if(num + size >= capacity) {
        printf("put resize %ld capacity%ld\n", queue.size(), queue.capacity());
        capacity += 1048608;
        queue.resize(capacity * store_width);
        printf("put resize %ld capacity%ld\n", queue.size(), queue.capacity());
    }
    for(int i = 0, k = 0; i < logic_width; i++) {
        if(mask[i]) {
            for(int j = 0; j < num; j++) {
                data[(size + j) * store_width + k]  = rays[src + j + i * rays_capacity];
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
        ids = (int*)data;
        for(int i = 0; i < 10; i ++) {
            os<<"||"<<ids[i * store_width]<<" "<< ids[i * store_width + 9];
        }
        os<<"\n";
        ids = (int*)data;
        os<<"list 0 " << num<<" "<<size<<"\n";
        for(int i = 0; i < 10; i ++) {
            os<<"#"<<ids[(i + size) * store_width]<<" "<<ids[(i + size) * store_width + 9];
        }
        os<<"\n\n";
    size += num;
}

int Rays::copy_to_buffer(float *rays, int buffer_size, size_t rays_capacity, int rank, bool primary) {

    std::ofstream os; 
    os.open("out/proc_buffer_" + std::to_string(rank), std::ios::out | std::ios::app ); 

    int copy_size = std::min(size, buffer_size); 
    if(primary)
        os<<"copy to primary "<<copy_size<<" size "<<size<<"\n";
    else 
        os<<"copy to secondary "<<copy_size<<" size "<<size<<"\n";

//    size = 0;
    size -= copy_size;
    int src = size; 
    for(int i = 0, k = 0; i < logic_width; i++) {
        if(mask[i]) {
            for(int j = 0; j < copy_size; j++) {
                rays[j + i * rays_capacity] = data[(src + j) * store_width + k];
         //       data[(src + j) * store_width + k] = 0;
            }
            k++;
        } else {
            for(int j = 0; j < copy_size; j++) {
                rays[j + i * rays_capacity] = 0 ; //data[(src + j) * store_width + k];
            }
        }
    }
        
        int* ids = (int*)data;
        for(int i = 0; i < 10; i ++) {
            os<<"#"<<ids[i * store_width]<<" "<< ids[i * store_width + 9];
        }
        os<<"\n";
        ids = (int*)rays;
        for(int i = 0; i < 10; i ++) {
            os<<"#"<<ids[i]<<" "<< ids[i + rays_capacity * 9];
        }
        os<<"\n";
        
        ids = (int*)data;
        os<<"copy to buffer inlist 2 " << copy_size<<" "<<size<<"\n";
        for(int i = 0; i < 10; i ++) {
            os<<"#"<<ids[(i + size) * store_width]<<" "<<ids[(i + size) * store_width + 9];
        }
        os<<"\n\n";
        return copy_size;
}

void Rays::copy_rays(struct Rays* a, int st){
    if(size >= capacity) {
        printf("copy resize %ld capacity%ld\n", queue.size(), queue.capacity());
        capacity += 1048608;
        queue.resize(capacity * store_width);
        printf("copy resize %ld capacity%ld\n", queue.size(), queue.capacity());
    }
    if(a->store_width != store_width) printf("width not match");
    for(int i = 0; i < store_width; i++){
        data[size * store_width + i] = a->data[st * store_width + i];
    }
    size ++;
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

void RayList::copy(RayList *inList, int rank) {
    if(type!="buffer") {
        std::cerr << "only buffer can copy\n";
    }
    if(primary->empty() && !inList->get_primary()->empty()) {
        primary->size = inList->get_primary()->copy_to_buffer(primary->data, logic_capacity, store_capacity, rank, true);
    }
    if(secondary->empty() && !inList->get_secondary()->empty()) {
        secondary->size = inList->get_secondary()->copy_to_buffer(secondary->data, logic_capacity, store_capacity, rank, false);
    }

}

void RayList::classification(RayList ** raylist, Rays *raybuffer) {
    int *bufferdata   = (int*) raybuffer->get_data();
    int recv_size     = raybuffer->size;
    int width   = raybuffer->store_width;
    if( raybuffer->logic_width == 21) {
        for(int i = 0; i < recv_size; i++) {
            int chunk = bufferdata[i * width + 9]  >> 12;
            raylist[chunk]->primary->copy_rays(raybuffer, i);
            if(raylist[chunk]->primary->full()){
                printf("error queue is is full %d %d \n", chunk, raylist[chunk]->primary->get_size());   
            }
        }
    } else {
        for(int i = 0; i < recv_size; i++){
            int chunk = bufferdata[i * width + 9]  >> 12;
            raylist[chunk]->secondary->copy_rays(raybuffer, i);
            if(raylist[chunk]->secondary->full()){
                printf("error queue is is full %d %d \n", chunk, raylist[chunk]->secondary->get_size());   
            }
        }
    }
    raybuffer->clear();
}

void RayList::classification(RayList ** raylist, float *raybuffer,  size_t size, size_t capacity, bool primary, int rank) {
    int width       = primary ? 21 : 14;
    int* chunkIds   = (int*)(raybuffer + 9 * capacity);
    printf("classification\n"); 
    int chunk = chunkIds[0] >> 12;
    
    Rays *list = primary ? raylist[chunk]->primary : raylist[chunk]->secondary;
    list->read_dev_rays(raybuffer, 0, size, capacity, rank); 

//    for(int i = 0; i < size; i++) {
//        int chunk = chunkIds[i] >> 12;
////        printf("chunk %d \n", chunkIds[i]);
//        Rays *list = primary ? raylist[chunk]->primary : raylist[chunk]->secondary;
//        list->read_dev_rays(raybuffer, i, 1, capacity, rank); 
//        if(list->full()){
//            printf("error queue is is full %d %d \n", chunk, raylist[chunk]->primary->get_size());   
//        }
//    }
}
