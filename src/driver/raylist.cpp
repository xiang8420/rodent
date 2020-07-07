#include "raylist.h"
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
    
    queue.reserve(4 * capacity * store_width);
    queue.resize(capacity * store_width);
    data = queue.data();
    size = 0;
}
void Rays::read_dev_rays(float *rays, int src, int num, int rays_capacity) {
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
    size += num;
}

int Rays::write_dev_rays(float *rays, int buffer_size, size_t rays_capacity) {
    int copy_size = std::min(size, buffer_size); 
    size -= copy_size;
    int src = size; 
    for(int i = 0, k = 0; i < logic_width; i++) {
        if(mask[i]) {
            for(int j = 0; j < copy_size; j++) {
                rays[j + i * rays_capacity] = data[(src + j) * store_width + k];
            }
            k++;
        } else {
        //    for(int j = 0; j < copy_size; j++) {
        //        rays[j + i * rays_capacity] = 0 ; //data[(src + j) * store_width + k];
        //    }
        }
    }
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

void RayList::classification(RayList ** raylist, float *raybuffer,  size_t size, size_t capacity, bool primary) {
    int width       = primary ? 21 : 14;
    int* chunkIds   = (int*)(raybuffer + 9 * capacity);
    printf("classification\n"); 
    for(int i = 0; i < size; i++) {
        int chunk = chunkIds[i] >> 12;
//        printf("chunk %d \n", chunkIds[i]);
        Rays *list = primary ? raylist[chunk]->primary : raylist[chunk]->secondary;
        list->read_dev_rays(raybuffer, i, 1, capacity); 
        if(list->full()){
            printf("error queue is is full %d %d \n", chunk, raylist[chunk]->primary->get_size());   
        }
    }
}
