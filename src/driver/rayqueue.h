
struct Ray {
    float *data;
    int size, capacity, width;
  
    Ray(int capacity, int width): capacity(capacity), width(width) {
        data = new float[capacity * width];
        size = 0;
    }
 
    ~Ray(){
        delete[] data;
    }

    float& operator[](int id) const {
       return data[id]; 
    }
   
    bool isfull() { return size == capacity; }
   
    bool isempty() {
        if(size < 0){printf("ray size < 0 %d\n", size);}
        return size == 0;
    }
   
    float* rays() {return data; }
   
    int* ids() {return (int*)data;}
   
    float* datas() {return &data[capacity];}
   
    void clear() {size = 0;}
   
    int get_size() { return size;}
   
    int get_capacity() {return capacity;}
 
    // queue store is different from orign ray
    void put(float *queue, int src, int num) {
        for(int i = 0; i < num; i++){
            for(int j = 0; j < width; j++) {
                data[size + i + j * capacity] = queue[(src + i) * width + j] ;
            }
        }
        size += num;
    }
};


struct RayQueue {
    float *data;
    int size, capacity, width;
 
    RayQueue(int capacity, int width): capacity(capacity), width(width){
        data = new float[capacity * width];
        size = 0;
    }

    ~RayQueue(){
        delete[] data;  
    }

  // queue store is different from orign ray

    void put(float *rays, int src, int num, int rays_capacity) {
        for(int i = 0; i < num; i++){
            for(int j = 0; j < width; j++) {
                data[(size + i) * width + j]  = rays[src + i + j * rays_capacity];
            }
        }
        size += num;
    }

    void copy(struct RayQueue* a, int n){
        for(int i = 0; i < width; i++){
            data[size * width + i] = a->data[n * width + i];
        }
        size ++;
    }

    float& operator[](int id) const {
        return data[id]; 
    }

    bool isfull() { return size == capacity; }

    bool isempty() {return size <= 0;}

    float* rays() {return data; }

    int* ids() {return (int*)data;}

    float* datas() {return &data[capacity];}

    void clear() {size = 0;}

    int get_size() { return size;}

    int get_capacity() {return capacity;}
};

inline void swap(struct RayQueue* &a, struct RayQueue* &b) {
    struct RayQueue* tmp;
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
         
inline void swap(float *&a, float *&b) {
    float* tmp;
    tmp = a;
    a = b;
    b = tmp;
}

