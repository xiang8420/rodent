#pragma once

#include <vector>
#include <mutex>

struct Rays {
    std::vector<float> queue;
    int size, capacity, logic_width, store_width;
    bool compact; 
    bool *mask;

    Rays(int capacity, int width, bool compact); 

    ~Rays(){queue.clear();}

    int check_capacity(int);

    void read_device_buffer(float *rays, int src, int num, int rays_capacity, int rank);

    int write_to_device_buffer(Rays * , int, int, int, bool ); 
    
    void read_from_rays(struct Rays* a, int st);
    
    void read_from_ptr(char *src_ptr, int copy_size);

    float& operator[](int id) { return queue[id];}

    bool full() { return size == capacity; }

    bool empty() {return size <= 0;}

    float* get_data() {return queue.data(); }

    int clear() {
        int s = size;
        size = 0;
        return s;
    }

    int get_size() { return size;}

    int get_capacity() {return capacity;}
};

struct RayList {
    struct Rays *primary;
    struct Rays *secondary;
    size_t logic_capacity, store_capacity;
    std::string type;
    std::mutex mutex;

    RayList(int n, std::string type, bool compact):logic_capacity(n), type(type) {
        store_capacity = (n & ~((1 << 5) - 1)) + 32; // round to 32
        primary   = new struct Rays(store_capacity, 21, compact);
        secondary = new struct Rays(store_capacity, 14, compact);
    }
    
    ~RayList() {
        delete primary;
        delete secondary;
    }

    struct Rays* get_primary() {return primary;}
    struct Rays* get_secondary() {return secondary;}

    int primary_size(){return primary->size;}    
    
    int secondary_size(){return secondary->size;}    
   

    int size() { return primary->size + secondary->size; }
    bool empty() { return primary->empty() && secondary->empty(); } 

    void clear() {
        primary->size = 0 ;
        secondary->size = 0;
    }

    void lock() {mutex.lock();}
    void unlock() {mutex.unlock();}
    
    static void read_from_device_buffer(RayList ** raylists, float *raybuffer, size_t size, size_t capacity, bool primary, int rank);
    
    void write_to_device_buffer(RayList *, int);
    
    void read_from_message(char*, int, int, int); 
    
    static void read_from_message(RayList **, char*, int, int);
     
    static void classification(RayList ** raylists, Rays *raybuffer);
};

