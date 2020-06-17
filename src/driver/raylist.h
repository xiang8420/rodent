#include <vector>
#include <mutex>

struct Rays {
    std::vector<float> queue;
    std::mutex mutex;
    float* data;
    int size, capacity, logic_width, store_width;
    bool compact; 
    bool *mask;

    Rays(int capacity, int width, bool compact = true); 

    ~Rays(){queue.clear();}

    void read_dev_rays(float *rays, int src, int num, int rays_capacity);

    int write_dev_rays(float *rays, int num, size_t rays_capacity); 
    
    void copy_rays(struct Rays* a, int st);

    float& operator[](int id) { return data[id];}

    bool full() { return size == capacity; }

    bool empty() {return size <= 0;}

    float* get_data() {return data; }

    int clear() {
        int s = size;
        size = 0;
        return s;
    }

    int get_size() { return size;}

    int get_capacity() {return capacity;}

    void lock() {mutex.lock();}
    void unlock() {mutex.unlock();}
};

struct RayList {
    struct Rays *camera;
    struct Rays *primary;
    struct Rays *secondary;
    int capacity;
    std::string type;
    std::mutex mutex;

    RayList(int n, std::string type):capacity(n), type(type) {
        camera    = new struct Rays(capacity, 21);
        primary   = new struct Rays(capacity, 21);
        secondary = new struct Rays(capacity, 14);
    }
    
    ~RayList() {
        delete camera;
        delete primary;
        delete secondary;
    }

    struct Rays* get_primary() {return primary;}
    struct Rays* get_secondary() {return secondary;}

    int primary_size(){return primary->size;}    
    
    int secondary_size(){return secondary->size;}    
   
    int camera_ray_size() {return camera->size;}

    int size() { return primary->size + secondary->size + camera->size; }
    bool empty() { return primary->empty() && secondary->empty() && camera->empty(); } 

    void clear() {
        primary->size = 0 ;
        secondary->size = 0;
        camera->size = 0;
    }

    void lock() {mutex.lock();}
    void unlock() {mutex.unlock();}

    static void classification(RayList ** raylists, Rays *raybuffer);
    static void classification(RayList ** raylists, float *raybuffer, size_t size, size_t capacity, bool primary);
};

