#include <vector>
#include <memory>
#include <dirent.h>
#include <unistd.h>
#include "bbox.h"
#include "obj.h"

#include <sys/file.h>
#ifdef WIN32
#include <direct.h>
#define create_directory(d) _mkdir(d)
#else
#include <sys/stat.h>
#define create_directory(d) { umask(0); mkdir(d, 0777); }
#endif

#define MAX_CHUNK_NUM 8
#define CAPACITY 1048608
struct Ray_Queue {
    std::vector<float> *ray;
    int num;
    Ray_Queue(){}
    void setup(std::vector<float> *r, int n) {
        ray = r;
        num = n;
    }
	bool empty() {
		return num == 0;
	}
};

void swap(struct Ray_Queue& a, struct Ray_Queue& b){
    Ray_Queue tmp;
    tmp = a;
    a = b;
    b = tmp;
}

void create(std::string path) {
	std::fstream of;
	of.open(path.c_str(), std::ios::out);
	if(!of){
		of << "a "; 
	}
	of.close();		
}

struct Scheduler {
    //only mpi 0 can create
    int mpi_id, mpi_size, num_chunks, chunk_id;
    bool mpi = false;
    std::string lock_path;
    std::vector<float> rays[MAX_CHUNK_NUM];
    struct Ray_Queue first, second;
    int    num_saved_rays;  //num of saved ray files
    int    num_rays[MAX_CHUNK_NUM];

    Scheduler(int mpi_id, int mpi_size, int num_chunks) : mpi_id(mpi_id), mpi_size(mpi_size), num_chunks(num_chunks) {
        num_saved_rays = 0;
        ////////// 
        //Temporarily set chunk id to mpi id
        chunk_id = mpi_id;
        if(mpi_id == 0){
            std::remove("data/rays");
            create_directory("data/rays");
            for(int i = 0; i < num_chunks; i++) {
                std::string path("data/rays/r" + std::to_string(i));
                create_directory(path.c_str());
                std::vector<float> lock(1);
                lock_path = "data/rays/lock" + std::to_string(i);
                write_buffer(lock_path.c_str(), lock);
            }
            create_directory("data/rays/rn/");
        }
        lock_path = "data/rays/lock" + std::to_string(chunk_id);
        for(int i = 0; i<= num_chunks; i++) {
            num_rays[i] = 0;
            rays[i].reserve(CAPACITY * 20);
            rays[i].resize(CAPACITY * 20);
            rays[i].push_back(0);
        }
        first.setup(&rays[chunk_id], 0);
        second.setup(&rays[num_chunks], 0);
        create("data/rays/syn" + std::to_string(chunk_id));
    }

    
    void reset() {
        num_saved_rays = 0;
        for(int i = 0; i <= num_chunks; i++) {
            num_rays[i] = 0;
            rays[i].resize(CAPACITY * 20);
            rays[i].push_back(0);
        }
        first.setup(&rays[chunk_id], 0);
        second.setup(&rays[num_chunks], 0);
        create("data/rays/syn" + std::to_string(chunk_id));
    }

    void copy_rays(float* a, float* b, int src, int dst, int copy_num, size_t capacity) {
//        printf("dst %d src %d \n", dst, src);
        for(int i = 0; i < 20; i++) {
            for(int j = 0; j < copy_num; j++){
                b[dst + j + i * capacity] = a[src + j + i * capacity];
            }
        }
    }
    
    bool mpi_rays_transfer(float *primary) {
        return false; 
    }
    bool send_rays(int chunk) {
        printf("send rays %d %d\n", chunk, num_rays[chunk]);
        rays[chunk].resize(CAPACITY * 20);
        FILE *lock = fopen(("data/rays/lock" + std::to_string(chunk)).c_str(), "r+b");
        flock(lock->_fileno, LOCK_EX);
        
        write_buffer("data/rays/r" + std::to_string(chunk) + "/" 
                                   + std::to_string(chunk_id)  + "_"
                                   + std::to_string(num_saved_rays)  
                                   + ".bin", rays[chunk]);
        if(chunk_id == 1) {
              write_buffer("data/rays/rn/"
               + std::to_string(chunk_id)  + "_"
               + std::to_string(num_saved_rays)
               + ".bin", rays[chunk]);
        } 

        printf("send ray %s %d %d\n", ("data/rays/r" + std::to_string(chunk) + "/"
                                       + std::to_string(chunk_id)  + "_"
                                       + std::to_string(num_saved_rays)
                                       + ".bin").c_str(), chunk, num_rays[chunk]);
        fclose(lock);
        flock(lock->_fileno, LOCK_UN);
        printf("chunk %d rays %d \n", chunk, num_rays[chunk]);
        num_rays[chunk] = 0;
        num_saved_rays++; 
    }
    
    bool recv_rays(struct Ray_Queue &queue, int capacity){
        DIR *dp;
        struct dirent *dirp;
        if((dp = opendir(("data/rays/r" + std::to_string(chunk_id)).c_str())) == NULL) {
            printf("Can't open dir\n");
        }
        while((dirp = readdir(dp)) != NULL) {
            if(dirp->d_type == 8){
                
                FILE *lock = fopen(lock_path.c_str(), "r+b");
                flock(lock->_fileno, LOCK_EX);
                
                std::string path("data/rays/r" + std::to_string(chunk_id) + "/" + std::string(dirp->d_name));
                printf("ray path %s\n", path.c_str());
                read_buffer(path, *(queue.ray)); 
                int *rays_id = (int*)&(*queue.ray)[0];
                for(int i = 0; i < capacity; i++) {               
                    if(rays_id[i] == -1) {
                        queue.num = i;
                        break;
                    }
                }
                if(queue.num == 0) { queue.num = capacity; }
                std::remove(path.c_str());
                fclose(lock);
                flock(lock->_fileno, LOCK_UN); 
                return true;
            } 
        }    
        return false;
    }
    
    bool any_waiting() {
        // add all_done file fast judge
        // //return if other node is done
        for(int i = 0; i < mpi_size; i ++) {
            if(i != chunk_id) {
                std::fstream of;
                of.open(("data/rays/syn" + std::to_string(i)).c_str(), std::ios::in);
                if(!of) {
                    return true;
                }
                of.close();
            }
        }
        return false;
    }    
    
    bool any_working() {
        // add all_done file fast judge
        //return if other node is done
        for(int i = 0; i < mpi_size; i ++) {
            if(i != chunk_id) {
                std::fstream of;
                of.open(("data/rays/syn" + std::to_string(i)).c_str(), std::ios::in);
                if(of) {
                    return true;
                }
                of.close();
            }
        }
        return false;
    }    

    int file_rays_transfer(float *primary, size_t size, size_t capacity) {
        //cluster rays if ray bins is full write to queue
        int  primary_size = 0;
        bool load_local = false;
        std::vector<int> record;
        int *rays_chunk_id = (int*)primary + 10 * capacity;
        int *rays_id = (int*)primary;

        printf("start ray cluster\n");

        for(int i = 0; i < size; i++) {
            int cid = rays_chunk_id[i] / 1000;
            if(cid == chunk_id){ 
                copy_rays(primary, primary, i, primary_size, 1, capacity);
                primary_size++;
            } else {
                if(cid < num_chunks && cid >= 0){
                    copy_rays(primary, &rays[cid][0], i, num_rays[cid]++, 1, capacity);
                    if(num_rays[cid] == capacity){
                        send_rays(cid);
                    }
                }
            }
        }
		if(first.empty()) {
            printf("empty read\n");
			recv_rays(first, capacity);
		}
		printf("0 %d primary size %d capacity %d  %d\n", chunk_id, primary_size, capacity, first.num);
		while( primary_size < capacity && !first.empty()) {			
			int sent_num = primary_size + first.num < capacity ? first.num : capacity - primary_size;
            printf("primary size %d send num %d first num %d\n", primary_size, sent_num, first.num);
            first.num -= sent_num;
            copy_rays(&(*first.ray)[0], primary, first.num, primary_size, sent_num, capacity);
            primary_size += sent_num;
            if(first.empty()) { recv_rays(first, capacity); }
		}
        
        printf("2 scheduler primary size %d first: %d\n", primary_size, first.num);
        printf("%d end ray transfer\n", chunk_id);
        std::string syn_path = "data/rays/syn" + std::to_string(chunk_id);
        
		printf("%d num_rays:", mpi_id);
        for(int i = 0; i < num_chunks; i++) {
            if(i == chunk_id)
                printf("%d ", primary_size);
            else
                printf("%d ", num_rays[i]);
        }
        
        printf("\n");
        if(any_waiting() || primary_size == 0) {
            for(int i = 0; i < num_chunks; i++){
                if(num_rays[i] !=0 && i != chunk_id) {
                    int *id = (int*)&rays[i][0];
                    id[num_rays[i]] = -1;
                    printf("\n\nprimary %d  %d %d\n\n", primary_size, i, num_rays[i]);
                    send_rays(i);
                }
            }
        }

        if(primary_size == 0) {
            remove(syn_path.c_str());
            printf("%d is waiting\n", chunk_id);
            while(1){
                usleep(300);
                if(any_working()){   
                    if(recv_rays(first, capacity))
                        break;
                } else {
                    if(recv_rays(first, capacity))
                        break;
                    return 0;
                }
            }
            printf("%d first num %d\n", chunk_id, first.num);
            copy_rays(&(*first.ray)[0], primary, 0, 0, first.num, capacity);
            primary_size = first.num;
            first.num = 0;
            create(syn_path);
        }
		return primary_size;
		
    }
    
    int mpi_rays_transfer(float *primary, size_t size, size_t capacity) {
        
    }
};

static std::unique_ptr<Scheduler> scheduler;

void setup_scheduler(size_t mpi_id, size_t mpi_num, size_t num_chunks) {
    scheduler.reset(new Scheduler(mpi_id, mpi_num, num_chunks));
}

void reset_scheduler(){
    scheduler->reset();
}

void cleanup_scheduler() {
    scheduler.reset();
}

int rays_transfer(float *primary, size_t size, size_t capacity){
    if(scheduler->mpi){
        return scheduler->mpi_rays_transfer(primary, size, capacity);
    } else {
        return scheduler->file_rays_transfer(primary, size, capacity);
    }
}


