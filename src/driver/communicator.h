#include "mpi.h"
#include "message.h"
#include "process_status.h"
#include <fstream>
#define MSG_SIZE 5

struct mpi_send_buffer {
// If p2p, use left
    MPI_Request lrq;
    MPI_Request rrq;

    int n;
    char *send_buffer;
    int total_size;
};

struct Communicator {
    MPI_Status  sta[3];
    MPI_Request req[3];
    MPI_Request lrq, rrq;

    MPI_Comm Client_Comm;
    int rank, size, master;
    int group_rank, group_size;
    bool first;
    bool pause, quit;
    std::ofstream os;
   
    int send_ray_count, recv_ray_count;
    int send_msg_count, recv_msg_count;
        

    std::mutex mutex;
    std::condition_variable cond; // primary, secondary buffer size < max

    std::vector<mpi_send_buffer*> mpi_in_flight;

    Communicator();
    ~Communicator(); 
    
    void lock() { mutex.lock(); }
    void unlock() { mutex.unlock(); }
    void signal() { cond.notify_all(); }
    void wait() {
        std::unique_lock <std::mutex> lk (mutex);
        cond.wait(lk); 
    } 

    void reduce_image(float* film, float *reduce_buffer, int pixel_num, bool server);
    
    void Isend_rays(struct Rays* buffer, int size, int dst, int tag); 
    void mpi_wait(int tag);
    int recv_rays(struct Rays* buffer, int src);
    
    // msg[0] where msg from; msg[1] 1 has work 0 no work
    void send_message(int dst, int msg[MSG_SIZE]); 
    void send_request(bool has_work, int dst, int request_size, bool primary);
    void send_noray(int dst);
    void send_end(int dst);

    void send_msg(int dst, int* msg);
    void recv_msg(int dst, int *msg);
    void recv_rays(int src, int recv_size, struct Rays* buffer);
    void send_rays(int dst, int send_size, struct Rays* buffer);

    // broadcast or p2p. return sent number 1, in p2p, 0, 1 or 2 in bcast
    int  Export(Message * m, ProcStatus *rs); 

    bool recv_message(RayList** list, ProcStatus *rs); 
    
    bool send_message(Message* msg, ProcStatus *rs); 

    void collective(ProcStatus *rs); 

    bool barrier(struct RayList* out);

    void purge_completed_mpi_buffers(); 

};



