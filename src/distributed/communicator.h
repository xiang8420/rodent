#include "mpi.h"
#include "Message.h"
#include "ProcStatus.h"
#include <fstream>

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
    std::ofstream os;
   
    int send_ray_count, recv_ray_count;
    int send_msg_count, recv_msg_count;
        
    std::vector<mpi_send_buffer*> mpi_in_flight;

    Communicator();
    ~Communicator(); 
   
    int get_comm_size() { return size; } 
    int get_comm_rank() { return rank; }
    int isMaster(){return rank == master;} 

    void reduce_image(float* film, float *reduce_buffer, int pixel_num);

    void all_gather_float(float *a, float *res, int size); 

    // broadcast or p2p. return sent number 1, in p2p, 0, 1 or 2 in bcast
    int  Export(Message * m, ProcStatus *rs); 

    bool recv_message(RayList** List, RayStreamList * inList, ProcStatus *ps); 
    
    void send_message(Message* msg, ProcStatus *rs); 

    void purge_completed_mpi_buffers(); 

};



