#include "mpi.h"

#define MSG_SIZE 3

struct Communicator {
    MPI_Status  sta[3];
    MPI_Request req[3];
    MPI_Comm Client_Comm;
    int rank, size, master;
    int group_rank, group_size;
    int width;
    bool first;

    Communicator(int width);
    
    ~Communicator(); 
    void Reduce_image(float* film, float *reduce_buffer, int pixel_num, bool server);
    void Send_rays(struct RayQueue* buffer, int size, int dst); 
    void Isend_rays(struct RayQueue* buffer, int size, int dst, int tag); 
    void Wait(int tag);
    int Recv_rays(struct RayQueue* buffer, int src);

    // msg[0] where msg from; msg[1] 1 has work 0 no work
    void Send_request(int id, bool has_work, int dst, int request_size);
    void Send_noray(int dst);
    void Send_end(int dst);
    bool Server_recv(struct RayQueue* buffer, int *msg); 
};


