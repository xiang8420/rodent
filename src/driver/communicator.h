#include "mpi.h"

#define MSG_SIZE 5

struct Communicator {
    MPI_Status  sta[3];
    MPI_Request req[3];
    MPI_Comm Client_Comm;
    int rank, size, master;
    int group_rank, group_size;
    bool first;

    Communicator();
    ~Communicator(); 
    void Reduce_image(float* film, float *reduce_buffer, int pixel_num, bool server);
    
    void Send_rays(struct RayQueue* buffer, int size, int dst, bool has_work, bool primary); 
    void Isend_rays(struct RayQueue* buffer, int size, int dst, int tag); 
    void Wait(int tag);
    int Recv_rays(struct RayQueue* buffer, int src);
    
    // msg[0] where msg from; msg[1] 1 has work 0 no work
    void Send_message(int dst, int msg[MSG_SIZE]); 
    void Send_request(bool has_work, int dst, int request_size, bool primary);
    void Send_noray(int dst);
    void Send_end(int dst);
    bool Server_recv(struct RayQueue* buffer, int *msg); 

    void Send_msg(int dst, int* msg);
    void Recv_msg(int dst, int *msg);
    void Recv_rays(int src, bool primary, int recv_size, struct RayQueue* buffer);
    void Send_rays(int dst, bool primary, int send_size, struct RayQueue* buffer);

};


