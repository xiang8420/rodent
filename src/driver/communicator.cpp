#include "communicator.h"
#include "rayqueue.h"

Communicator::Communicator() {
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    master = size - 1;
    first = true;
    
    int group = rank == master ? 1 : 0;
    MPI_Comm_split(MPI_COMM_WORLD, group, rank, &Client_Comm);
    MPI_Comm_size(Client_Comm, &group_size);
    MPI_Comm_rank(Client_Comm, &group_rank);
}
    
Communicator::~Communicator() {
    MPI_Finalize();
}
    
void Communicator::Send_rays(struct RayQueue* buffer, int size,  int dst, bool has_work, bool primary) {
//    printf("send rays %d %d\n", rank, dst);
    int* a = (int*)buffer->data;
    int width = buffer->width;
    int buffer_size = buffer->get_size();
    int msg[MSG_SIZE];
    msg[0] = has_work ? 2 : 1; // -1 send ray 
    msg[1] = rank;
    msg[2] = primary? buffer_size : 0;
    msg[3] = !primary? buffer_size : 0;
    float *data = buffer->rays();
    int    send_size = buffer_size > size ? size : buffer_size;
    MPI_Send(&msg[0], MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD);
    MPI_Send(&data[(buffer_size - send_size) * width], send_size * width, MPI_FLOAT, dst, 1, MPI_COMM_WORLD);
    buffer->size -= send_size;
}

void Communicator::Isend_rays(struct RayQueue* buffer, int size, int dst, int tag) {
    int buffer_size = buffer->get_size();
    int width = buffer->width;
    float *data = buffer->rays();
    int    send_size = buffer_size > size ? size : buffer_size;
 //   printf("send ray send size %d buffer_size%d size %d\n", send_size, buffer_size, size);    
    MPI_Isend(&data[(buffer_size - send_size) * width], send_size * width, MPI_FLOAT, dst, 1, MPI_COMM_WORLD, &req[tag]);
    buffer->size -= send_size;
}

void Communicator::Wait(int tag){
     MPI_Wait(&req[tag],sta);
}

// need mark for secondaty or primary
int Communicator::Recv_rays(struct RayQueue* buffer, int src) {
//    printf("recv rays %d\n", rank);
    // return 0 end, 1 no ray, 2 get rays
    int msg[MSG_SIZE];
    int recv_num;
    MPI_Recv(msg, MSG_SIZE, MPI_INT, src, 1, MPI_COMM_WORLD, sta);
    if(msg[0] == -1) {
        return -1;
    } else if (msg[1] > 0) {
        int src = msg[1];
        int width, size; 
        if (msg[2] > 0) {
            width = 21; size = msg[2];
        } else {
            width = 14; size = msg[3];
        }
        float *rays = &buffer->rays()[buffer->get_size() * width];
        buffer->size += size;
        MPI_Recv(rays, size * width, MPI_FLOAT, src, 1, MPI_COMM_WORLD, sta); 
//        printf("get rays from %d recv_num %d size %d\n", src, recv_num / width, buffer->get_size());
        return 2;
    }
}

// msg[0] message status [1] src [2] primary size [3] secondary size  
void Communicator::Send_request(bool has_work, int dst, int request_size, bool primary){
    int msg[MSG_SIZE];
    msg[0] = has_work ? 0 : -1;
    msg[1] = rank;
    msg[2] = primary? request_size : 0;
    msg[3] = !primary? request_size : 0;
//    printf("send request %d %d %d %d \n", msg[0], msg[1], msg[2], msg[3]);
    MPI_Send(&msg[0], MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD);
}

void Communicator::Send_noray(int dst) {
//    printf("send norays %d %d\n", rank, dst);
    int msg[MSG_SIZE];
    msg[0] = 0;
    MPI_Send(&msg[0], MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD);
}

void Communicator::Send_end(int dst) {
//    printf("send end %d %d\n",rank, dst);
    int msg[MSG_SIZE];
    msg[0] = -1;
    MPI_Send(&msg[0], MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD);
}

bool Communicator::Server_recv(struct RayQueue* buffer, int *msg) {
//    printf("server recv  %d\n",rank);
    int recv_num, src;
    MPI_Recv(msg, MSG_SIZE, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, sta);
    if(msg[0] > 0) {
//        printf("server recv rays %d %d %d %d\n", msg[0], msg[1], msg[2], msg[3]);
        if(msg[2] > 0) {
            buffer->width = 21;
            buffer->size = msg[2];
            MPI_Recv(buffer->rays(), msg[2] * 21, MPI_FLOAT, msg[1], 1, MPI_COMM_WORLD, sta);
        } 
        if(msg[3] > 0) {
            buffer->width = 14;
            buffer->size = msg[3];
            MPI_Recv(buffer->rays(), msg[3] * 14, MPI_FLOAT, msg[1], 1, MPI_COMM_WORLD, sta);
        } 
        return false;
    } else {
        return true;
    }
}

void Communicator::Reduce_image(float* film, float *reduce_buffer, int pixel_num, bool server){
    if(server){
        MPI_Reduce(film, reduce_buffer, pixel_num, MPI_FLOAT, MPI_SUM, 0, Client_Comm); 
    } else {
        MPI_Reduce(film, reduce_buffer, pixel_num, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
    }
}

void Communicator::Send_msg(int dst, int* msg){
    MPI_Send(&msg[0], MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD);
}

void Communicator::Recv_msg(int dst, int *msg) {
    MPI_Recv(msg, MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD, sta);
}

void Communicator::Recv_rays(int src, bool primary, int recv_size, struct RayQueue* buffer) {
    if(recv_size == 0) return;
    int width = primary ? 21 : 14;
    printf("recv from %d %d %d buffer size%d\n", src, recv_size, width, buffer->size);
    buffer->width = width;
    float *rays = &buffer->rays()[buffer->get_size() * width];
    buffer->size += recv_size;
    MPI_Recv(rays, recv_size * width, MPI_FLOAT, src, 1, MPI_COMM_WORLD, sta); 
    
//    int* ids = (int*)rays;
//    for(int i = 0; i < 10; i ++) {
//        printf("|r %d %d ", ids[i * width], ids[i * width + 9]);
//    }
//    printf("\n");
}

void Communicator::Send_rays(int dst, bool primary, int send_size, struct RayQueue* buffer) {
    if(send_size == 0) return;
    int width = primary ? 21 : 14;
    float *rays = &buffer->rays()[(buffer->get_size() - send_size) * width];
    buffer->size -= send_size;
    printf("send to %d %d %d buffer size%d\n", dst, send_size, width, buffer->size);
    int* ids = (int*)rays;
//    for(int i = 0; i < 10; i ++) {
//        printf("|s %d %d ", ids[i * width], ids[i * width + 9]);
//    }
//    printf("\n");
    MPI_Send(rays, send_size * width, MPI_FLOAT, dst, 1, MPI_COMM_WORLD);
}

