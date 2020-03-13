#include "communicator.h"
#include "rayqueue.h"
#include "buffer.h"

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
    printf("%d <- %d |recv %d %d\n", rank, src, recv_size, width);
    buffer->width = width;
    float *rays = &buffer->rays()[buffer->get_size() * width];
    buffer->size += recv_size;

    MPI_Recv(rays, recv_size * width, MPI_FLOAT, src, 1, MPI_COMM_WORLD, sta); 

//    int out_size;
//    MPI_Status status;
//    MPI_Probe(src, 1, MPI_COMM_WORLD, &status);
//    MPI_Get_count(&status, MPI_CHAR, &out_size);
//    std::vector<char> in(out_size);
//    MPI_Recv(in.data(), out_size, MPI_CHAR, src, 1, MPI_COMM_WORLD, sta);
//    LZ4_decompress_safe(in.data(), (char*)rays, in.size(), recv_size * width * sizeof(float));

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
    printf("%d -> %d |send %d %d \n", rank, dst, send_size, width);

    MPI_Send(rays, send_size * width, MPI_FLOAT, dst, 1, MPI_COMM_WORLD);

//    std::vector<char> out;
//    compress(rays, send_size * width, out);
//    size_t in_size  = sizeof(rays[0]) * send_size * width;
//    size_t out_size = out.size();
//    printf("%d after compress %d\n", in_size, out_size);
//    MPI_Send(out.data(), out_size, MPI_CHAR, dst, 1, MPI_COMM_WORLD);
    
//    int* ids = (int*)rays;
//    for(int i = 0; i < 10; i ++) {
//        printf("|s %d %d ", ids[i * width], ids[i * width + 9]);
//    }
//    printf("\n");
}

