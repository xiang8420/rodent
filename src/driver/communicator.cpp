#include "communicator.h"
#include "rayqueue.h"
#define MSG_SIZE 3

Communicator::Communicator(int width): width(width) {
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
    
void Communicator::Send_rays(struct RayQueue* buffer, int size, int dst) {
//    int* a = (int*)buffer->data;
//    for(int i = 0; i < 10; i ++) {
//        printf("%d %d *", a[i * width], a[i * width + 10]);
//    } 
//    printf("\n");
    int buffer_size = buffer->get_size();
    float *data = buffer->rays();
    int    send_size = buffer_size > size ? size : buffer_size;
//    printf("send ray send size %d buffer_size%d size %d", send_size, buffer_size, size);    
    MPI_Send(&data[(buffer_size - send_size) * width], send_size * width, MPI_FLOAT, dst, 1, MPI_COMM_WORLD);
    buffer->size -= send_size;
}

void Communicator::Isend_rays(struct RayQueue* buffer, int size, int dst, int tag) {
    int buffer_size = buffer->get_size();
    float *data = buffer->rays();
    int    send_size = buffer_size > size ? size : buffer_size;
 //   printf("send ray send size %d buffer_size%d size %d\n", send_size, buffer_size, size);    
    MPI_Isend(&data[(buffer_size - send_size) * width], send_size * width, MPI_FLOAT, dst, 1, MPI_COMM_WORLD, &req[tag]);
    buffer->size -= send_size;
}

void Communicator::Wait(int tag){
     MPI_Wait(&req[tag],sta);
}

int Communicator::Recv_rays(struct RayQueue* buffer, int src) {
    // return 0 end, 1 no ray, 2 get rays
    int msg[3];
    int recv_num;
    MPI_Probe(src, 1, MPI_COMM_WORLD, sta);
    MPI_Get_count(&sta[0], MPI_FLOAT, &recv_num);
    if(recv_num == 1) {
        MPI_Recv(&msg, 1, MPI_FLOAT, src, 1, MPI_COMM_WORLD, sta);        
        return 0;
    } else if (recv_num == 2) {
        MPI_Recv(&msg, 2, MPI_FLOAT, src, 1, MPI_COMM_WORLD, sta);        
        return 1;
    } else {
        float *rays = &buffer->rays()[buffer->get_size() * width];
        MPI_Recv(rays, recv_num, MPI_FLOAT, src, 1, MPI_COMM_WORLD, sta); 
        buffer->size += recv_num / width;
//        printf("get rays from %d recv_num %d size %d\n", src, recv_num / width, buffer->get_size());
        return 2;
    }
}

// msg[0] where msg from; msg[1] 1 has work 0 no work
void Communicator::Send_request(int id, bool has_work, int dst, int request_size){
    int msg[3];
    msg[0] = has_work ? 1 : 0;
    msg[1] = id;
    msg[2] = request_size;
    MPI_Send(&msg[0], MSG_SIZE, MPI_INT, dst, 1, MPI_COMM_WORLD);
}

void Communicator::Send_noray(int dst) {
    float a[2];
    MPI_Send(&a, 2, MPI_FLOAT, dst, 1, MPI_COMM_WORLD);
}


void Communicator::Send_end(int dst) {
    float a = 0;
    MPI_Send(&a, 1, MPI_FLOAT, dst, 1, MPI_COMM_WORLD); 
}

bool Communicator::Server_recv(struct RayQueue* buffer, int *msg) {
    int recv_num, src;
    MPI_Probe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, sta);
    MPI_Get_count(&sta[0], MPI_INT, &recv_num);
    src = sta[0].MPI_SOURCE;
    if(recv_num == 3) {
//        printf("master recv request\n\n");
        MPI_Recv(msg, 3, MPI_INT, src, 1, MPI_COMM_WORLD, sta);
        return true;
    } else {
        MPI_Get_count(&sta[0], MPI_FLOAT, &recv_num);
        MPI_Recv(buffer->rays(), recv_num, MPI_FLOAT, src, 1, MPI_COMM_WORLD, sta);
//    int* a = (int*)buffer->data;
//    for(int i = 0; i < 10; i ++) {
//        printf("%d %d *", a[i * width], a[i * width + 10]);
//    } 
//    printf("\n");
        buffer->size = recv_num / width;
        return false;
    }
}

void Communicator::Reduce_image(float* film, float *reduce_buffer, int pixel_num, bool server){
    if(server){
        MPI_Reduce(film, reduce_buffer, pixel_num, MPI_FLOAT, MPI_SUM, 0, Client_Comm); 
    } else {
        MPI_Reduce(film, reduce_buffer, pixel_num, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
    }
}
