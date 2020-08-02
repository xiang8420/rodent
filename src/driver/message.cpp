#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include "message.h"

Message::Message() {
    printf("construct message void\n");
    content.reserve(1024 * 1024 * 21 * 3);
   //??? 
}
//recv msg
RecvMsg::RecvMsg(MPI_Status &status, MPI_Comm comm) {
    int count;
    MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &count);
    MPI_Status s0;
    header = new MessageHeader();
    if (count == sizeof(*header)) {
        MPI_Recv((char *)header, count, MPI_UNSIGNED_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm, &s0);
    } else {
        char *buffer = (char *)malloc(count);
        MPI_Recv(buffer, count, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm, &s0);
        memcpy(header, buffer, sizeof(*header));
        content.resize(header->content_size);
        memcpy(content.data(), buffer+sizeof(*header), header->content_size);
        free(buffer);
    }
}

//recv message, if it's ray msg load it to list 
RecvMsg::RecvMsg(RayList** List, MPI_Status &status, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    std::ofstream os; 
    os.open("out/proc_buffer_" + std::to_string(rank), std::ios::out | std::ios::app ); 
    os<<"master read from message\n"; 
    printf("construct message recv ray from proc %d\n", status.MPI_SOURCE);
    
    int count;
    MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &count);
    MPI_Status s0;
    header = new MessageHeader();
    if (count == sizeof(*header)) {
        MPI_Recv((char *)header, count, MPI_UNSIGNED_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm, &s0);
        if(header->collective)
            os<<"recv collective\n";
    } else {
        char *buffer = (char *)malloc(count);
        MPI_Recv(buffer, count, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm, &s0);
        memcpy(header, buffer, sizeof(*header));
        printf("recv msg size count %d \n", count);
        os<<"recv msg sender "<<header->sender <<" root "<<header->root<<"\n";
        if(header->primary == 0 && header->secondary ==0 ) {
            //status
            os<<"recv status "<<header->content_size<<"\n";
            content.resize(header->content_size);
            memcpy(content.data(), buffer+sizeof(*header), header->content_size);
            os<<"recv status over\n";
        } else {
            //rays 
            if(header->chunk == -1) {
                os<<"recv mixed rays "<< header->sender <<" "<<header->content_size<<"\n";
                RayList::read_from_message(List, (char*)buffer+sizeof(MessageHeader), header->primary, header->secondary);
            } else {
                os<<"recv normal rays "<< header->sender <<" "<<get_ray_size()<<"\n";
                RayList *in = List[header->chunk];
                in->read_from_message((char*)buffer+sizeof(MessageHeader), header->primary, header->secondary, rank);
            }
        }
        free(buffer);
    }
}

int Message::serialize(struct RayList* outList) {
    if(outList->empty()) return 0;
  
    Rays* primary = outList->get_primary(); 
    int width = primary->store_width; 
    int* ids = (int*)(primary->get_data());
    for(int i = 0; i < 5; i ++) {
        printf(" %d %d$", ids[i*width], ids[i*width + 9]);
    }
    printf("\n");


    header->primary = outList->get_primary()->size; 
    header->secondary = outList->get_secondary()->size; 

    printf("serialize size %d %d:", header->primary, header->secondary);
    
    int primary_length = header->primary * outList->get_primary()->store_width * sizeof(float);
    int secondary_length = header->secondary * outList->get_secondary()->store_width * sizeof(float);
    content.resize(primary_length + secondary_length);
    printf("Semd Message size %d\n", content.size());  

    char* ptr = content.data();
    memcpy(ptr, outList->get_primary()->get_data(), primary_length);
    memcpy(ptr + primary_length, outList->get_secondary()->get_data(), secondary_length);
    
    return primary_length + secondary_length;
}

RayMsg::RayMsg(RayList** List, int src, int dst, int chunk_size, bool idle) {
    std::ofstream os; 
    os.open("out/proc_buffer_worker", std::ios::out | std::ios::app ); 
    
    header = new MessageHeader(-1, src, MsgType::Ray, false, 0, idle, -1);
    destination = dst;
    
    os<<"RayMsg :\n";
    for(int i = 0; i < chunk_size; i++) {
        if(List[i]->type == "out") {
            header->primary += List[i]->primary_size();
            header->secondary += List[i]->secondary_size();
            os<<"list "<<i<<"primary copy to "<<header->primary<<" secondary"<<header->secondary<<"\n";
        } 
    }
    int primary_length = header->primary * List[0]->get_primary()->store_width * sizeof(float);
    printf("new RayMsg primary %d secondary %d\n", header->primary, header->secondary);
    int secondary_length = header->secondary * List[0]->get_secondary()->store_width * sizeof(float);
    content.resize(primary_length + secondary_length);
    header->content_size = primary_length + secondary_length; 
    
    char* primary_ptr   = content.data();
    char* secondary_ptr = primary_ptr + primary_length;

    for(int i = 0; i < chunk_size; i++) {
        if(List[i]->type == "out" && !List[i]->empty()) {
            RayList * out = List[i];
            int length = out->primary_size() * out->get_primary()->store_width * sizeof(float);
            memcpy(primary_ptr,   out->get_primary()->get_data(), length);
            primary_ptr += length;
            length = out->secondary_size() * List[i]->get_secondary()->store_width * sizeof(float);
            memcpy(secondary_ptr, out->get_secondary()->get_data(), length);
            secondary_ptr += length;  
            out->clear();
        } 
    }
}

RayMsg::RayMsg(RayList* outList, int src, int dst, int chunk, bool idle) {
    printf("construct message ray\n");
    header = new MessageHeader(-1, src, MsgType::Ray, false, 0, idle, chunk);
    destination = dst;
//    Rays* primary = outList->get_primary(); 
//    int width = primary->store_width; 
//    int* ids = (int*)(primary->get_data());
//    for(int i = 0; i < 5; i ++) {
//        printf(" %d %d$", ids[i*width], ids[i*width + 9]);
//    }
//    printf("\n");


    header->primary = outList->get_primary()->size; 
    header->secondary = outList->get_secondary()->size; 
    
    int primary_length = header->primary * outList->get_primary()->store_width * sizeof(float);
    int secondary_length = header->secondary * outList->get_secondary()->store_width * sizeof(float);
    content.resize(primary_length + secondary_length);
    printf("Semd Message size %d\n", content.size());  

    char* ptr = content.data();
    memcpy(ptr, outList->get_primary()->get_data(), primary_length);
    memcpy(ptr + primary_length, outList->get_secondary()->get_data(), secondary_length);
    
    header->content_size = primary_length + secondary_length; 
    
    outList->clear();
}

QuitMsg::QuitMsg(int src) {
    printf("construct message quit \n");
    header = new MessageHeader(src, -1, MsgType::Quit, true, 0, true, -1);
    printf("new Message root %d collective %d\n", header->root, header->collective);
    destination = -1;
}

StatusMsg::StatusMsg(int src, int dst, int* status, int chunk, int proc_size) {
    printf("construct message status src %d\n", src);
    if(dst == -1) 
        header = new MessageHeader(src, -1, MsgType::Status, false, 0, true, chunk);
    else 
        header = new MessageHeader(-1, src, MsgType::Status, false, 0, true, chunk);

    destination = dst;
   
    int length = proc_size * proc_size * sizeof(int);
    content.resize(length);
    memcpy(content.data(), status, length);
    header->content_size = length;
    
    printf("new Message root %d collective %d\n", header->root, header->collective);
}

//broad cast schedule
ScheduleMsg::ScheduleMsg(int src, int* chunkStatus, int chunk, int chunk_size) {
    printf(" construct message schedule src %d\n", src);
    header = new MessageHeader(src, -1, MsgType::Schedule, false, 0, false, chunk);
    destination = -1;
    int length = chunk_size * sizeof(int);
    content.resize(length);
    memcpy(content.data(), chunkStatus, length);
    header->content_size = length;
    printf("new Message root %d collective %d\n", header->root, header->collective);
}
