#include <cstdlib>
#include <cstring>
#include "raylist.h"
#include "message.h"

Message::Message() {
    printf("construct message void\n");
    content.reserve(1024 * 1024 * 21 * 3); 
}

Message::Message(MPI_Status &status, MPI_Comm comm) {
    int count;
    MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &count);
    MPI_Status s0;
    header = new MessageHeader();
    if (count == sizeof(*header)) {
        printf("only header %d\n", count);
        MPI_Recv((char *)header, count, MPI_UNSIGNED_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm, &s0);
        printf("Message count %d root %d collective %d\n", count, header->broadcast_root, header->collective);
        if(header->collective)
            printf("is collective\n");
    } else {
        char *buffer = (char *)malloc(count);
        printf("content\n");
        printf("before recv \n");
        MPI_Recv(buffer, count, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm, &s0);
        memcpy(header, buffer, sizeof(*header));
        printf("recv Message::Message %d\n", header->content_size);
        content.resize(header->content_size);
        memcpy(content.data(), buffer+sizeof(*header), header->content_size);
        printf("Message copy success\n");
        free(buffer);
    }
}

void Message::deserialize(struct RayList* inList) {
    char* ptr = content.data();
    ///////////////ptr = content.data
    printf("message primary %d, secondary %d\n", header->primary, header->secondary);
    if(header->primary > 0 ) {
        Rays *primary = inList->get_primary();
        char* primary_copy_ptr = (char*)(primary->get_data() + primary->size * primary->store_width * sizeof(float));
        int   primary_length = header->primary * inList->get_primary()->store_width * sizeof(float);
        memcpy(primary_copy_ptr, ptr, primary_length);
        primary->size += header->primary; 
        
        int width = primary->store_width; 
        int* ids = (int*)(primary->get_data());
        printf("deserialize size %d %d:", header->primary, header->secondary);
        for(int i = 0; i < 5; i ++) {
            printf("%d %d*", ids[i*width], ids[i*width + 9]);
        }
        printf("\n");
        ptr += primary_length;
    } 
    if(header->secondary > 0) { 
        Rays *secondary = inList->get_secondary();
        char* secondary_copy_ptr = (char*)(secondary->get_data() + secondary->size * secondary->store_width * sizeof(float));
        int secondary_length = header->secondary * inList->get_secondary()->store_width * sizeof(float);
        printf("secondary copy \n");
        memcpy(secondary_copy_ptr, ptr, secondary_length);
        printf("after secondary coyp \n");
        secondary->size += header->secondary; 
        
        int width = secondary->store_width; 
        int* ids = (int*)(secondary->get_data());
        printf("deserialize:");
        for(int i = 0; i < 5; i ++) {
            printf("%d %d*", ids[i*width], ids[i*width + 9]);
        }
        printf("\n");

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

RayMsg::RayMsg(RayList* outList, int src, int dst, int chunk, bool idle) {
    printf("construct message ray\n");
    header = new MessageHeader(-1, src, 0, false, 0, idle, chunk);
    destination = dst;
    printf("serialize %d\n", chunk);
    header->content_size = serialize(outList);
    printf("end serialize %d\n", get_chunk());
    outList->clear();
}

QuitMsg::QuitMsg(int src) {
    printf("construct message collective\n");
    header = new MessageHeader(src, -1, 0, true, 0, true, -1);
    printf("new Message root %d collective %d\n", header->broadcast_root, header->collective);
    destination = -1;
}

StatusMsg::StatusMsg(int src, int* status, int proc_size) {
    printf("construct message status\n");
    header = new MessageHeader(src, -1, 0, false, 0, true, -1);
   
    int status_length = proc_size * proc_size * sizeof(int);
    content.resize(status_length);
    memcpy(content.data(), status, status_length);
    header->content_size = status_length;
    
    printf("new Message root %d collective %d\n", header->broadcast_root, header->collective);
    destination = -1;
}
