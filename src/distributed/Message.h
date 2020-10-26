#include <vector>
#include <mutex>
#include "RayList.h"
#include "RayStreamList.h"
#include <condition_variable>
#include "mpi.h"

enum MsgType {Status, Ray, Schedule, Quit};

struct MessageHeader {
    int      root;                   // will be -1 for point-to-point
    int      sender;                 // will be -1 for broadcast
    int      tag;
    MsgType  type; //0 status, 1 rays , 2 quit
    int      content_size;
    int      primary, secondary;    
    int      chunk;   //-1 means mixed rays
    bool     collective;
    bool     idle;      //if this process is idle
    int      destination;
//    int      padding[2]; 

    MessageHeader(){}

    MessageHeader(int root, int sender, int dst,  MsgType type, bool collective, int size, bool idle, int chunk, int tag) 
        : root(root), sender(sender), destination(dst), type(type), collective(collective), content_size(size), idle(idle), chunk(chunk), tag(tag)
    {
       primary   = 0; 
       secondary = 0; 
    }
    
    bool HasContent() { return content_size > 0; }
};

class Message { 
protected:

    MessageHeader *header;
    std::vector<char> content;

public:
    Message();
    
    Message(MessageHeader *h): header(h) { }

    // write content to a raylists
    int write_raylist(struct RayList *outlist);
  
    // read inlist to message 
    int read_raylist(struct RayList *inlist);

    int get_tag() {return header->tag;}

    bool exit_msg() {return header->type == 2; }
    bool ray_msg() {return header->type == 1; }
    bool status_msg() {return header->type == 0; }

	bool is_busy() { return !header->idle; }

    int  get_type() { return header->type; }
    
    void set_type(int t) { header->type = MsgType(t); }

    int get_chunk() { 
        if (header->chunk < 0)
            printf("header chunk == -1\n");
        return header->chunk;
    }

    void set_chunk(int t) { header->chunk = t; }
    
    char *get_content() { return content.data(); }
    
    size_t get_ray_size() {return header->primary + header->secondary;}

    size_t primary_size(){return header->primary;}
    size_t secondary_size(){return header->secondary;}

    size_t get_size() { return content.size(); }

    int get_root() { return header->root; }
    
    int get_sender() { return header->sender; }

    bool is_broadcast() { return header->root != -1; }
    
    bool isP2P() { return ! is_broadcast(); }

	bool has_content() { return header->HasContent(); }

    void set_destination(int i) { header->destination = i; }
    
    int get_destination() { return header->destination; }

    unsigned char *get_header() { return (unsigned char *)header; }
      
    int get_header_size() { return sizeof(*header); }

    //! Wait for a blocking message to be handled.  
    /*! This will return when the message has been sent, and if the message is a broadcast
     * message (that is, runs in the messaging thread), for the message's
     * action to happen locally.  
     * \warning only valid for blocking messages
     */
    void wait();

    //! is this message collective (i.e. synchrnoizing across processes)?
	bool is_collective() { return header->collective; }
    
    void deserialize(char* ptr, struct RayList* outList); 
    
    int serialize(struct RayList* outList); 
};

class RayMsg : public Message{

public:
    RayMsg(RayList &outList, int src, int dst, int chunk, bool idle, int tag);
    RayMsg(RayList* List, int src, int dst, int chunk_size, bool idle, int tag);
};

class RecvMsg : public Message {

public:
    RecvMsg(MPI_Status &status, MPI_Comm comm);
    RecvMsg(RayList* List, RayStreamList* inList, int local_chunk, MPI_Status &status, MPI_Comm comm); 
};

class QuitMsg : public Message {

public:
    QuitMsg(int src, int tag);
};

class StatusMsg : public Message {

public:
    StatusMsg(int , int, int* , int, int, int);
};

class ScheduleMsg : public Message {

public:
    ScheduleMsg(int , int*, int, int, int);
};

Message::Message() {
    printf("construct message void\n");
    content.reserve(1024 * 1024 * 21);
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
RecvMsg::RecvMsg(RayList* List, RayStreamList *inList, int local_chunk, MPI_Status &status, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
//    std::ofstream os; 
//    os.open("out/proc_buffer_" + std::to_string(rank), std::ios::out | std::ios::app ); 
//    os<<"master read from message\n"; 
//    printf("construct message recv ray from proc %d\n", status.MPI_SOURCE);
    
    int count;
    MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &count);
    MPI_Status s0;
    header = new MessageHeader();
    if (count == sizeof(*header)) {
        MPI_Recv((char *)header, count, MPI_UNSIGNED_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm, &s0);
//        if(header->collective)
//            os<<"recv collective\n";
    } else {
        char *buffer = (char *)malloc(count);
        MPI_Recv(buffer, count, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, comm, &s0);
        memcpy(header, buffer, sizeof(*header));
//        printf("recv msg size count %d \n", count);
//        os<<"recv msg sender "<<header->sender <<" root "<<header->root<<"\n";
        if(header->primary == 0 && header->secondary ==0 ) {
            //status
//            os<<"recv status "<<header->content_size<<"\n";
            content.resize(header->content_size);
            memcpy(content.data(), buffer+sizeof(*header), header->content_size);
//            os<<"recv status over\n";
        } else {
            //rays 
            if(header->chunk != local_chunk) {
//                os<<"recv mixed rays "<< header->sender <<" "<<header->content_size<<"\n";
                if(header->chunk < 0) 
                    error("header chunk < 0 ", header->chunk);
                List[header->chunk].read_from_message((char*)buffer+sizeof(MessageHeader), header->primary, header->secondary);

                //RayList::read_from_message(List, (char*)buffer+sizeof(MessageHeader), header->primary, header->secondary);
            } else {
//                os<<"recv normal rays "<< header->sender <<" "<<get_ray_size()<<"\n";
                inList->lock();
                inList->read_from_message((char*)buffer+sizeof(MessageHeader), header->primary, header->secondary, rank); 
                inList->unlock();
            }
        }
        free(buffer);
    }
}

int Message::serialize(struct RayList* outList) {
    if(outList->empty()) return 0;
  
    RaysArray &primary   = outList->get_primary();
    RaysArray &secondary = outList->get_secondary();
    
    
    int width = primary.get_store_width(); 
    int* ids = (int*)(primary.get_data());
    for(int i = 0; i < 5; i ++) {
        printf(" %d %d$", ids[i*width], ids[i*width + 9]);
    }
    printf("\n");

    


    header->primary = primary.get_size(); 
    header->secondary = secondary.get_size(); 

    printf("serialize size %d %d:", header->primary, header->secondary);
    
    int primary_length = header->primary * primary.get_store_width() * sizeof(float);
    int secondary_length = header->secondary * secondary.get_store_width() * sizeof(float);
    content.resize(primary_length + secondary_length);
    printf("Semd Message size %ld\n", content.size());  

    char* ptr = content.data();
    memcpy(ptr, primary.get_data(), primary_length);
    memcpy(ptr + primary_length, secondary.get_data(), secondary_length);
    
    return primary_length + secondary_length;
}

RayMsg::RayMsg(RayList* List, int src, int dst, int chunk_size, bool idle, int tag) {
    std::ofstream os; 
    os.open("out/proc_buffer_worker", std::ios::out | std::ios::app ); 
    
    header = new MessageHeader(-1, src, dst, MsgType::Ray, false, 0, idle, -1, tag);
    
    os<<"RayMsg :\n";
    for(int i = 0; i < chunk_size; i++) {
        header->primary += List[i].primary_size();
        header->secondary += List[i].secondary_size();
        os<<"list "<<i<<"primary copy to "<<header->primary<<" secondary"<<header->secondary<<"\n";
    }
    int primary_length = header->primary * List[0].primary_store_width() * sizeof(float);
    printf("new RayMsg primary %d secondary %d\n", header->primary, header->secondary);
    int secondary_length = header->secondary * List[0].secondary_store_width() * sizeof(float);
    content.resize(primary_length + secondary_length);
    header->content_size = primary_length + secondary_length; 
    
    char* primary_ptr   = content.data();
    char* secondary_ptr = primary_ptr + primary_length;

    for(int i = 0; i < chunk_size; i++) {
        if(!List[i].empty()) {
            RayList &out = List[i];
            int length = out.primary_size() * out.primary_store_width() * sizeof(float);
            memcpy(primary_ptr,   out.get_primary().get_data(), length);
            
            
            int* i_ptr = (int*)primary_ptr;
            printf("test102 %d src %d dst %d\n", i_ptr[9], src, dst);


            primary_ptr += length;
            length = out.secondary_size() * List[i].secondary_store_width() * sizeof(float);
            memcpy(secondary_ptr, out.get_secondary().get_data(), length);
            
            i_ptr = (int*)secondary_ptr;
            printf("test102 %d src %d dst %d\n", i_ptr[9], src, dst);
            
            secondary_ptr += length;  
            out.clear();
        } 
    }
}

RayMsg::RayMsg(RayList &outList, int src, int dst, int chunk, bool idle, int tag) {
    printf("construct message ray\n");
    header = new MessageHeader(-1, src, dst, MsgType::Ray, false, 0, idle, chunk, tag);
    
    RaysArray &primary = outList.get_primary(); 
    int width = primary.get_store_width(); 
    int* ids = (int*)(primary.get_data());
    printf(" ray msg raylist primary size %d  ", primary.get_size());
    for(int i = 0; i < std::min(5, primary.get_size()); i ++) {
        printf(" %d %d$", ids[i*width], ids[i*width + 9]);
    }
    printf("\n");

    RaysArray &secondary = outList.get_secondary(); 
    printf(" ray msg raylist secondary size %d  ", secondary.get_size());
    width = secondary.get_store_width(); 
    ids = (int*)(secondary.get_data());
    for(int i = 0; i < std::min(5, secondary.get_size()); i ++) {
        printf(" %d %d$", ids[i*width], ids[i*width + 9]);
    }
    printf("\n");


    header->primary = primary.get_size(); 
    header->secondary = secondary.get_size(); 
    
    int primary_length = header->primary * primary.get_store_width() * sizeof(float);
    int secondary_length = header->secondary * secondary.get_store_width() * sizeof(float);
    content.resize(primary_length + secondary_length);
    printf("Semd Message size %ld  dst %d primary %d secondary %d\n", content.size(), dst, header->primary, header->secondary);  

    char* ptr = content.data();
    memcpy(ptr, outList.get_primary().get_data(), primary_length);
    memcpy(ptr + primary_length, outList.get_secondary().get_data(), secondary_length);
    
    header->content_size = primary_length + secondary_length; 
    
    outList.clear();
//    printf("outlist clear src %d dst %d\n", src, dst);
}

QuitMsg::QuitMsg(int src, int tag) {
    printf("construct message quit \n");
    header = new MessageHeader(src, -1, -1, MsgType::Quit, true, 0, true, -1, tag);
    printf("new Message root %d collective %d\n", header->root, header->collective);
}

StatusMsg::StatusMsg(int src, int dst, int* status, int chunk, int proc_size, int tag) {
    printf("construct message status src %d\n", src);
    if(dst == -1) 
        header = new MessageHeader(src, -1, dst, MsgType::Status, false, 0, true, chunk, tag);
    else 
        header = new MessageHeader(-1, src, dst, MsgType::Status, false, 0, true, chunk, tag);

   
    int length = proc_size * proc_size * sizeof(int);
    content.resize(length);
    memcpy(content.data(), status, length);
    header->content_size = length;
    
    printf("new Message root %d collective %d\n", header->root, header->collective);
}

//broad cast schedule
ScheduleMsg::ScheduleMsg(int src, int* chunkStatus, int chunk, int chunk_size, int tag) {
    printf(" construct message schedule src %d\n", src);
    header = new MessageHeader(src, -1, -1, MsgType::Schedule, false, 0, false, chunk, tag);
    int length = chunk_size * sizeof(int);
    content.resize(length);
    memcpy(content.data(), chunkStatus, length);
    header->content_size = length;
    printf("new Message root %d collective %d\n", header->root, header->collective);
}
