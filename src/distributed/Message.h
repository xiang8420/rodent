#include <vector>
#include <mutex>
#include "RayArrayList.h"
#include "RayStreamList.h"
#include <condition_variable>
#include "mpi.h"

enum MsgType {Status, ArrayRay, StreamRay, Schedule, Quit};

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
    char *content;

public:
    Message(){};
    
    Message(MessageHeader *h): header(h) { }

    ~Message() { delete header; }

    // write content to a raylists
    int write_raylist(struct RayArrayList *outlist);
  
    // read inlist to message 
    int read_raylist(struct RayArrayList *inlist);

    int get_tag() {return header->tag;}
	bool is_busy() { return !header->idle; }
    int  get_type() { return header->type; }
    void set_type(int t) { header->type = MsgType(t); }

    int get_chunk() { 
        if (header->chunk < 0)
            printf("header chunk == -1\n");
        return header->chunk;
    }

    void set_chunk(int t) { header->chunk = t; }
    
    char *get_content() { return content; }
    
    size_t get_ray_size() {return header->primary + header->secondary;}

    size_t primary_size(){return header->primary;}
    size_t secondary_size(){return header->secondary;}

    size_t get_content_size() {return header->content_size; }

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
    
    void deserialize(char* ptr, struct RayArrayList* outList); 
    
    int serialize(struct RayArrayList* outList) {
        if(outList->empty()) return 0;
      
        RaysArray &primary   = outList->get_primary();
        RaysArray &secondary = outList->get_secondary();
        
        
        int width = primary.get_store_width(); 
        int* ids = (int*)(primary.get_data());

        header->primary = primary.get_size(); 
        header->secondary = secondary.get_size(); 
        
        int primary_length = header->primary * primary.get_store_width() * sizeof(float);
        int secondary_length = header->secondary * secondary.get_store_width() * sizeof(float);
        content = new char[primary_length + secondary_length];

        char* ptr = content;
        memcpy(ptr, primary.get_data(), primary_length);
        memcpy(ptr + primary_length, secondary.get_data(), secondary_length);
        
        return primary_length + secondary_length;
    }
    
};

class RayMsg : public Message{

public:
    RayMsg(RayArrayList* List, int src, int dst, int chunk_size, bool idle, int tag) {
        std::ofstream os; 
        os.open("out/proc_buffer_worker", std::ios::out | std::ios::app ); 
        
        header = new MessageHeader(-1, src, dst, MsgType::ArrayRay, false, 0, idle, -1, tag);
        
        os<<"RayMsg :\n";
        for(int i = 0; i < chunk_size; i++) {
            header->primary += List[i].primary_size();
            header->secondary += List[i].secondary_size();
            os<<"list "<<i<<"primary copy to "<<header->primary<<" secondary"<<header->secondary<<"\n";
        }
        int primary_length = header->primary * List[0].primary_store_width() * sizeof(float);
        printf("new RayMsg primary %d secondary %d\n", header->primary, header->secondary);
        int secondary_length = header->secondary * List[0].secondary_store_width() * sizeof(float);
        content = new char[primary_length + secondary_length];
        header->content_size = primary_length + secondary_length; 
        
        char* primary_ptr   = content;
        char* secondary_ptr = primary_ptr + primary_length;

        for(int i = 0; i < chunk_size; i++) {
            if(!List[i].empty()) {
                RayArrayList &out = List[i];
                int length = out.primary_size() * out.primary_store_width() * sizeof(float);
                memcpy(primary_ptr,   out.get_primary().get_data(), length);
                
                primary_ptr += length;
                length = out.secondary_size() * List[i].secondary_store_width() * sizeof(float);
                memcpy(secondary_ptr, out.get_secondary().get_data(), length);
                
                secondary_ptr += length;  
                out.clear();
            } 
        }
    }
    RayMsg(RayArrayList &outList, int src, int dst, int chunk, bool idle, int tag) {
        printf("construct message ray\n");
        header = new MessageHeader(-1, src, dst, MsgType::ArrayRay, false, 0, idle, chunk, tag);
        
        RaysArray &primary = outList.get_primary(); 
        int width = primary.get_store_width(); 
        int* ids = (int*)(primary.get_data());
        printf(" dst %d ray msg raylist primary size %d width %d %d : ", dst, primary.get_size(), width, primary.get_store_width());
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
        content = new char[primary_length + secondary_length];
        printf("Semd Message size %ld  dst %d primary %d secondary %d\n",primary_length + secondary_length, dst, header->primary, header->secondary);  

        char* ptr = content;
        memcpy(ptr, outList.get_primary().get_data(), primary_length);
        memcpy(ptr + primary_length, outList.get_secondary().get_data(), secondary_length);
        
        header->content_size = primary_length + secondary_length; 
        
        outList.clear();
    //    printf("outlist clear src %d dst %d\n", src, dst);
    }

    RayMsg(RayStreamList &outList, int src, int dst, int chunk, bool idle, int tag) {
        printf("construct message ray\n");
        header = new MessageHeader(-1, src, dst, MsgType::StreamRay, false, 0, idle, chunk, tag);

        header->primary = outList.primary_size(); 
        header->secondary = outList.secondary_size(); 
       
        // first element is size 
        int primary_length = header->primary * (1 + 21 * outList.get_store_capacity()) * sizeof(float);
        int secondary_length = header->secondary * (1 + 14 * outList.get_store_capacity()) * sizeof(float);
        content = new char[primary_length + secondary_length];
        header->content_size = primary_length + secondary_length; 
       
        char* ptr = content;
        printf("m %d > %d send RayMsg primary %d secondary %d\n", src, dst, header->primary, header->secondary);
        int primary_size = outList.primary_size();
        for(int i = 0; i < primary_size; i++) {
            int * iptr = (int*) ptr;
            struct RaysStream *primary = outList.get_primary();
            iptr[0] = primary->get_size();
            printf("primary size %d :", iptr[0]);
            ptr += sizeof(int);
            int length = 21 * outList.get_store_capacity() * sizeof(float);
            memcpy(ptr, primary->get_data(), length);
            int * irayptr = (int*) primary->get_data();
            printf("size %d ", iptr[0]);
            if(iptr[0] > 6) {
                //for(int i = iptr[0] - 5; i < iptr[0]; i++) {
                for(int i = 0; i < 5; i++) {
                    printf("%d %d|%d|", irayptr[i], irayptr[i + 9 * outList.get_store_capacity()], src);
                }
                printf("write to msg %d : ", iptr[0]);
                //for(int i = iptr[0] - 5; i < iptr[0]; i++) {
                for(int i = 0; i < 5; i++) {
                    printf("%d %d|%d|", iptr[1 + i], iptr[1 + i + 9 * outList.get_store_capacity()], src);
                }
                printf("\n");
            }
            delete primary;
            ptr += length;
        } 
            printf("end primary size\n");
        int secondary_size = outList.secondary_size();
        for(int i = 0; i < secondary_size; i++) {
            int * iptr = (int*) ptr;
            struct RaysStream *secondary = outList.get_secondary();
            iptr[0] = secondary->get_size();
            printf("primary size %d :", iptr[0]);
            ptr += sizeof(int);
            int length = 14 * outList.get_store_capacity() * sizeof(float);
            memcpy(ptr, secondary->get_data(), length);
            int * irayptr = (int*) secondary->get_data();
            if(iptr[0] > 6) {
                //for(int i = iptr[0] - 5; i < iptr[0]; i++) {
                for(int i = 0; i < 5; i++) {
                    printf("%d %d||", irayptr[i], irayptr[i + 9 * outList.get_store_capacity()]);
                }
                printf("write to msg %d : ", iptr[0]);
                //for(int i = iptr[0] - 5; i < iptr[0]; i++) {
                for(int i = 0; i < 5; i++) {
                    printf("%d %d||", iptr[1 + i], iptr[1 + i + 9 * outList.get_store_capacity()]);
                }
                printf("\n");
               // for(int i = iptr[0]-6; i < iptr[0]; i++) {
               //     printf("%d %d||", irayptr[i], irayptr[i + 9 * outList.get_store_capacity()]);
               // }
               // printf("write to msg %d : ", iptr[0]);
               // for(int i = iptr[0]-6; i < iptr[0]; i++) {
               //     printf("%d %d||", iptr[1 + i], iptr[1 + i + 9 * outList.get_store_capacity()]);
               // }
               // printf("\n");
            }
            delete secondary;
            ptr += length;
        } 
    }
    ~RayMsg() { delete[] content; }
};

class RecvMsg : public Message {

public:
    //recv msg
    RecvMsg(MPI_Status &status, MPI_Comm comm) {
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
            content = new char[header->content_size];
            memcpy(content, buffer+sizeof(*header), header->content_size);
            free(buffer);
        }
    }

    //recv message, if it's ray msg load it to list 
    RecvMsg(RayArrayList* outArray, RayStreamList * outStream, RayStreamList *inList, int local_chunk, MPI_Status &status, MPI_Comm comm) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        
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
            if(header->primary == 0 && header->secondary ==0 ) {
                content = new char[header->content_size];
                memcpy(content, buffer+sizeof(*header), header->content_size);
            } else {
                int cur_chk = header->chunk;
                if(cur_chk != local_chunk) {
                    if(header->type == MsgType::ArrayRay) {
                        if(cur_chk < 0) 
                            error("header chunk < 0 ", cur_chk);
                        
                        std::unique_lock <std::mutex> lock(outArray[0].mutex); 
                        outArray[cur_chk].get_primary().read_from_ptr((char*)buffer+sizeof(MessageHeader), header->primary);
                        char* secondary_ptr = (char*)buffer + sizeof(MessageHeader) + header->primary * 21 * sizeof(float); 
                        outArray[cur_chk].get_secondary().read_from_ptr(secondary_ptr, header->secondary);
                    } else {
                        std::unique_lock <std::mutex> lock(outStream[0].mutex); 
                        outStream[cur_chk].read_from_stream_message((char*)buffer+sizeof(MessageHeader), header->primary, header->secondary, rank); 
                    }
                } else {
                    statistics.start("run => message_thread => recv_message => RecvMsg => read_from_message");
                    
                    std::unique_lock <std::mutex> lock(inList->mutex); 
                    while (inList->full()) {
                        inList->cond_empty.wait(lock);
                    }
                    if(header->type == MsgType::ArrayRay) {
                        inList->read_from_array_message((char*)buffer+sizeof(MessageHeader), header->primary, header->secondary, rank); 
                    } else {
                        inList->read_from_stream_message((char*)buffer+sizeof(MessageHeader), header->primary, header->secondary, rank); 
                    }
                    os<<"inList size "<< inList->size()<<"\n";
                    lock.unlock();
                    statistics.end("run => message_thread => recv_message => RecvMsg => read_from_message");
                }
            }
            free(buffer);
        }
    }

    ~RecvMsg() { delete[] content; }
};

class QuitMsg : public Message {

public:
    QuitMsg(int src, int tag) {
        header = new MessageHeader(src, -1, -1, MsgType::Quit, true, 0, true, -1, tag);
    }
};

class StatusMsg : public Message {

public:
    StatusMsg(int src, int dst, int* status, int chunk, int proc_size, int tag) {
        printf("construct message status src %d\n", src);
        if(dst == -1) 
            header = new MessageHeader(src, -1, dst, MsgType::Status, false, 0, true, chunk, tag);
        else 
            header = new MessageHeader(-1, src, dst, MsgType::Status, false, 0, true, chunk, tag);
       
        int length = proc_size * proc_size * sizeof(int);
        content = new char[length];
        memcpy(content, status, length);
        header->content_size = length;
        
        printf("new Message root %d collective %d\n", header->root, header->collective);
    }
};

class ScheduleMsg : public Message {

public:
    //broad cast schedule
    ScheduleMsg(int src, int* chunkStatus, int chunk, int chunk_size, int tag) {
        printf(" construct message schedule src %d\n", src);
        header = new MessageHeader(src, -1, -1, MsgType::Schedule, false, 0, false, chunk, tag);
        int length = chunk_size * sizeof(int);
        content = new char[length];
        memcpy(content, chunkStatus, length);
        header->content_size = length;
        printf("new Message root %d collective %d\n", header->root, header->collective);
    }
    ~ScheduleMsg(){
        delete[] content;
    }
};




