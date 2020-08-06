#include <vector>
#include <mutex>
#include "RayList.h"
#include <condition_variable>
#include "mpi.h"

enum MsgType {Status, Ray, Schedule, Quit};

struct MessageHeader {
    int      root;                   // will be -1 for point-to-point
    int      sender;                 // will be -1 for broadcast
    MsgType  type; //0 status, 1 rays , 2 quit
    int      content_size;
    int      primary, secondary;    
    int      chunk;   //-1 means mixed rays
    bool     collective;
    bool     idle;      //if this process is idle
//    int      padding[2]; 

    MessageHeader(){}

    MessageHeader(int root, int sender, MsgType type, bool collective, int size, bool idle, int chunk) 
        : root(root), sender(sender), type(type), collective(collective), content_size(size), idle(idle), chunk(chunk)
    {
       primary   = 0; 
       secondary = 0; 
    }
    
    bool HasContent() { return content_size > 0; }
    
};

class Message { 
protected:

    int destination;
    MessageHeader *header;
    std::vector<char> content;
    std::mutex mutex;
    std::condition_variable cond; 

public:
    Message();
    
    Message(MessageHeader *h): header(h) { }

    // write content to a raylists
    int write_raylist(struct RayList *outlist);
  
    // read inlist to message 
    int read_raylist(struct RayList *inlist);
    
    //! notify to a blocking thread that it should resume processing
	void notify() { 
		mutex.lock();
		header->idle = true;
        cond.notify_one();
	    mutex.unlock();
    }

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

    void set_destination(int i) { destination = i; }
    
    int get_destination() { return destination; }

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
    RayMsg(RayList* outList, int src, int dst, int chunk, bool idle);
    RayMsg(RayList** List, int src, int dst, int chunk_size, bool idle);
};

class RecvMsg : public Message {

public:
    RecvMsg(MPI_Status &status, MPI_Comm comm);
    
    RecvMsg(RayList** List, RayStreamList* inList, int local_chunk, MPI_Status &status, MPI_Comm comm); 
};

class QuitMsg : public Message {

public:
    QuitMsg(int src);
};

class StatusMsg : public Message {

public:
    StatusMsg(int , int, int* , int, int);
};

class ScheduleMsg : public Message {

public:
    ScheduleMsg(int , int*, int, int);
};
