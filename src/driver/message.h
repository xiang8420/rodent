#include <vector>
#include <mutex>
#include <condition_variable>
#include "mpi.h"

struct MessageHeader {
    int  broadcast_root; // will be -1 for point-to-point
    int  sender;                 // will be -1 for broadcast
    int  type;
    bool collective;
    int  content_size;
    int  in_size[3];       //logic size of 3 types of rays primary, secondary, camera ray
    bool idle;      //if this process is idle
    int  chunk;

    MessageHeader(){}

    MessageHeader(int a, int b, int c, bool d, int e, int *in, bool f) 
        : broadcast_root(a), sender(b), type(c), collective(d), content_size(e), idle(f){
       in_size[0] = in[0]; 
       in_size[1] = in[1]; 
    }

    MessageHeader(int a, int b, int c, bool d, int e, bool f, int g) 
        : broadcast_root(a), sender(b), type(c), collective(d), content_size(e), idle(f), chunk(g){}
    
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
    
    Message(MPI_Status &status, MPI_Comm comm);
    
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

	bool is_busy() { return !header->idle; }

    int  get_type() { return header->type; }
    
    void set_type(int t) { header->type = t; }

    int get_chunk() { 
        if (header->chunk < 0)
            printf("header chunk == -1\n");
        return header->chunk;
    }

    void set_chunk(int t) { header->chunk = t; }
    
    char *get_content() { return content.data(); }
    
    size_t get_size() { return content.size(); }

    int get_root() { return header->broadcast_root; }
    
    int get_sender() { return header->sender; }

    bool is_broadcast() { return header->broadcast_root != -1; }
    
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
    
    void deserialize(struct RayList* outList); 
    
    int serialize(struct RayList* outList); 
};

class RayMsg : public Message{

public:
    RayMsg(RayList* outList, int src, int dst, int chunk, bool idle);
};

class BroadcastMsg : public Message {

public:
    BroadcastMsg(int src, int idle);
};

class RecvMsg : public Message {

public:
    
};

class CollectiveMsg : public Message {

public:
    CollectiveMsg(int src);
};

