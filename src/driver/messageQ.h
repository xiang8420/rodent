#include <mutex>
#include <condition_variable>
#include <deque>
#include "message.h"

class MessageQ {
public:
    
    MessageQ(const char *n) : name(n) {
        running = true;
    }

    //! destructor
    ~MessageQ() {
        running = false;
        cond.notify_all(); 
    }

    void kill();

    void enqueue(Message *w); 
    
    Message *dequeue();
    
    int is_ready(); 

    int size() { return workq.size(); }

  //! print the messages pending on this queue
	void printContents(); 

private:
  const char *name;

  std::mutex mutex;
  std::condition_variable cond; // primary, secondary buffer size < max
  bool running;

  std::deque<Message *> workq;
};

