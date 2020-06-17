#include "messageQ.h"

void MessageQ::enqueue(Message *w) {
    std::unique_lock <std::mutex> lock(mutex); 
    workq.push_back(w);
    cond.notify_all();
    mutex.unlock();
}


Message *MessageQ::dequeue() {
    
    std::unique_lock <std::mutex> lock(mutex); 
    while (workq.empty() && running)
        cond.wait(lock); 

    Message *r = NULL;
    if (! workq.empty()) {
        r = workq.front();
        workq.pop_front();
    }
    return r;

}

int MessageQ::is_ready() {
    mutex.lock();
    int t = (workq.empty() && running) ? 0 : 1;
    mutex.unlock();
    return t;
}

void MessageQ::printContents() {
   // for(auto a = workq.begin(); a != workq.end(); ++a)
   //     GetTheApplication()->Identify(*a);
}

void MessageQ::kill() {
    std::unique_lock <std::mutex> lock(mutex); 
    running = false;
	cond.notify_all();
}
	
