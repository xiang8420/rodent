#include "interface.h"
#include <thread>
#include <mutex>
#include <string.h>
#include "rayqueue.h"
#include "communicator.h"
#include <sys/file.h>
#include <dirent.h>
#include <unistd.h>

#define PRIMARY_WIDTH 21
#define SECONDARY_WIDTH 14
#define MAX_CLIENT_CHUNK 3

#include "worker.h"
#include "master.h"


static std::unique_ptr<Master> master;

void setup_master(struct Communicator *comm, int chunk_size) {
    master.reset(new Master(comm, chunk_size));
}

void master_run(){
    master->run();
}

static std::unique_ptr<Worker> worker;

void worker_run(struct Settings settings, int dev_num){
    worker->run(settings, dev_num);
}

void setup_worker(struct Communicator *comm, bool master) {
    worker.reset(new Worker(comm, master));
}

void worker_send_rays(float *rays, size_t size, size_t capacity, bool primary, bool send_all){
    worker->Write_outcome_buffer(rays, size, capacity, primary, send_all);
}

int  worker_recv_rays(float **rays, size_t size, bool idle, bool primary){
    return worker->Read_income_buffer(rays, size, idle, primary);
}
