#include "interface.h"
#include <thread>
#include <mutex>
#include "rayqueue.h"
#include "communicator.h"
#include "communicator_f.h"
#include "client.h"
#include "server.h"

static std::unique_ptr<Server> server;

void setup_server(struct Communicator *comm) {
    server.reset(new Server(comm));
}

void server_run(){
    server->run();
}

void server_send_rays(){
}

bool server_recv_rays(){
}


static std::unique_ptr<Client> client;

void setup_client(struct Communicator *comm, bool server) {
    client.reset(new Client(comm, server));
}

void client_run(struct Settings settings, int dev_num){
    client->run(settings, dev_num);
}

void client_send_rays(float *rays, size_t size, size_t capacity, bool send_all){
    client->send_rays(rays, size, capacity, send_all);
}

int  client_recv_rays(float *rays, size_t size){
    return client->recv_rays(rays, size);
}
