#include "../driver/interface.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <string.h>
#include <sys/file.h>
#include <dirent.h>
#include <unistd.h>

#include "RayList.h"
#include "communicator.h"
#include "decomposition.h"
#include "ProcStatus.h"

#include "Node.h"
#include "P2PNode.h"
#include "MasterWorker.h"

#define PRIMARY_WIDTH 21
#define SECONDARY_WIDTH 14

static void save_image(float *result, const std::string& out_file, size_t width, size_t height, uint32_t iter) {
    ImageRgba32 img;
    img.width = width;
    img.height = height;
    img.pixels.reset(new uint8_t[width * height * 4]);

    auto inv_iter = 1.0f / iter;
    auto inv_gamma = 1.0f / 2.2f;
    for (size_t y = 0; y < height; ++y) {
        for (size_t x = 0; x < width; ++x) {
            auto r = result[(y * width + x) * 3 + 0];
            auto g = result[(y * width + x) * 3 + 1];
            auto b = result[(y * width + x) * 3 + 2];

            img.pixels[4 * (y * width + x) + 0] = clamp(std::pow(r * inv_iter, inv_gamma), 0.0f, 1.0f) * 255.0f;
            img.pixels[4 * (y * width + x) + 1] = clamp(std::pow(g * inv_iter, inv_gamma), 0.0f, 1.0f) * 255.0f;
            img.pixels[4 * (y * width + x) + 2] = clamp(std::pow(b * inv_iter, inv_gamma), 0.0f, 1.0f) * 255.0f;
            img.pixels[4 * (y * width + x) + 3] = 255;
        }
    }

    if (!save_png(out_file, img))
        error("Failed to save PNG file '", out_file, "'");
}
//

struct DistributedFrameWork {
    Node *node;
    Communicator *comm;  
    ProcStatus *ps;

    std::string type;
    
    DistributedFrameWork(std::string type, int chunk, int dev):  type(type) 
    {
        comm = new Communicator();
        ps = new ProcStatus(comm->rank, comm->size, chunk, dev);
        // mpi
        if(type == "P2PNode") {
            node = new P2PNode(comm, ps);
        } else if(type == "MWNode") {
            node = new MWNode(comm, ps);
        } else {
            std::cerr << "Undifined node type\n";
        } 
    }

    ~DistributedFrameWork() {
        printf("delete distributed frame work\n");
        delete node;
        delete comm;
        delete ps;
    }
    
    void run(ImageDecomposition *camera) {
        printf("dis frame worker run\n");
        int worker_size = comm->get_comm_size();
        printf("before all gather %d\n", comm->rank);
        
        camera->decomposition(ps->get_chunk_map(), true, comm->rank, comm->size); 
        ps->updata_local_chunk();

        node->run(camera);
        
    }

    void gather_image(const std::string& out_file, float* film, int frame, int width, int height) {
        //film 
        printf("dfw reduce %d %d %d %d\n", comm->rank, comm->master, width, height); 
        int pixel_num = width * height * 3;
        float *reduce_buffer = new float[pixel_num];
        printf("%d before reduce\n", comm->rank);
        
        comm->reduce_image(film, reduce_buffer, pixel_num);
        
        std::string out = out_file + "f" + std::to_string(frame) + "w" + std::to_string(comm->get_comm_size()) + "g" + std::to_string(ps->get_chunk_size());
        out += ".png";
        if (comm->rank == 0 && out_file != "") {
            save_image(reduce_buffer, out, width, height, 1 /* iter*/ );
        }
        printf("%d end\n", comm->rank);  
    }
};

static std::unique_ptr<DistributedFrameWork> dfw;

void setup_distributed_framework(std::string type, int chunk, int dev) {
    dfw.reset(new DistributedFrameWork(type, chunk, dev));
}

void cleanup_distributed_framework() {
    dfw.reset();
}

void dfw_run(ImageDecomposition *camera) {
    dfw->run(camera);
}

void dfw_save_image(const std::string& out_file, float* film, int frame, int width, int height) {
    dfw->gather_image(out_file, film, frame, width, height);
}

void send_rays(float *rays, size_t size, size_t capacity, bool isPrimary){
    dfw->node->save_outgoing_buffer(rays, size, capacity, isPrimary);
}

int recv_rays(float **rays, size_t size, bool isPrimary, int thread_id, bool thread_wait){
    return dfw->node->load_incoming_buffer(rays, size, isPrimary, thread_id, thread_wait);
}

void master_save_ray_batches(float *rays, size_t size, size_t capacity, size_t thread_id) {
//    dfw->worker->save_ray_batches(rays, size, capacity, thread_id);
}

int32_t thread_buffer_size() {
    return dfw->node->proc_status()->get_buffer_size();
}
