#include "../driver/interface.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <string.h>
#include <sys/file.h>
#include <dirent.h>
#include <unistd.h>
#include <assert.h>

extern "C" {
void rodent_present(int32_t dev);
void rodent_unload_chunk_data(int32_t dev ); 
}

#include "MemoryPool.h"
#include "statistic.h"
#include "ProcStatus.h"
#include "communicator.h"
#include "decomposition.h"

#include "Node.h"
#include "SingleNode.h"
#include "SyncNode.h"
#include "AllCopyNode.h"
#include "MasterWorker.h"

#define PRIMARY_WIDTH 21
#define SECONDARY_WIDTH 14

#include "AsyncNode.h"

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
    ImageDecomposition *splitter;
    Node *node;
    Communicator *comm;
    ProcStatus *ps;
    bool rough_trace;

    std::string type;

    DistributedFrameWork(std::string dis_type, int chunk, int dev, int width, int height, int spp):  type(dis_type) 
    {
        splitter = new ImageDecomposition(width, height, spp);
        comm = new Communicator();
        rough_trace = false; //true || (comm->get_size() % 2 == 1 && comm->get_size() > 1);
        ps = new ProcStatus(comm->get_rank(), comm->get_size(), chunk, dev, rough_trace);
        // mpi
        if(chunk == 1 && comm->get_size() == 1 ) 
            type = "Single";
        std::cout<<"new type "<<type<<" "<<chunk<<" "<<comm->get_size()<<"\n";

        if(type == "Single") node = new SingleNode(comm, ps);
        else if(type == "Sync") node = new SyncNode(comm, ps);
        else if(type == "MasterWorker") node = new MWNode(comm, ps);
        else if(type == "Async") node = new AsyncNode(comm, ps);
        else if(type == "AllCopy") node = new AllCopyNode(comm, ps); 
        else error("Unknown node type");
        
    }

    ~DistributedFrameWork() {
        printf("delete distributed frame work\n");
        delete node;
        delete ps;
        delete comm;
    }
    
    void run(Camera *camera) {
        statistics.start("run");

        printf("dis frame worker run\n");
        int proc_rank = comm->get_rank();
        int proc_size = comm->get_size(); 
        /*block size equels proc size*/
        if( type == "AllCopy" || type == "Single") {
            int block_count = comm->get_size() == 1 ? 1 : comm->get_size() * 2;
            splitter->split_image_block(camera, block_count, comm->get_rank(), comm->get_size());
        } else {
            int block_count = comm->get_size();
            splitter->image_domain_decomposition(camera, block_count, proc_rank, proc_size, rough_trace);
        }
        splitter->write_chunk_proc(ps->get_chunk_proc());
        ps->updata_local_chunk();
        node->run(splitter);
        
        statistics.end("run");
        statistics.print(comm->os);
        ps->reset();
    }

    void gather_image(const std::string& out_file, float* film, int fid, int frame, int width, int height) {

        Communicator * comm = node->get_communicator();
        ProcStatus * ps = node->get_proc_status();
        
        for(int i = 0; i < ps->get_dev_num(); i++) 
            rodent_present(i);

        //film 
        printf("dfw reduce %d %d %d %d\n", comm->get_rank(), comm->get_master(), width, height); 
        int pixel_num = width * height * 3;
        float *reduce_buffer = new float[pixel_num];
        printf("%d before reduce\n", comm->get_rank());
        
        comm->reduce(film, reduce_buffer, pixel_num);
        
        std::string out = "picture/" + out_file + "_f_" + std::to_string(fid) + "_w_" + std::to_string(comm->get_size()) + "_g_" + std::to_string(ps->get_chunk_size());

        save_image(film, out + "_rank_" + std::to_string(comm->get_rank()) + ".png", width, height, 1 /* iter*/ );
        printf("%d end\n", comm->get_rank());  

        if (comm->get_rank() == 0 && out_file != "") {
            if(fid == frame)
                save_image(reduce_buffer, out + ".png", width, height, frame);
            else 
                save_image(reduce_buffer, out + ".png", width, height, 1);
        }
    }
};

static std::unique_ptr<DistributedFrameWork> dfw;

void setup_distributed_framework(std::string type, int chunk, int dev, int w, int h, int spp) {
    dfw.reset(new DistributedFrameWork(type, chunk, dev, w, h, spp));
}

void cleanup_distributed_framework() {
    dfw.reset();
}

void dfw_run(Camera *camera) {
    dfw->run(camera);
}

void dfw_save_image(const std::string& out_file, float* film, int fid, int frame, int width, int height) {
    dfw->gather_image(out_file, film, fid, frame, width, height);
}

void send_rays(float *rays, size_t size, size_t capacity, bool isPrimary){
    statistics.start("run => work_thread => send");
    dfw->node->save_outgoing_buffer(rays, size, isPrimary);
    statistics.end("run => work_thread => send");
}

int recv_rays(float **rays, size_t size, bool isPrimary, int thread_id, bool thread_wait){
    statistics.start("run => work_thread => load_incoming_buffer");
    int res = dfw->node->load_incoming_buffer(rays, size, isPrimary, thread_id, thread_wait);
    statistics.end("run => work_thread => load_incoming_buffer");
    return res;
}

int32_t dfw_stream_logic_capacity() {
    return dfw->ps->get_stream_logic_capacity();
}

int32_t dfw_stream_store_capacity() {
    return dfw->ps->get_stream_store_capacity();
}

int32_t dfw_out_stream_capacity() {
    return dfw->ps->get_out_stream_capacity();
}

int32_t dfw_thread_num() {
    return dfw->ps->get_thread_num(); 
}


