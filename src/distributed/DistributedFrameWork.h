#include "../driver/interface.h"
#include <thread>
#include <mutex>
#include <condition_variable>
#include <string.h>
#include <sys/file.h>
#include <dirent.h>
#include <unistd.h>
#include <assert.h>
#include <cstdlib>

extern "C" {
void rodent_present(int32_t dev);
void rodent_unload_chunk_data(int32_t chk, int32_t dev); 
int* rodent_get_chunk_hit();
void rodent_update_render_chunk_hit(int32_t*, int32_t);
int* rodent_get_render_chunk_hit_data(); 
}
#include "../driver/common.h"
#include "statistic.h"
#include "Message.h"
#include "scheduler.h"
#include "ProcStatus.h"
#include "communicator.h"

#include "Node.h"
#include "SingleNode.h"
#include "SyncNode.h"
#include "AllCopyNode.h"
//#include "MasterWorker.h"
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

struct DistributedFrameWork {
    Scheduler *scheduler;
    Node *node;
    Communicator *comm;
    ProcStatus *ps;

    std::string type;

    DistributedFrameWork(std::string dis_type, int width, int height, int spp):  type(dis_type) 
    {
        int chunk_size = get_chunk_num();
        int dev = get_dev_num();

        comm          = new Communicator();
        scheduler     = new Scheduler(width, height, spp, comm->get_rank(), comm->get_size());
        ps            = new ProcStatus(comm->get_rank(), comm->get_size(), chunk_size, dev);
        ps->chunk_manager = scheduler->chunk_manager;


        // mpi
        if(chunk_size == 1 && comm->get_size() == 1 ) 
            type = "Single";
        std::cout<<"new type "<<type<<" "<<chunk_size<<" "<<comm->get_size()<<"\n";

        if(type == "Single") node = new SingleNode(comm, ps, scheduler);
        else if(type == "Sync") node = new SyncNode(comm, ps, scheduler);
//        else if(type == "MasterWorker") node = new MWNode(comm, ps);
        else if(type == "Async") node = new AsyncNode(comm, ps, scheduler);
        else if(type == "AllCopy") node = new AllCopyNode(comm, ps, scheduler); 
        else error("Unknown node type");
    }

    ~DistributedFrameWork() {
        printf("delete distributed frame work\n");
        delete node;
        delete ps;
        delete comm;
        delete scheduler;
    }
    
    void run(Camera *camera) {
        printf("dis frame worker run\n");
        int proc_rank = comm->get_rank();
        int proc_size = comm->get_size(); 
        /*block size equels proc size*/
        statistics.start("schedule");
        std::cout<<"run type "<<type<<"\n";
        if( type == "AllCopy" || type == "Single") {
            int block_count = comm->get_size() == 1 ? 1 : comm->get_size() * 2;
            scheduler->split_image_block(camera, block_count, comm->get_rank(), comm->get_size());
        } else {
            int block = comm->get_size();
            scheduler->image_domain_decomposition(camera, block, proc_rank, proc_size, ps->get_simple_trace(), type == "Sync");
        }
        scheduler->generate_chunk_manager();
        statistics.end("schedule");

        statistics.start("run");
        node->run();
        
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
        delete[] reduce_buffer;
    }

};

static std::unique_ptr<DistributedFrameWork> dfw;

void setup_distributed_framework(std::string type, int w, int h, int spp) {
    dfw.reset(new DistributedFrameWork(type, w, h, spp));
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
    printf("dfw save rays\n");
    dfw->node->save_outgoing_buffer(rays, size, isPrimary);
    statistics.end("run => work_thread => send");
}

int recv_rays(float **rays, bool isPrimary, int thread_id){
    statistics.start("run => work_thread => load_incoming_buffer");
    int res = dfw->node->load_incoming_buffer(rays, isPrimary, thread_id);
    statistics.end("run => work_thread => load_incoming_buffer");
    return res;
}

int32_t dfw_chunk_size() { return dfw->ps->get_chunk_size(); };

int32_t dfw_stream_logic_capacity() { return dfw->ps->get_stream_logic_capacity(); }

int32_t dfw_stream_store_capacity() { return dfw->ps->get_stream_store_capacity(); }

int32_t dfw_out_stream_capacity() { return dfw->ps->get_out_stream_capacity(); }

int32_t dfw_thread_num() { return dfw->ps->get_thread_num(); }

int32_t dfw_mpi_rank() { return dfw->comm->get_rank(); }

void dfw_print_memory(int i) { return dfw->node->loop_check(i); }

void dfw_time_start(std::string s) { statistics.start(s); }

void dfw_time_end(std::string s) { statistics.start(s); }

