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
void rodent_unload_chunk_data(int32_t dev ); 
int* rodent_get_light_field();
}
#include "../driver/common.h"
#include "MemoryPool.h"
#include "statistic.h"
#include "ProcStatus.h"
#include "communicator.h"
#include "decomposition.h"

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

static void save_image_its(int* reduce_buffer, int chunk_size, float spp) {
    int res = LIGHT_FIELD_RES;
    ImageRgba32 img;
    img.width = res * 6 + 20;
    img.height = res * chunk_size;
    img.pixels.reset(new uint8_t[img.width * img.height * 4]);
    
    int width = img.width; 
    int height = img.height;
    
    std::vector<rgb> color;
    for(int i = 0; i < chunk_size; i++) {
        rgb col(rand() % 255, rand() % 255, rand() % 255 ); 
        color.emplace_back(col);
        for(int h = 0; h < 128; h++) {
            int ph = i * 128 + h; 
            for(int w = 0; w < 19; w ++) {

                img.pixels[4 * (ph * width + w) + 0] = col.x; 
                img.pixels[4 * (ph * width + w) + 1] = col.y;
                img.pixels[4 * (ph * width + w) + 2] = col.z;
                img.pixels[4 * (ph * width + w) + 3] = 255;
            }
        }
    }


    int all_face_size = res * res * 6;
    int face_size = res * res;
    int inv_spp = 1;
    for(int i = 0; i < chunk_size; i++) {
        int h_st = res * i; 
        for(int j = 0; j < 6; j++) {
            int w_st = 20 + res * j;
            for(int u = 0; u < res; u++) {
                int pu = w_st + u;
                for(int v = 0; v < res; v++) {
                    int pv = h_st + v;
                    
                    if(v == res - 1 || u == res - 1) {
                        img.pixels[4 * (pv * width + pu) + 0] = 0;
                        img.pixels[4 * (pv * width + pu) + 1] = 255;
                        img.pixels[4 * (pv * width + pu) + 2] = 0;
                        img.pixels[4 * (pv * width + pu) + 3] = 255;
                    } else { 
                        int its = (reduce_buffer[all_face_size * i + face_size * j + (res - v) * res + u]) >> 8;
                        if(its >= 255 || its == 0) continue;
                        its = its & 0x7F;
                        rgb& col = color[its];
                        img.pixels[4 * (pv * width + pu) + 0] = col.x;
                        img.pixels[4 * (pv * width + pu) + 1] = col.y;
                        img.pixels[4 * (pv * width + pu) + 2] = col.z;
                        img.pixels[4 * (pv * width + pu) + 3] = 255;
                    }
                }
            }
        } 
    }
    if (!save_png(std::string("picture/light_field_its.png"), img))
        error("Failed to save PNG file light_field.png");
}

static void save_image_ctrb(int* reduce_buffer, int chunk_size, float spp) {
    int res = LIGHT_FIELD_RES;
    ImageRgba32 img;
    img.width = res * 6;
    img.height = res * chunk_size;
    img.pixels.reset(new uint8_t[img.width * img.height * 4]);

    int width = img.width; 
    int height = img.height;
    int all_face_size = res * res * 6;
    int face_size = res * res;
    int inv_spp = 1;//5 / spp;

    for(int i = 0; i < chunk_size; i++) {
        int h_st = res * i; 
        for(int j = 0; j < 6; j++) {
            int w_st = res * j;
            for(int u = 0; u < res; u++) {
                int pu = w_st + u;
                for(int v = 0; v < res; v++) {
                    int pv = h_st + v;
                    
                    int lu = reduce_buffer[all_face_size * i + face_size * j + (res - v) * res + u];
                    int lu_q0 = lu & 0xFF;
                    int lu_q1 = lu >> 8  & 0xFF;
                    int lu_q2 = lu >> 16 & 0xFF;
                    int lu_q3 = lu >> 24 & 0xFF;
                    
                    img.pixels[4 * (pv * width + pu) + 0] = lu_q3; 
                    img.pixels[4 * (pv * width + pu) + 1] = lu_q3;
                    img.pixels[4 * (pv * width + pu) + 2] = lu_q3;
                    img.pixels[4 * (pv * width + pu) + 3] = 255;
                    if(v == res - 1 || u == res - 1) {
                        img.pixels[4 * (pv * width + pu) + 0] = 255;
                        img.pixels[4 * (pv * width + pu) + 1] = 0;
                        img.pixels[4 * (pv * width + pu) + 2] = 0;
                        img.pixels[4 * (pv * width + pu) + 3] = 255;
                    } 
                } 
            }
        } 
    }
    if (!save_png(std::string("picture/light_field_ctrb.png"), img))
        error("Failed to save PNG file light_field.png");
}
//

struct DistributedFrameWork {
    Scheduler *scheduler;
    Node *node;
    Communicator *comm;
    ProcStatus *ps;
    bool rough_trace;

    std::string type;

    DistributedFrameWork(std::string dis_type, int chunk, int dev, int width, int height, int spp):  type(dis_type) 
    {
        scheduler = new Scheduler(width, height, spp);
        comm = new Communicator();
        rough_trace = false; //true || (comm->get_size() % 2 == 1 && comm->get_size() > 1);
        ps = new ProcStatus(comm->get_rank(), comm->get_size(), chunk, dev, rough_trace);
        // mpi
        if(chunk == 1 && comm->get_size() == 1 ) 
            type = "Single";
        std::cout<<"new type "<<type<<" "<<chunk<<" "<<comm->get_size()<<"\n";

        if(type == "Single") node = new SingleNode(comm, ps);
        else if(type == "Sync") node = new SyncNode(comm, ps);
//        else if(type == "MasterWorker") node = new MWNode(comm, ps);
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
    
    void save_light_field(int* light_field, float spp) {
        int size = LIGHT_FIELD_RES * LIGHT_FIELD_RES * 2 * 6 * ps->get_chunk_size();
        int *reduce_buffer = new int[size];
        comm->reduce(light_field, reduce_buffer, size);
        int not_zero = 0;
        for(int i = 0; i < size; i++) {
            if(light_field[i] != 0){
                not_zero ++;
                comm->os<<i<<" "<<light_field[i]<<" | "; 
            } 
        }
        printf("not zero %d\n", not_zero);
        if(comm->get_rank() == 0 && VIS_LIGHT_FIELD) { 
            save_image_its(reduce_buffer + LIGHT_FIELD_RES * LIGHT_FIELD_RES * 6 * ps->get_chunk_size(), ps->get_chunk_size(), spp);
            save_image_ctrb(reduce_buffer, ps->get_chunk_size(), spp);
        }
    }

    void run(Camera *camera) {
        statistics.start("run");

        printf("dis frame worker run\n");
        int proc_rank = comm->get_rank();
        int proc_size = comm->get_size(); 
        /*block size equels proc size*/
        std::cout<<"run type "<<type<<"\n";
        if( type == "AllCopy" || type == "Single") {
            int block_count = comm->get_size() == 1 ? 1 : comm->get_size() * 2;
            scheduler->split_image_block(camera, block_count, comm->get_rank(), comm->get_size());
        } else {
            int block_count = comm->get_size();
            scheduler->image_domain_decomposition(camera, block_count, proc_rank, proc_size, rough_trace, node->sync);
        }
        scheduler->write_chunk_proc(ps->get_chunk_proc());
        ps->updata_local_chunk();
        node->run(scheduler);
        
        statistics.end("run");
        statistics.print(comm->os);
       
        save_light_field(rodent_get_light_field(), scheduler->get_spp());

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
    printf("dfw save rays\n");
    dfw->node->save_outgoing_buffer(rays, size, isPrimary);
    statistics.end("run => work_thread => send");
}

int recv_rays(float **rays, size_t size, bool isPrimary, int thread_id, bool thread_wait){
    statistics.start("run => work_thread => load_incoming_buffer");
    int res = dfw->node->load_incoming_buffer(rays, size, isPrimary, thread_id, thread_wait);
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
