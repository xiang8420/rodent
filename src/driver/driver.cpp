#include <memory>
#include <sstream>
#include <algorithm>
#include <string>
#include <cstring>
#include <chrono>
#include <cmath>
#include <mpi.h>

#ifndef DISABLE_GUI
#include <SDL2/SDL.h>
#endif

#include "interface.h"
#include "float3.h"
#include "common.h"
#include "image.h"
#include "process_status.h"
#include "communicator.h"
#include "multiproc.h"
#if defined(__x86_64__) || defined(__amd64__) || defined(_M_X64)
#include <x86intrin.h>
#endif


void setup_interface(size_t, size_t);
float* get_pixels();
void clear_pixels();
void cleanup_interface();
float* get_first_primary();
void setup_master(struct Communicator *, struct ProcStatus *);
void cleanup_master();
void master_run();

void setup_worker(struct Communicator *, struct ProcStatus *, bool);
void cleanup_worker();
void worker_run(float*);

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

static inline void check_arg(int argc, char** argv, int arg, int n) {
    if (arg + n >= argc)
        error("Option '", argv[arg], "' expects ", n, " arguments, got ", argc - arg);
}

static inline void usage() {
    std::cout << "Usage: rodent [options]\n"
              << "Available options:\n"
              << "   --help              Shows this message\n"
              << "   --width  pixels     Sets the viewport horizontal dimension (in pixels)\n"
              << "   --height pixels     Sets the viewport vertical dimension (in pixels)\n"
              << "   --eye    x y z      Sets the position of the camera\n"
              << "   --dir    x y z      Sets the direction vector of the camera\n"
              << "   --up     x y z      Sets the up vector of the camera\n"
              << "   --fov    degrees    Sets the horizontal field of view (in degrees)\n"
              << "   --bench  iterations Enables benchmarking mode and sets the number of iterations\n"
              << "   -o       image.png  Writes the output image to a file" << std::endl;
}


int main(int argc, char** argv) {
    std::string out_file;
    size_t bench_iter = 0;
    size_t width  = 1024;
    size_t height = 1024;
    float fov = 60.0f;
    bool imageDecompose = false;
    int granularity = 1; 
    float3 eye(0.0f, 1.0f, 2.5f),dir(0.0f, 0.0f, -1.0f), up(0.0f, 1.0f, 0.0f);   //cbox
//    float3 eye(-3.0f, 2.0f, -7.0f), dir(0.0f, 0.0f, 1.0f), up(0.0f, 1.0f, 0.0f);
    printf("strat\n");
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            if (!strcmp(argv[i], "--width")) {
                check_arg(argc, argv, i, 1);
                width = strtoul(argv[++i], nullptr, 10);
            } else if (!strcmp(argv[i], "--height")) {
                check_arg(argc, argv, i, 1);
                height = strtoul(argv[++i], nullptr, 10);
            } else if (!strcmp(argv[i], "--eye")) {
                check_arg(argc, argv, i, 3);
                eye.x = strtof(argv[++i], nullptr);
                eye.y = strtof(argv[++i], nullptr);
                eye.z = strtof(argv[++i], nullptr);
            } else if (!strcmp(argv[i], "--dir")) {
                check_arg(argc, argv, i, 3);
                dir.x = strtof(argv[++i], nullptr);
                dir.y = strtof(argv[++i], nullptr);
                dir.z = strtof(argv[++i], nullptr);
            } else if (!strcmp(argv[i], "--up")) {
                check_arg(argc, argv, i, 3);
                up.x = strtof(argv[++i], nullptr);
                up.y = strtof(argv[++i], nullptr);
                up.z = strtof(argv[++i], nullptr);
            } else if (!strcmp(argv[i], "--fov")) {
                check_arg(argc, argv, i, 1);
                fov = strtof(argv[++i], nullptr);
            } else if (!strcmp(argv[i], "--bench")) {
                check_arg(argc, argv, i, 1);
                bench_iter = strtoul(argv[++i], nullptr, 10);
            } else if (!strcmp(argv[i], "-o")) {
                check_arg(argc, argv, i, 1);
                out_file = argv[++i];
            } else if (!strcmp(argv[i], "-g")) {
                check_arg(argc, argv, i, 1);
                granularity = strtoul(argv[++i], nullptr, 10);
            } else if (!strcmp(argv[i], "--help")) {
                usage();
                return 0;
            } else if (!strcmp(argv[i], "--grid")){
                check_arg(argc, argv, i, 1);
                imageDecompose = strtoul(argv[++i], nullptr, 10);
            } else {
                error("Unknown option '", argv[i], "'");
            }
            continue;
        }
        error("Unexpected argument '", argv[i], "'");
    }
   

    // Force flush to zero mode for denormals
#if defined(__x86_64__) || defined(__amd64__) || defined(_M_X64)
    _mm_setcsr(_mm_getcsr() | (_MM_FLUSH_ZERO_ON | _MM_DENORMALS_ZERO_ON));
#endif
    // render 
    printf("interface\n");
    setup_interface(width, height);
    auto film = get_pixels();
    auto spp = get_spp();
    uint64_t timing = 0;
    uint32_t frames = 0;
    uint32_t iter = 0;
    
    // mpi
    struct Communicator comm;
    bool   use_master = comm.size % 2;
    printf("comm size %d\n", comm.size); 
    int    worker_size = use_master ? comm.size - 1: comm.size;
    int    master_id = worker_size; 
    int    rank = comm.rank;
    ProcStatus rs(eye, dir, up, fov, width, height, spp, rank, comm.size, 
            imageDecompose, get_chunk_num(), false/*ray queuing*/, get_dev_num(),  use_master);
    printf("after construct rendersettings \n "); 
    if(rank == master_id) {
        setup_master(&comm, &rs);
    } else {
        setup_worker(&comm, &rs, use_master); 
    }
    
    // time 
    std::vector<double> samples_sec;
    float *processTime = new float[worker_size];
    for(int i = 0; i < worker_size; i++){
        processTime[i] = 1;
    }
    float elapsed_ms = 1;  //avoid div 0
    
    //film 
    int frame = 0;
    int pixel_num = width * height * 3;
    float *reduce_buffer = new float[pixel_num];
    printf("beforether %d\n", rank);
    
    while(frame < 1) {
        clear_pixels();
        MPI_Allgather(&elapsed_ms, 1, MPI_FLOAT, processTime, 1, MPI_FLOAT, MPI_COMM_WORLD);
//        rs.camera_rotate(0.3f, 0.0f);
        auto ticks = std::chrono::high_resolution_clock::now();
        if(rank == master_id) {
            printf("start master\n");
            master_run();
        } else {
            worker_run(processTime);
        }
        printf("end init\n");
        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - ticks).count();
        samples_sec.emplace_back(1000.0 * double(spp * width * height) / double(elapsed_ms));
        film = get_pixels();
        if(rank != master_id) {
            printf("%d before reduce\n", rank);
            comm.reduce_image(film, reduce_buffer, pixel_num, use_master);
        }
        printf("%d end reduce\n", rank);

        printf("process %d time %f", rank, float(elapsed_ms));
        
        frames++;
        std::string out = out_file + "f" + std::to_string(frames) + "w" + std::to_string(worker_size) + "g" + std::to_string(get_chunk_num());
        out += ".png";
        if (rank == 0 && out_file != "") {
            float total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - ticks).count();
            save_image(reduce_buffer, out, width, height, 1 /* iter*/ );
            printf("0 node proc time %f \n", total_ms);
        }
        frame ++;
    }
    cleanup_interface();
    auto inv = 1.0e-6;
    std::sort(samples_sec.begin(), samples_sec.end());
    info("# ", samples_sec.front() * inv,
         "/", samples_sec[samples_sec.size() / 2] * inv,
         "/", samples_sec.back() * inv,
         " (min/med/max Msamples/s)");
    printf("%d end\n", rank);  
    if(rank == master_id) 
        cleanup_master();
    else 
        cleanup_worker();
    
    return 0;
}
