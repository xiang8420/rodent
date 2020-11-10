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
#include "../distributed/DistributedFrameWork.h"
#include "sdtree.h"
#if defined(__x86_64__) || defined(__amd64__) || defined(_M_X64)
#include <x86intrin.h>
#endif


void setup_interface(size_t, size_t);
float* get_pixels();
void clear_pixels();
void cleanup_interface();
float* get_first_primary();


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
    std::string distributedMode = "P2PNode";
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
            } else if (!strcmp(argv[i], "--DMode")){
                check_arg(argc, argv, i, 1);
                distributedMode = argv[++i];
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
    auto spp = get_spp();
    uint64_t timing = 0;
    uint32_t frames = 0;
    uint32_t iter = 0;
    
    // inital process setting and status  
    printf("after construct rendersettings \n "); 
    setup_distributed_framework(distributedMode, get_chunk_num(), get_dev_num(), width, height, spp); 
    
    Camera camera(eye, dir, up, fov, width, height);
    // time statistic 
    std::vector<double> samples_sec;
    float elapsed_ms = 1;  //avoid div 0
   
    bool animation = false; 
    int frame = 1;
    for(int i = 0; i < frame; i++) {
        auto ticks = std::chrono::high_resolution_clock::now();
        
        dfw_run(&camera);
        
        printf("after dfw run\n");

        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - ticks).count();
        samples_sec.emplace_back(1000.0 * double(spp * width * height) / double(elapsed_ms));
        
        float total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - ticks).count();
        printf("end reduce process time %f",  float(elapsed_ms));
        
        if(animation) {
            auto film = get_pixels();
            dfw_save_image(out_file, film, i, frame, width, height);
            clear_pixels();
            camera.rotate(0.1f, 0.0f);
        }
        camera.iter++;
    }
    if(!animation) {
        auto film = get_pixels();
        dfw_save_image(out_file, film, frame, frame, width, height);
    }

    auto inv = 1.0e-6;
    std::sort(samples_sec.begin(), samples_sec.end());
    info("# ", samples_sec.front() * inv,
         "/", samples_sec[samples_sec.size() / 2] * inv,
         "/", samples_sec.back() * inv,
         " (min/med/max Msamples/s)");
    
    cleanup_interface();
    cleanup_distributed_framework(); 
    return 0;
}
