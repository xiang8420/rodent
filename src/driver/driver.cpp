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

#include <thread>

#include "interface.h"
#include "float3.h"
#include "common.h"
#include "image.h"
#include "grid.h"
#if defined(__x86_64__) || defined(__amd64__) || defined(_M_X64)
#include <x86intrin.h>
#endif

static constexpr float pi = 3.14159265359f;

struct Camera {
    float3 eye;
    float3 dir;
    float3 right;
    float3 up;
    float w, h;
    ImageGridTree *tree;

    Camera(const float3& e, const float3& d, const float3& u, float fov, float ratio, float width, float height, int proc_num) {
        eye = e;
        dir = normalize(d);
        right = normalize(cross(dir, u));
        up = normalize(cross(right, dir));

        w = std::tan(fov * pi / 360.0f);
        h = w / ratio;
        tree = new ImageGridTree(width, height, proc_num);
    }

    void rotate(float yaw, float pitch) {
        dir = ::rotate(dir, right,  -pitch);
        dir = ::rotate(dir, up,     -yaw);
        dir = normalize(dir);
        right = normalize(cross(dir, up));
        up = normalize(cross(right, dir));
    }

    void move(float x, float y, float z) {
        eye += right * x + up * y + dir * z;
    }
    
};

void setup_interface(size_t, size_t);
float* get_pixels();
void clear_pixels();
void cleanup_interface();

#ifndef DISABLE_GUI
static bool handle_events(uint32_t& iter, Camera& cam) {
    static bool camera_on = false;
    static bool arrows[4] = { false, false, false, false };
    static bool speed[2] = { false, false };
    const float rspeed = 0.005f;
    static float tspeed = 0.1f;

    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        bool key_down = event.type == SDL_KEYDOWN;
        switch (event.type) {
            case SDL_KEYUP:
            case SDL_KEYDOWN:
                switch (event.key.keysym.sym) {
                    case SDLK_ESCAPE:   return true;
                    case SDLK_KP_PLUS:  speed[0] = key_down; break;
                    case SDLK_KP_MINUS: speed[1] = key_down; break;
                    case SDLK_UP:       arrows[0] = key_down; break;
                    case SDLK_DOWN:     arrows[1] = key_down; break;
                    case SDLK_LEFT:     arrows[2] = key_down; break;
                    case SDLK_RIGHT:    arrows[3] = key_down; break;
                }
                break;
            case SDL_MOUSEBUTTONDOWN:
                if (event.button.button == SDL_BUTTON_LEFT) {
                    SDL_SetRelativeMouseMode(SDL_TRUE);
                    camera_on = true;
                }
                break;
            case SDL_MOUSEBUTTONUP:
                if (event.button.button == SDL_BUTTON_LEFT) {
                    SDL_SetRelativeMouseMode(SDL_FALSE);
                    camera_on = false;
                }
                break;
            case SDL_MOUSEMOTION:
                if (camera_on) {
                    cam.rotate(event.motion.xrel * rspeed, event.motion.yrel * rspeed);
                    iter = 0;
                }
                break;
            case SDL_QUIT:
                return true;
            default:
                break;
        }
    }

    if (arrows[0]) cam.move(0, 0,  tspeed);
    if (arrows[1]) cam.move(0, 0, -tspeed);
    if (arrows[2]) cam.move(-tspeed, 0, 0);
    if (arrows[3]) cam.move( tspeed, 0, 0);
    if (arrows[0] | arrows[1] | arrows[2] | arrows[3]) iter = 0;
    if (speed[0]) tspeed *= 1.1f;
    if (speed[1]) tspeed *= 0.9f;
    return false;
}

static void update_texture(uint32_t* buf, SDL_Texture* texture, size_t width, size_t height, uint32_t iter) {
    auto film = get_pixels();
    auto inv_iter = 1.0f / iter;
    auto inv_gamma = 1.0f / 2.2f;
    for (size_t y = 0; y < height; ++y) {
        for (size_t x = 0; x < width; ++x) {
            auto r = film[(y * width + x) * 3 + 0];
            auto g = film[(y * width + x) * 3 + 1];
            auto b = film[(y * width + x) * 3 + 2];

            buf[y * width + x] =
                (uint32_t(clamp(std::pow(r * inv_iter, inv_gamma), 0.0f, 1.0f) * 255.0f) << 16) |
                (uint32_t(clamp(std::pow(g * inv_iter, inv_gamma), 0.0f, 1.0f) * 255.0f) << 8)  |
                 uint32_t(clamp(std::pow(b * inv_iter, inv_gamma), 0.0f, 1.0f) * 255.0f);
        }
    }
    SDL_UpdateTexture(texture, nullptr, buf, width * sizeof(uint32_t));
}
#endif

static void save_image(float *result, const std::string& out_file, size_t width, size_t height, uint32_t iter) {
    ImageRgba32 img;
    img.width = width;
    img.height = height;
    img.pixels.reset(new uint8_t[width * height * 4]);

//    auto film = get_pixels();
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
    bool splitimage = false;
    int dev_num = 1;
    float3 eye(0.0f, 1.0f, 2.5f),dir(0.0f, 0.0f, -1.0f), up(0.0f, 1.0f, 0.0f);   //cbox
//    float3 eye(-3.0f, 2.0f, -7.0f), dir(0.0f, 0.0f, 1.0f), up(0.0f, 1.0f, 0.0f);
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
            } else if (!strcmp(argv[i], "--help")) {
                usage();
                return 0;
            } else if (!strcmp(argv[i], "--grid")){
                check_arg(argc, argv, i, 1);
                splitimage = strtoul(argv[++i], nullptr, 10);
            } else if (!strcmp(argv[i], "--dev")){
                check_arg(argc, argv, i, 1);
                dev_num = strtoul(argv[++i], nullptr, 10);
            } else {
                error("Unknown option '", argv[i], "'");
            }
            continue;
        }
        error("Unexpected argument '", argv[i], "'");
    }
    int  world_rank = -1, world_size = -1;
    MPI_Init(NULL, NULL); 
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    Camera cam(eye, dir, up, fov, (float)width / (float)height, (float)width, (float)height, world_size);

#ifdef DISABLE_GUI
//    info("Running in console-only mode (compiled with -DDISABLE_GUI).");
//    if (bench_iter == 0) {
//        warn("Benchmark iterations no set. Defaulting to 1.");
//        bench_iter = 1;
//    }
#else
    if (SDL_Init(SDL_INIT_VIDEO) != 0)
        error("Cannot initialize SDL.");

    auto window = SDL_CreateWindow(
        "Rodent",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        width,
        height,
        0);
    if (!window)
        error("Cannot create window.");

    auto renderer = SDL_CreateRenderer(window, -1, 0);
    if (!renderer)
        error("Cannot create renderer.");

    auto texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STATIC, width, height);
    if (!texture)
        error("Cannot create texture");

    std::unique_ptr<uint32_t> buf(new uint32_t[width * height]);
#endif

    setup_interface(width, height);

    auto film = get_pixels();
    // Force flush to zero mode for denormals
#if defined(__x86_64__) || defined(__amd64__) || defined(_M_X64)
    _mm_setcsr(_mm_getcsr() | (_MM_FLUSH_ZERO_ON | _MM_DENORMALS_ZERO_ON));
#endif

    auto spp = get_spp();
    bool done = false;
    uint64_t timing = 0;
    uint32_t frames = 0;
    uint32_t iter = 0;
    std::vector<double> samples_sec;
    
    float *proc_time = new float[world_size];
    int spp_old = spp / world_size;
    float elapsed_ms = 1;  //avoid div 0
//    cam.rotate(-0.25f, 0.0f);
    float *reduce_buffer = NULL;
    int pixel_num = width * height * 3;
    if (world_rank==0)
        reduce_buffer = new float[pixel_num];
    if(splitimage){
        printf("distributed image grid\n");
    } else {
        printf("distributed spp \n");
    }
    std::thread *thread[5];
    while (frames < 7) {
        if (iter == 0)
            clear_pixels();
        cam.rotate(0.05f, 0.0f);
        MPI_Allgather(&elapsed_ms, 1, MPI_FLOAT, proc_time, 1, MPI_FLOAT, MPI_COMM_WORLD);
        float range[] = {0.0, 0.0, width, height}; 
        int spp_cur;
        if(splitimage){
            spp_cur = spp;
            cam.tree->getgrid(world_rank, world_size, proc_time, range); 
        }
        else{
            float total = 0;
            for(int i = 0; i < world_size; i++){
                total += proc_time[i];
            }
            float average = total / world_size; //average speed
            spp_cur = spp_old / proc_time[world_rank] * average;
        }
        Settings settings {
            Vec3 { cam.eye.x, cam.eye.y, cam.eye.z },
            Vec3 { cam.dir.x, cam.dir.y, cam.dir.z },
            Vec3 { cam.up.x, cam.up.y, cam.up.z },
            Vec3 { cam.right.x, cam.right.y, cam.right.z },
            cam.w,
            cam.h,
            Vec4 { range[0], range[1], range[2], range[3]},
            spp_cur
        };

#ifdef DEBUG    
        printf("range %f %f %f %f %d \n", range[0], range[1], range[2], range[3], spp_cur);
        for (int i = 0; i < world_size; i++)
                printf("%f|",proc_time[i]);
#endif
              
        auto ticks = std::chrono::high_resolution_clock::now();
        
        for(int i = 0; i < dev_num; i++){
            thread[i] = new std::thread(render, &settings, world_rank, i);    
        } 
        for(int i = 0; i < dev_num; i++){
            thread[i]->join();
        }
        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - ticks).count();
        
        film = get_pixels(); 
        MPI_Reduce(film, reduce_buffer, pixel_num, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
#ifdef DEBUG        
        printf("proc %d time %f", world_rank, float(elapsed_ms));
#endif
        samples_sec.emplace_back(1000.0 * double(spp * width * height) / double(elapsed_ms));
        frames++;
        timing += elapsed_ms;
    
        std::string out = out_file + std::to_string(frames);
        out += ".png";
        if (world_rank == 0 && out_file != "") {
//            info("0 image wirting\n");
            float total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - ticks).count();
       //     film = get_pixels();
            save_image(reduce_buffer, out, width, height, 1 /* iter*/ );
//            info("Image saved to '", out, "'");
            printf("0 node proc time %f \n", total_ms);
        }
    }
    delete [] proc_time; 
    MPI_Finalize(); 
    cleanup_interface();

    auto inv = 1.0e-6;
    std::sort(samples_sec.begin(), samples_sec.end());
    info("# ", samples_sec.front() * inv,
         "/", samples_sec[samples_sec.size() / 2] * inv,
         "/", samples_sec.back() * inv,
         " (min/med/max Msamples/s)");
    return 0;
}
