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
#include "distribution.h"
#include "communicator.h"
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

    Camera(const float3& e, const float3& d, const float3& u, float fov, float ratio) {
        eye = e;
        dir = normalize(d);
        right = normalize(cross(dir, u));
        up = normalize(cross(right, dir));

        w = std::tan(fov * pi / 360.0f);
        h = w / ratio;
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
float* get_first_primary();
void clear_pixels();
void cleanup_interface();
void setup_server(struct Communicator *);
void server_run();
void setup_client(struct Communicator *, bool);
void client_run(Settings, int);

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
                splitimage = strtoul(argv[++i], nullptr, 10);
            } else {
                error("Unknown option '", argv[i], "'");
            }
            continue;
        }
        error("Unexpected argument '", argv[i], "'");
    }
   
    Camera cam(eye, dir, up, fov, (float)width / (float)height);

    // Force flush to zero mode for denormals
#if defined(__x86_64__) || defined(__amd64__) || defined(_M_X64)
    _mm_setcsr(_mm_getcsr() | (_MM_FLUSH_ZERO_ON | _MM_DENORMALS_ZERO_ON));
#endif
    // render 
    printf("interface\n");
    setup_interface(width, height);
    auto chunk_num = get_chunk_num();
    auto film = get_pixels();
    auto spp_g = get_spp();
    auto dev_num = get_dev_num();
    uint64_t timing = 0;
    uint32_t frames = 0;
    uint32_t iter = 0;
    
    // mpi
    struct Communicator comm;
    bool c_s = true; 
    int  client_size = c_s ? comm.size - 1: comm.size;
    int server_id = client_size; 
     
    if(comm.rank == server_id) {
        setup_server(&comm);
    } else {
        setup_client(&comm, c_s); 
    }
    
    // time 
    std::vector<double> samples_sec;
    float *proc_time = new float[client_size];
    for(int i = 0; i < client_size; i++){
        proc_time[i] = 1;
    }
    float elapsed_ms = 1;  //avoid div 0
    
    //film 
    ImageTile tiles  = ImageTile(width, height, client_size);
    int frame = 0;
    float *reduce_buffer = NULL;
    int pixel_num = width * height * 3;
    int spp_cur = spp_g;
    if (comm.rank==0)  // !!!
        reduce_buffer = new float[pixel_num];
    while(frame < 1) {
        clear_pixels();
        MPI_Allgather(&elapsed_ms, 1, MPI_FLOAT, proc_time, 1, MPI_FLOAT, MPI_COMM_WORLD);
        auto ticks = std::chrono::high_resolution_clock::now();
        if(comm.rank == server_id) {
            printf("start server\n");
            server_run();
        } else {
            float range[] = {0.0, 0.0,float(width), float(height)}; 
            if(splitimage){
                spp_cur = spp_g;
                tiles.get_tile(comm.rank, client_size, proc_time, range); 
            }
            else{
                float total = 0;
                for(int i = 0; i < client_size; i++){
                    total += proc_time[i];
                    printf("proc [] %d ", proc_time[i]);
                }
                float average = total / client_size; 
                spp_cur = spp_cur / client_size;//* average / proc_time[comm.rank];
                printf("spp_cur %d %d average%d total%d \n", spp_cur, spp_g, average, total);
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
            client_run(settings, dev_num);
            
            cam.rotate(-0.1f, 0.0f);
        }
        film = get_pixels();
        if(comm.rank != server_id) {
            printf("before reduce\n");
            comm.Reduce_image(film, reduce_buffer, pixel_num, c_s);
        }

        printf("proc %d time %f", comm.rank, float(elapsed_ms));
        
        samples_sec.emplace_back(1000.0 * double(spp_g * width * height) / double(elapsed_ms));
        frames++;
        std::string out = out_file + std::to_string(frames);
        out += ".png";
        if (comm.rank == 0 && out_file != "") {
            float total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - ticks).count();
            printf("2\n");
            save_image(reduce_buffer, out, width, height, 1 /* iter*/ );
            printf("0 node proc time %f \n", total_ms);
        }
        elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - ticks).count();
        frame ++;
    }
    delete [] proc_time; 
    cleanup_interface();
    auto inv = 1.0e-6;
    std::sort(samples_sec.begin(), samples_sec.end());
    info("# ", samples_sec.front() * inv,
         "/", samples_sec[samples_sec.size() / 2] * inv,
         "/", samples_sec.back() * inv,
         " (min/med/max Msamples/s)");
    return 0;
}
