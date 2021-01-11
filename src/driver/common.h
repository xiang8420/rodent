#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <random>

// Round to the integer above the division.
inline uint32_t round_up(uint32_t val, uint32_t div) {
    auto mod = val % div;
    return val + (mod ? div - mod : 0);
}

/// Clamps a between b and c.
template <typename T>
inline T clamp(T a, T b, T c) {
    return (a < b) ? b : ((a > c) ? c : a);
}

/// Returns the integer that is greater or equal to the logarithm base 2 of the argument.
template <typename T>
inline T closest_log2(T i) {
    T p = 1, q = 0;
    while (i > p) p <<= 1, q++;
    return q;
}

/// Reinterprets a floating point number as an integer.
inline int32_t float_as_int(float f) {
    union { float vf; int32_t vi; } v;
    v.vf = f;
    return v.vi;
}

/// Reinterprets an integer as a floating point number.
inline float int_as_float(int32_t i) {
    union { float vf; int32_t vi; } v;
    v.vi = i;
    return v.vf;
}

inline void error [[noreturn]] () {
    std::cerr << std::endl;
    abort();
}

/// Outputs an error message in the console.
template <typename T, typename... Args>
inline void error [[noreturn]] (T t, Args... args) {
#if COLORIZE
    std::cerr << "\033[1;31m";
#endif
    std::cerr << t;
#if COLORIZE
    std::cerr << "\033[0m";
#endif
    error(args...);
}

inline void info() {
    std::cout << std::endl;
}

/// Outputs an information message in the console.
template <typename T, typename... Args>
inline void info(T t, Args... args) {
    std::cout << t;
    info(args...);
}

inline void warn() {
    std::clog << std::endl;
}

/// Outputs an warning message in the console.
template <typename T, typename... Args>
inline void warn(T t, Args... args) {
#if COLORIZE
    std::clog << "\033[1;33m";
#endif
    std::clog << t;
#if COLORIZE
    std::clog << "\033[0m";
#endif
    warn(args...);
}

inline size_t physical_memory_used_by_process() {
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != nullptr) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            int len = strlen(line);
            const char* p = line;
            for (; std::isdigit(*p) == false; ++p) {}
            line[len - 3] = 0;
            result = atoi(p);
            break;
        }
    }
    fclose(file);
    return result / 1024;
}

#define MAX_CHUNK 3
#define MAX_PROC  128
#define PRIMARY_WIDTH 22
#define SECONDARY_WIDTH 17
#define CHUNK_HIT_RES 128
#define INTERSECT_FIELD_RES 128
#define INTERSECT_FIELD_RES 128
#define VIS_CHUNK_HIT true
#define STREAM_CAPACITY 256 * 256
#define OUT_BUFFER_CAPACITY 1024 * 1024
#define MIN_SEND_STREAM_SIZE 10
#define SIMPLE_TRACE true

#endif // COMMON_H
