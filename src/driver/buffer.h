#ifndef BUFFER_H
#define BUFFER_H

#include <ostream>
#include <istream>
#include <fstream>

#include <lz4.h>

static void skip_buffer(std::istream& is) {
    size_t in_size = 0, out_size = 0;
    is.read((char*)&in_size,  sizeof(uint32_t));
    is.read((char*)&out_size, sizeof(uint32_t));
    is.seekg(out_size, std::ios::cur);
}

template <typename Array>
static void decompress(const std::vector<char>& in, Array& out) {
    LZ4_decompress_safe(in.data(), (char*)out.data(), in.size(), out.size() * sizeof(out[0]));
}

template <typename Array>
static void read_buffer(std::istream& is, Array& array) {
    size_t size;
    is.read((char*)&size, sizeof(size_t));

    // decompress out 
    array = std::move(Array(size / sizeof(array[0])));
    char* out = (char*)array.data();
    size_t st = 0;
    while(st < size) {
        size_t compress_size, decompress_size;
       
        is.read((char*)&decompress_size, sizeof(size_t));
        is.read((char*)&compress_size, sizeof(size_t));
        
        std::vector<char> in;
        in.resize(compress_size); 
        is.read(in.data(), in.size());
        
        LZ4_decompress_safe(in.data(), &out[st], in.size(), decompress_size);

        st += decompress_size;
    }
}

template <typename Array>
static void read_buffer(const std::string& file_name, Array& array) {
    std::ifstream is(file_name, std::ios::binary);
    read_buffer(is, array);
}

template <typename Array>
static void compress(const Array& in, std::vector<char>& out) {
    size_t in_size = sizeof(in[0]) * in.size();
    out.resize(LZ4_compressBound(in_size));
    out.resize(LZ4_compress_default((const char*)in.data(), out.data(), in_size, out.size()));
}

template <typename Array>
static void compress(const Array& in, size_t size, std::vector<char>& out) {
    size_t in_size = sizeof(in[0]) * size;
    out.resize(LZ4_compressBound(in_size));
    out.resize(LZ4_compress_default((const char*)in, out.data(), in_size, out.size()));
}

template <typename Array>
static void write_buffer(std::ostream& os, const Array& array) {
    char* in = (char*)array.data();

    int st = 0;
    size_t size = array.size() * sizeof(array[0]); 
    os.write((char*)&size, sizeof(size_t));

    size_t lz4_max = LZ4_MAX_INPUT_SIZE;
    while(st < size) {
        size_t origin_size = size - st > lz4_max ? lz4_max : size - st; 
        printf("before compress Memory %ld array size %ld array[0 ]%ld %d\n", physical_memory_used_by_process() / 1024, array.size(), sizeof(array[0]), LZ4_compressBound(origin_size));
        std::vector<char> out;
        out.resize(LZ4_compressBound(origin_size));
        out.resize(LZ4_compress_default((const char*)&in[st], (char*)out.data(), origin_size, out.size()/*LZ4_compressBound(origin_size)*/));
        
        size_t compress_size = out.size();

        os.write((char*)&origin_size, sizeof(size_t));
        os.write((char*)&compress_size, sizeof(size_t));
        os.write(out.data(), out.size());

        st += origin_size;
    }
}

template <typename Array>
static void write_buffer(const std::string& file_name, const Array& array) {
    //std::ofstream of(file_name, std::ios::binary);
    remove(file_name.c_str());
    std::ofstream of(file_name, std::ios::binary | std::ios::app);
    write_buffer(of, array);
}

#endif // BUFFER_H
