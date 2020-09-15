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
    printf("before read buffer %d mb\n", physical_memory_used_by_process() / 1024);
    size_t in_size = 0, out_size = 0;
    is.read((char*)&in_size,  sizeof(size_t));
    is.read((char*)&out_size, sizeof(size_t));
    std::vector<char> in(out_size);
    is.read(in.data(), in.size());
    array = std::move(Array(in_size / sizeof(array[0])));
    printf("before decompress memory use %d mb array %ld insize %ld outsize %ld\n", physical_memory_used_by_process() / 1024, in_size / sizeof(array[0]), in_size, out_size);
    decompress(in, array);
    printf("end decompress \n");
}

template <typename Array>
static void read_buffer(const std::string& file_name, Array& array) {
    std::ifstream is(file_name, std::ios::binary);
    read_buffer(is, array);
}

template <typename Array>
static void compress(const Array& in, std::vector<char>& out) {
    size_t in_size = sizeof(in[0]) * in.size();
    printf("in_size0 %ld\n", in_size);
    size_t out_size = LZ4_compressBound(in_size);
    printf("out size1 %ld\n", out_size);
    out.resize(out_size);
    out_size = LZ4_compress_default((const char*)in.data(), out.data(), in_size, out.size());
    printf("out size2 %ld out %ld\n", out_size, out.size());
    out.resize(out_size);
}

template <typename Array>
static void compress(const Array& in, size_t size, std::vector<char>& out) {
    size_t in_size = sizeof(in[0]) * size;
    out.resize(LZ4_compressBound(in_size));
    out.resize(LZ4_compress_default((const char*)in, out.data(), in_size, out.size()));
}

template <typename Array>
static void write_buffer(std::ostream& os, const Array& array) {
    std::vector<char> out;
    printf("before compress %ld %ld", array.size(), out.size());
    compress(array, out);
    printf("compress %ld %ld", array.size(), out.size());
    size_t in_size  = sizeof(array[0]) * array.size();
    size_t out_size = out.size();
    printf("before compress memory use %d mb array in %ld out %ld sizeof size_t %d\n", physical_memory_used_by_process() / 1024, in_size, out_size, sizeof(size_t));
    os.write((char*)&in_size,  sizeof(size_t));
    os.write((char*)&out_size, sizeof(size_t));
    os.write(out.data(), out.size());
}

template <typename Array>
static void write_buffer(const std::string& file_name, const Array& array) {
    std::ofstream of(file_name, std::ios::binary);
    write_buffer(of, array);
}

#endif // BUFFER_H
