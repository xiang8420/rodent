#pragma once
#include <vector>
#include <memory>
#include <mutex>
#include "float3.h"
#include "bbox.h"
#include "obj.h"

struct Tile{
    int children[2];
    float w;   //area weight 
    float t;    // time 
    float axis;  // axis

    Tile(float w, float t, float axis): w(w), t(t), axis(axis) {}
    bool isleaf(){return children[0] == -1 && children[1] == -1;}
};

// Grid hierarchy for scheduling
struct TileScheduler {
    std::vector<Tile> grids;
    Tile *data;
    int count;
    int worker_mpi_rank, worker_mpi_size;
    int width, height, dep;
    

    int build(int axis, int iter){
        if(iter == -1) // or w  
           return -1; 
        Tile node(0.5, 0.0, axis);
        grids.emplace_back(node);
        int cur = count;
        ++count;
        grids.at(cur).children[0] = build(1 - axis, iter - 1);
        grids.at(cur).children[1] = build(1 - axis, iter - 1);
        return cur;
    
    }
    TileScheduler(){}
    TileScheduler(int w, int h, int mpi_rank, int mpi_size)
            : width(w), height(h), worker_mpi_rank(mpi_rank), worker_mpi_size(mpi_size) {

        //scale grid and build gridtree
        count = 0;
        grids.reserve(100);
        data = grids.data();
        if(worker_mpi_size == 1) dep = 0;
        else if(worker_mpi_size == 4) dep = 2;
        else if(worker_mpi_size == 2) dep = 1;
        else dep = 4;
        build(0, dep);
        data[0].w = 1.0;
    }
    float treeupdate(int cur, int &leafnum, float* proc_time, int axis, int iter){
        Tile &node = grids.at(cur);
        if (iter == 0){
            node.t = proc_time[leafnum++];  //privious w 
        }
        else{
            Tile &child1 = grids.at(node.children[0]);
            Tile &child2 = grids.at(node.children[1]);

            //time
            float t1 = treeupdate(node.children[0], leafnum, proc_time, 1 - axis, iter - 1);
            float t2 = treeupdate(node.children[1], leafnum, proc_time, 1 - axis, iter - 1);

            // updata child weight
            child1.w = child1.w + (t2 - t1) / (t2 + t1) * 0.5;
//            child1.w = child1.w * (1 + (t2 - t1) / t1 * 0.5);
            child1.w = (child1.w + 0.5) * 0.5;
            
            child2.w = 1.0 - child1.w;
            node.t = child1.t + child2.t;
        }
        return node.t;
    }
    void get_tile(float* proc_time, int *region){
        float search[] = {0, float(worker_mpi_size)}; 
        int leaf = 0;
        treeupdate(0, leaf, proc_time, 0, dep);
        Tile *curnode = &grids.at(0);
        int axis = 0;
        while(1){
            float mid = (search[0] + search[1]) / 2;
            Tile* child;
            if(worker_mpi_rank < mid){
                child = &grids.at(curnode->children[0]);
                search[1] = mid;
                region[axis + 2] = region[axis] + (region[axis + 2] - region[axis]) * child->w;
            } else{
                child = &grids.at(curnode->children[1]);
                search[0] = mid;
                region[axis] = region[axis + 2] - (region[axis + 2] - region[axis]) * child->w;
            }
            curnode = child;
            axis = 1 - axis;
            printf("%d %d %d %d %d\n", worker_mpi_rank, region[0], region[1], region[2], region[3]);
            if(search[0] + 1 == search[1] ) break;
        }
    }

    int get_camera_info(int* region, int spp, float* processTime, bool imageDecompose) {
        int sppTask = spp;
        region[0] = 0;     region[1] = 0;
        region[2] = width; region[3] = height;
        if(imageDecompose) {
            sppTask = spp;
            get_tile(processTime, region); 
        } else {
            float total = 0;
            for(int i = 0; i < worker_mpi_size; i++){
                total += processTime[i];
                printf("proc [] %f ", processTime[i]);
            }
            float average = total / worker_mpi_size; 
            sppTask = sppTask / worker_mpi_size;//* average / processTime[comm.mpi_rank];
            printf("sppTask %d %d average%f total%f \n", sppTask, spp, average, total);
        }
        return sppTask;
    }
};

//Splitting bounding box for distributed computing
struct int6 {
    union{
        int x0, x1, y0, y1, z0, z1;
        int values[6];
    };
    int6(int *array){
        x0 = array[0]; x1 = array[1];
        y0 = array[2]; y1 = array[3];
        z0 = array[4]; z1 = array[5];
    }
    int6(){
        for(int i = 0; i < 6; i++){
            values[i] = -1;   
        }
    }
    int operator [] (size_t i) const { return values[i]; }
};
struct MeshChunk{
    std::vector<BBox>      list;
    std::vector<int6>      neighbors;  // 6 * chunk number
    BBox                   bbox;
    float3                 scale;
    float3                 step; 
    size_t size(){
        return list.size();
    }
    MeshChunk(){}
    MeshChunk(obj::File& file, float x, float y, float z) {
        chunk_division(file.bbox, x, y, z);
    }
    
    bool chunk_division(BBox bb, float scale_x, float scale_y, float scale_z){
        bbox = bb;
        printf("min %f %f %f max %f %f %f\n", bbox.min.x, bbox.min.y, bbox.min.z, bbox.max.x, bbox.max.y, bbox.max.z);
        
        // shortest axis and cut it
        float3 length = {bbox.max.x - bbox.min.x, bbox.max.y - bbox.min.y, bbox.max.z - bbox.min.z};
        int axis = 0;
        float minlength = length[0];
        for(int i = 1; i< 3; i++){
           if(minlength > length[i]){
               minlength = length[i];
               axis = i;
           }
        }
        scale[0] = scale_x;
        scale[1] = scale_y;
        scale[2] = scale_z;
        for(int i = 0; i < 3; i++){
            step[i] = length[i] / scale[i];
        }

        for(int i = 0; i < scale[0]; i++){
            float x = bbox.min.x + i * step.x;
            for(int j = 0; j < scale[1]; j++){
                float y = bbox.min.y + j * step.y;
                for(int k = 0; k < scale[2]; k++){
                    float z = bbox.min.z + k * step.z;
                    struct BBox bb;
                    bb.min.x = x - 0.00001f;
                    bb.min.y = y - 0.00001f;
                    bb.min.z = z - 0.00001f;
                    bb.max.x = ((i == scale[0] - 1)? bbox.max.x:x + step.x) + 0.0001f;
                    bb.max.y = ((j == scale[1] - 1)? bbox.max.y:y + step.y) + 0.0001f;
                    bb.max.z = ((k == scale[2] - 1)? bbox.max.z:z + step.z) + 0.0001f;
                    printf("min %f %f %f max %f %f %f\n", bb.min.x, bb.min.y, bb.min.z, bb.max.x, bb.max.y, bb.max.z);
                    list.push_back(bb);  //emplace_back
                }
            }
        }
        return true;
    }
};



