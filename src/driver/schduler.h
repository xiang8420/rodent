#include <vector>
#include "bbox.h"
#include "obj.h"
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
    int                    scale[3];
    
    size_t size(){
        return list.size();
    }
    MeshChunk(){}
    MeshChunk(obj::File& file) {
        chunk_division(file.bbox);
    }
    
    bool chunk_division(BBox bb){
        bbox = bb;
        printf("min %f %f %f max %f %f %f\n", bbox.min.x, bbox.min.y, bbox.min.z, bbox.max.x, bbox.max.y, bbox.max.z);
        
        // shortest axis and cut it
        float lenth[] = {bbox.max.x - bbox.min.x, bbox.max.y - bbox.min.y, bbox.max.z - bbox.min.z};
        int axis = 0;
        float minlenth = lenth[0];
        for(int i = 1; i< 3; i++){
           if(minlenth > lenth[i]){
               minlenth = lenth[i];
               axis = i;
           }
        }
        ///scale num per axis set shortest axis 2, so 8 bvh at least
        scale[0] = 1;
        scale[1] = 1;
        scale[2] = 1;

        axis = 0;
        scale[axis] = 2;
        float step = lenth[axis] / 2;
        //    for(int i = 0; i < 3; i++){
        //        scale[i] = lenth[i] / step;
        //        printf("scale i %d", scale[i]);
        //    }
        printf("\n");
        for(int i = 0; i < scale[0]; i++){
            float x = bbox.min.x + i * step;
            for(int j = 0; j < scale[1]; j++){
                float y = bbox.min.y + j * step;
                for(int k = 0; k < scale[2]; k++){
                    float z = bbox.min.z + k * step;
                    struct BBox bb;
                    bb.min.x = x;
                    bb.min.y = y;
                    bb.min.z = z;
                    bb.max.x = (i == scale[0] - 1)? bbox.max.x:x + step;
                    bb.max.y = (j == scale[1] - 1)? bbox.max.y:y + step;
                    bb.max.z = (k == scale[2] - 1)? bbox.max.z:z + step;
                    printf("min %f %f %f max %f %f %f\n", bb.min.x, bb.min.y, bb.min.z, bb.max.x, bb.max.y, bb.max.z);
                    list.push_back(bb);  //emplace_back
                }
            }
        }
        return true;
    }
};

struct Tile{
    int children[2];
    float w;   //area weight 
    float t;    // time 
    float axis;  // axis

    Tile(float w, float t, float axis): w(w), t(t), axis(axis) {}
    bool isleaf(){return children[0] == -1 && children[1] == -1;}
};

// Grid hierarchy for scheduling
struct ImageTile {
    std::vector<Tile> grids;
    int count;
    int width, height;
    int dep;

    int build(int axis, int iter){
        if(iter == -1) // or w  
           return -1; 
        Tile node(0.5, 0.0, axis);
        grids.emplace_back(node);
        int cur = count;
        ++count;
        grids[cur].children[0] = build(1 - axis, iter - 1);
        grids[cur].children[1] = build(1 - axis, iter - 1);
        return cur;
    
    }
    ImageTile(){}
    ImageTile(int w, int h, int proc_num): width(w), height(h){
        //scale grid and build gridtree
        count = 0;
        if(proc_num == 4) dep = 2;
        else if(proc_num == 2) dep = 1;
        else dep = 4;
        build(0, dep);
        grids.at(0).w = 1.0;
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
    void getgrid(int proc_id, int proc_num, float* proc_time, float range[4]){
        float search[] = {0, float(proc_num)}; 
        range[0] = 0; range[1] = 0;
        range[2] = width; range[3] = height;
        int leaf = 0;
        treeupdate(0, leaf, proc_time, 0, dep);
        Tile *curnode = &grids.at(0);
        int axis = 0;
        while(1){
            float mid = (search[0] + search[1]) / 2;
            Tile* child;
            if(proc_id < mid){
                child = &grids.at(curnode->children[0]);
                search[1] = mid;
                range[axis + 2] = range[axis] + (range[axis + 2] - range[axis]) * child->w;
            } else{
                child = &grids.at(curnode->children[1]);
                search[0] = mid;
                range[axis] = range[axis + 2] - (range[axis + 2] - range[axis]) * child->w;
            }
            curnode = child;
            axis = 1 - axis;
            printf("%d %f %f %f %f\n", proc_id, range[0], range[1], range[2], range[3]);
            if(search[0] + 1 == search[1] ) break;
        }
    }
};

struct Schduler {
    ImageTile   imagetile; 
    MeshChunk   meshchunk;
    
    Schduler(float width, float height, int proc_num) {
        imagetile = ImageTile(width, height, proc_num);
    
    }
};
