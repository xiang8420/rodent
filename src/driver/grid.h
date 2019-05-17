#include <vector>
#define DEP 2
//Splitting bounding box for distributed computing
struct Grid3{
    float3 min;
    float3 max;
    int children[2];    
    int count;  //primtive count
    float w;    //weight to adjust

    Grid3() {}
    Grid3(BBox bbox) : min(bbox.min), max(bbox.max){}
    Grid3(float3& min, float3& max) : min(min), max(max){}
    Grid3(float3& min, float3& max, float w) : min(min), max(max), w(w){}
    Grid3(float *value){
        min.x = value[0]; min.y = value[1]; min.z = value[2];
        max.x = value[3]; max.y = value[4]; max.z = value[5];
    } 
    bool is_leaf(){return children[0] == -1 && children[1] == -1;}
    int longaxis(){
        float lenth[] = {max.x - min.x, max.y - min.y, max.z - min.z};
        float maxlenth = -FLT_MAX;
        int axis;
        for(int i = 0; i < 3; i++){
            if(lenth[i] > maxlenth){
                axis = i;
                maxlenth = lenth[i];
            }
        }
       return axis;
    }
};
struct Grid2{
    int children[2];
    float w;   //area weight 
    float t;    // time 
    float axis;  // axis

    Grid2(float w, float t, float axis): w(w), t(t), axis(axis) {}
    bool is_leaf(){return children[0] == -1 && children[1] == -1;}
};

// Grid hierarchy for scheduling
struct ImageGridTree {
    std::vector<Grid2> grids;
    int count;
    int width, height;

    int build(int axis, int iter){
        if(iter == -1) // or w  
           return -1; 
        Grid2 node(0.5, 0.0, axis);
        grids.emplace_back(node);
        int cur = count;
        ++count;
        grids[cur].children[0] = build(1 - axis, iter - 1);
        grids[cur].children[1] = build(1 - axis, iter - 1);
        return cur;
    
    }
    ImageGridTree(int w, int h): width(w), height(h){
        //split grid and build gridtree
        count = 0;
        build(0, DEP);
        grids.at(0).w = 1.0;
    }
    float treeupdate(int cur, int &leafnum, float* proc_time, int axis, int iter){
        Grid2 &node = grids.at(cur);
        if (iter == 0){
            node.t = proc_time[leafnum++];  //privious w 
        }
        else{
            Grid2 &child1 = grids.at(node.children[0]);
            Grid2 &child2 = grids.at(node.children[1]);

            //time
            float t1 = treeupdate(node.children[0], leafnum, proc_time, 1 - axis, iter - 1);
            float t2 = treeupdate(node.children[1], leafnum, proc_time, 1 - axis, iter - 1);

            // updata child weight
     //       child1.w = child1.w + (t2 - t1) / (t2 + t1) * 0.5;
            child1.w = child1.w * (1 + (t2 - t1) / t1 * 0.5);
            child1.w = (child1.w + 0.5) * 0.5;
            
            child2.w = 1.0 - child1.w;
            node.t = child1.t + child2.t;
        }
        return node.t;
    }
    void getgrid(int proc_id, int proc_num, float* proc_time, float range[4]){
        float search[] = {0, proc_num}; 
        range[0] = 0; range[1] = 0;
        range[2] = width; range[3] = height;
        int leaf = 0;
        treeupdate(0, leaf, proc_time, 0, DEP);
        Grid2 *curnode = &grids.at(0);
        int axis = 0;
        while(1){
            float mid = (search[0] + search[1]) / 2;
            Grid2* child;
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


