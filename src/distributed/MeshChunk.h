int get_max_length(float* length, int d) {
    float max = 0;
    int id = -1;
    for(int i = 0; i < d; i++) { 
        if( i==1 ) continue;
        if(length[i] > max) { 
            id = i;
            max = length[i];
        }
    }
    return id;
}

void splat(size_t n, float* grid, float* length, int d) {
    for(int i = 0; i < d; i++) 
        grid[i] = 1; 

    int axit = get_max_length(length, d);
    int cur_n = n;
    //choose longest axit splat
    while(cur_n != 1){
        grid[axit] *= 2;
        length[axit] /= 2;
        axit = get_max_length(length, d);
        cur_n /= 2;
    }
}

//Splitting bounding box for distributed computing

struct MeshChunk{
    std::vector<BBox>      list;
    BBox                   bbox;
    float3                 scale;
    float3                 step; 
    int                    size;

    MeshChunk(BBox bbox, int grid_num):bbox(bbox) {
        float length[] = { bbox.max.x - bbox.min.x
                         , bbox.max.y - bbox.min.y
                         , bbox.max.z - bbox.min.z} ;
        splat(grid_num, &scale[0], &length[0], 3);
        chunk_division();
        size = scale[0] * scale[1] * scale[2];
    }
    
    MeshChunk() {
        bbox  = BBox(get_bbox());
        scale = float3(get_chunk());
        chunk_division();
        size = get_chunk_num();
    }
    
    bool chunk_division() {
//        printf("mesh div min %f %f %f max %f %f %f\n", bbox.min.x, bbox.min.y, bbox.min.z, bbox.max.x, bbox.max.y, bbox.max.z);
        
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
                    bb.min.x = x - 0.001f;
                    bb.min.y = y - 0.001f;
                    bb.min.z = z - 0.001f;
                    bb.max.x = std::min(bbox.max.x, x + step.x + 0.001f);
                    bb.max.y = std::min(bbox.max.y, y + step.y + 0.001f);
                    bb.max.z = std::min(bbox.max.z, z + step.z + 0.001f);
//                    printf("min %f %f %f max %f %f %f\n", bb.min.x, bb.min.y, bb.min.z, bb.max.x, bb.max.y, bb.max.z);
                    list.push_back(bb);  //emplace_back
                }
            }
        }
        return true;
    }
};
