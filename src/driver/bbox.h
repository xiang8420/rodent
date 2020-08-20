#ifndef BBOX_H
#define BBOX_H

#include <cfloat>
#include <algorithm>
#include "float3.h"

/// Bounding box represented by its two extreme points.
struct BBox {
    float3 min, max;

    BBox() {}
    BBox(const float3& f) : min(f), max(f) {}
    BBox(const float3& min, const float3& max) : min(min), max(max) {}
    BBox(const float *value){
        min.x = value[0]; min.y = value[1]; min.z = value[2];
        max.x = value[3]; max.y = value[4]; max.z = value[5];
    }

    BBox& extend(const BBox& bb) {
        min = ::min(min, bb.min);
        max = ::max(max, bb.max);
        return *this;
    }
    
    void extend(float n) {
        min.x -= n; min.y -= n; min.z -= n;
        max.x += n; max.y += n; max.z += n;
    }
    
    void ToFloat(float *value){
        value[0] = min.x; value[1] = min.y; value[2] = min.z;
        value[3] = max.x; value[4] = max.y; value[5] = max.z;
    }
    BBox& extend(const float3& v) {
        min = ::min(min, v);
        max = ::max(max, v);
        return *this;
    }
    float half_area() const {
        const float3 len = max - min;
        const float kx = std::max(len.x, 0.0f);
        const float ky = std::max(len.y, 0.0f);
        const float kz = std::max(len.z, 0.0f);
        return kx * (ky + kz) + ky * kz;
    }
    
    BBox& overlap(const BBox& bb) {
        min = ::max(min, bb.min);
        max = ::min(max, bb.max);
        return *this;
    }

    bool is_empty() const {
        return min.x > max.x ||
               min.y > max.y ||
               min.z > max.z;
    }


    bool is_inside(const float3& v) const {
        return v.x >= min.x && v.y >= min.y && v.z >= min.z &&
               v.x <= max.x && v.y <= max.y && v.z <= max.z;
    }

    bool line_intersect(const float3& v0, const float3& v1){
        float3 dir = v1 - v0;
        float3 t0  = (min - v0) / dir;
        float3 t1  = (max - v0) / dir;
        float tmin = 0, tmax = 1;
        
        tmin = std::max(std::max(std::min(t0.x, t1.x), std::min(t0.y, t1.y)), std::min(t0.z, t1.z));
        tmax = std::min(std::min(std::max(t0.x, t1.x), std::max(t0.y, t1.y)), std::max(t0.z, t1.z)); 
        if(tmin > tmax)
            return false;
        return (tmin >= 0 && tmin <= 1) || (tmax >= 0 && tmax <= 1);
    }

    bool is_overlapping(const BBox& bb) const {
        return min.x <= bb.max.x && max.x >= bb.min.x &&
               min.y <= bb.max.y && max.y >= bb.min.y &&
               min.z <= bb.max.z && max.z >= bb.min.z;
    }
    static BBox empty() { return BBox(float3(FLT_MAX), float3(-FLT_MAX)); }
    static BBox full() { return BBox(float3(-FLT_MAX), float3(FLT_MAX)); }
};

#endif // BBOX_H
