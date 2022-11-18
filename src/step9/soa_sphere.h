#pragma once

// generated file!

#include <stdint.h>
#include <assert.h>

#include "mymath.h"

class Material;



class SphereSOA
{
public:
    SphereSOA();
    ~SphereSOA();

    void cleanUp();

    uint32_t add(Vec3 center, float radius, Material *material);
    void copy(uint32_t dst_idx, uint32_t src_id);
    void remove(uint32_t idx);                      // ! swap with last and pop
    void reserve(uint32_t capacity);

    inline uint32_t getCount() const
    {
        return _data._count;
    }

    inline void setCount(uint32_t new_size)         // decrement allowed only; use add() to add an element
    {
        assert(new_size <= _data._count);
        _data._count = new_size;
    }

    struct InstanceData
    {
        enum { Size =  + sizeof(float) + sizeof(float) + sizeof(float) + sizeof(float) + sizeof(float) + sizeof(float) + sizeof(Material *) };

        void *_buffer;
        uint32_t _count;
        uint32_t _capacity;

        
        float * __restrict center_x;
        float * __restrict center_y;
        float * __restrict center_z;
        float * __restrict radius_sq;
        float * __restrict inv_radius;
        float * __restrict discriminant;
        Material ** __restrict material;
    };

    InstanceData *getData()
    {
        return &_data;
    }

    const InstanceData *getData() const
    {
        return &_data;
    }

private:
    void allocate(uint32_t new_count);

    InstanceData _data;
};
