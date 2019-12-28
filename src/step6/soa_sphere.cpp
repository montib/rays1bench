// generated file!

#include "soa_sphere.h"
#include "allocator.h"
#include <string.h>


static inline bool isPower2(size_t x)
{
    return x && !(x & (x - 1));
}

static inline char *align_ptr(char *ptr, size_t align)
{
    assert(align > 0 && isPower2((uint32_t)align));

    return (char*)(((uintptr_t)ptr + (align-1)) & ~(uintptr_t)(align-1));
}



SphereSOA::SphereSOA()
{
    memset(&_data, 0, sizeof(_data));
}


SphereSOA::~SphereSOA()
{
    cleanUp();
}

void SphereSOA::cleanUp()
{
    if (_data._buffer)
    {
        aligned_free(_data._buffer);
        _data._buffer = nullptr;
        _data.center = nullptr;
        _data.radius = nullptr;
        _data.discriminant = nullptr;
        _data.material = nullptr;
    }

    memset(&_data, 0, sizeof(_data));
}

void SphereSOA::copy(uint32_t dst_idx, uint32_t src_idx)
{
    _data.center[dst_idx] = _data.center[src_idx];
    _data.radius[dst_idx] = _data.radius[src_idx];
    _data.discriminant[dst_idx] = _data.discriminant[src_idx];
    _data.material[dst_idx] = _data.material[src_idx];
}

void SphereSOA::remove(uint32_t idx)
{
    --_data._count;
    uint32_t last_idx = _data._count;
    uint32_t del_idx = idx;

    if (del_idx != last_idx)
        copy(del_idx, last_idx);
}

uint32_t SphereSOA::add(Vec3 center, float radius, Material *material)
{
    reserve(_data._count+1);

    uint32_t idx = _data._count++;

    // init
    _data.center[idx] = center;
    _data.radius[idx] = radius;
    _data.discriminant[idx] = 0;
    _data.material[idx] = material;

    return idx;
}

void SphereSOA::reserve(uint32_t capacity)
{
    if (_data._capacity < capacity)
        allocate(1 + 2*capacity);
}

void SphereSOA::allocate(uint32_t new_capacity)
{
    assert(new_capacity > _data._count);    // it's a simple container :)

    const size_t Align = 32;
    InstanceData new_data;
    const size_t new_bytes = new_capacity * InstanceData::Size + Align * 4;

    new_data._buffer = aligned_allocate(new_bytes, Align);
    new_data._capacity = new_capacity;
    new_data._count = _data._count;
    
    new_data.center = (Vec3*)new_data._buffer;
    new_data.radius = (float*)(align_ptr((char*)&new_data.center[new_capacity], Align));
    new_data.discriminant = (float*)(align_ptr((char*)&new_data.radius[new_capacity], Align));
    new_data.material = (Material **)(align_ptr((char*)&new_data.discriminant[new_capacity], Align));

    if (_data._count)
    {
        memcpy(new_data.center, _data.center, _data._count * sizeof(Vec3));
        memcpy(new_data.radius, _data.radius, _data._count * sizeof(float));
        memcpy(new_data.discriminant, _data.discriminant, _data._count * sizeof(float));
        memcpy(new_data.material, _data.material, _data._count * sizeof(Material *));
    }

    if (_data._buffer)
       aligned_free(_data._buffer);

    _data = new_data;
}

