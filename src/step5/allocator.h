#pragma once

#include <stdlib.h>

#ifndef _WIN32
#include <mm_malloc.h>
#endif

inline void* aligned_allocate(size_t size, size_t alignment)
    {
#ifdef _WIN32
        return _aligned_malloc(size, alignment);
#else 
        return _mm_malloc(size, alignment);
#endif
    }

inline void aligned_free(void* p)
{
#ifdef _WIN32
    _aligned_free(p);
#else 
    _mm_free(p);
#endif 
}
