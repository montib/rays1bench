#pragma once

#include <math.h>
#include <immintrin.h>


#if defined(_MSC_VER)
#define VM_INLINE __forceinline
#else
#define VM_INLINE __attribute__((unused, always_inline)) inline
#endif



//static uint32_t xorshift_state = DEFAULT_RANDOM_STATE;

VM_INLINE static uint32_t XorShift32(uint32_t& state)
{
    uint32_t x = state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 15;
    state = x;
    return x;
}

VM_INLINE float myrand01(uint32_t& state)
{
    return (XorShift32(state) & 0xFFFFFF) * float(1.0 / 16777216.0);
}
                        
VM_INLINE float myrand02(uint32_t& state)
{
    return (XorShift32(state) & 0xFFFFFF) / (float)(0xFFFFFF/2 + 1);
}


//static __m128i xorshift_state4 = DEFAULT_RANDOM4_STATE;


inline __m128 myrand01_x4(__m128i &state4)
{
    // u32 random
    __m128i x = state4;
    x = _mm_xor_si128(x, _mm_slli_epi32(x, 13));
    x = _mm_xor_si128(x, _mm_srli_epi32(x, 17));
    x = _mm_xor_si128(x, _mm_slli_epi32(x, 15));
    state4 = x;

    // convert to float
    x = _mm_and_si128(x, _mm_set1_epi32(0xFFFFFF));
    __m128 xf = _mm_cvtepi32_ps(x);
    xf = _mm_mul_ps(xf, _mm_set1_ps((float)(1.0 / (0xFFFFFF + 1))));

    return xf;
}

VM_INLINE __m128 myrand02_x4(__m128i& state4)
{
    // u32 random
    __m128i x = state4;
    x = _mm_xor_si128(x, _mm_slli_epi32(x, 13));
    x = _mm_xor_si128(x, _mm_srli_epi32(x, 17));
    x = _mm_xor_si128(x, _mm_slli_epi32(x, 15));
    state4 = x;

    // convert to float
    x = _mm_and_si128(x, _mm_set1_epi32(0xFFFFFF));
    __m128 xf = _mm_cvtepi32_ps(x);
    xf = _mm_mul_ps(xf, _mm_set1_ps((float)(1.0 / (0xFFFFFF / 2 + 1))));

    return xf;
}


#if 1

// SSE/SIMD vector largely based on http://www.codersnotes.com/notes/maths-lib-2016/

// SHUFFLE3(v, 0,1,2) leaves the vector unchanged (v.xyz).
// SHUFFLE3(v, 0,0,0) splats the X (v.xxx).
#define SHUFFLE3(V, X,Y,Z) Vec3(_mm_shuffle_ps((V).m, (V).m, _MM_SHUFFLE(Z,Z,Y,X)))


struct Vec3
{
    __m128 m;

    VM_INLINE Vec3() {}
    VM_INLINE explicit Vec3(const float* p) { m = _mm_set_ps(p[2], p[2], p[1], p[0]); }
    VM_INLINE explicit Vec3(float x, float y, float z) { m = _mm_set_ps(z, z, y, x); }
    VM_INLINE explicit Vec3(float v) { m = _mm_set1_ps(v); }
    VM_INLINE explicit Vec3(__m128 v) { m = v; }

    void setX(float x)
    {
        m = _mm_move_ss(m, _mm_set_ss(x));
    }
    void setY(float y)
    {
        __m128 t = _mm_move_ss(m, _mm_set_ss(y));
        t = _mm_shuffle_ps(t, t, _MM_SHUFFLE(3, 2, 0, 0));
        m = _mm_move_ss(t, m);
    }
    void setZ(float z)
    {
        __m128 t = _mm_move_ss(m, _mm_set_ss(z));
        t = _mm_shuffle_ps(t, t, _MM_SHUFFLE(3, 0, 1, 0));
        m = _mm_move_ss(t, m);
    }

    //VM_INLINE float getX() const { return _mm_cvtss_f32(m); }
    //VM_INLINE float getY() const { return _mm_cvtss_f32(_mm_shuffle_ps(m, m, _MM_SHUFFLE(1, 1, 1, 1))); }
    //VM_INLINE float getZ() const { return _mm_cvtss_f32(_mm_shuffle_ps(m, m, _MM_SHUFFLE(2, 2, 2, 2))); }

    //VM_INLINE float getX() const { return m.m128_f32[0]; }
    //VM_INLINE float getY() const { return m.m128_f32[1]; }
    //VM_INLINE float getZ() const { return m.m128_f32[2]; }

    VM_INLINE float getX() const { return ((float*)&m)[0]; }
    VM_INLINE float getY() const { return ((float*)&m)[1]; }
    VM_INLINE float getZ() const { return ((float*)&m)[2]; }
    

    VM_INLINE Vec3 yzx() const { return SHUFFLE3(*this, 1, 2, 0); }
    VM_INLINE Vec3 zxy() const { return SHUFFLE3(*this, 2, 0, 1); }
};

typedef Vec3 bool3;

//VM_INLINE Vec3 operator+ (Vec3 a, Vec3 b) { a.m = _mm_add_ps(a.m, b.m); return a; }
//VM_INLINE Vec3 operator- (Vec3 a, Vec3 b) { a.m = _mm_sub_ps(a.m, b.m); return a; }
//VM_INLINE Vec3 operator* (Vec3 a, Vec3 b) { a.m = _mm_mul_ps(a.m, b.m); return a; }
//VM_INLINE Vec3 operator/ (Vec3 a, Vec3 b) { a.m = _mm_div_ps(a.m, b.m); return a; }
//VM_INLINE Vec3 operator* (Vec3 a, float b) { a.m = _mm_mul_ps(a.m, _mm_set1_ps(b)); return a; }
//VM_INLINE Vec3 operator/ (Vec3 a, float b) { a.m = _mm_div_ps(a.m, _mm_set1_ps(b)); return a; }
//VM_INLINE Vec3 operator* (float a, Vec3 b) { b.m = _mm_mul_ps(_mm_set1_ps(a), b.m); return b; }
//VM_INLINE Vec3 operator/ (float a, Vec3 b) { b.m = _mm_div_ps(_mm_set1_ps(a), b.m); return b; }
//VM_INLINE Vec3& operator+= (Vec3& a, Vec3 b) { a = a + b; return a; }
//VM_INLINE Vec3& operator-= (Vec3& a, Vec3 b) { a = a - b; return a; }
//VM_INLINE Vec3& operator*= (Vec3& a, Vec3 b) { a = a * b; return a; }
//VM_INLINE Vec3& operator/= (Vec3& a, Vec3 b) { a = a / b; return a; }
//VM_INLINE Vec3& operator*= (Vec3& a, float b) { a = a * b; return a; }
//VM_INLINE Vec3& operator/= (Vec3& a, float b) { a = a / b; return a; }
//VM_INLINE bool3 operator==(Vec3 a, Vec3 b) { a.m = _mm_cmpeq_ps(a.m, b.m); return a; }
//VM_INLINE bool3 operator!=(Vec3 a, Vec3 b) { a.m = _mm_cmpneq_ps(a.m, b.m); return a; }
//VM_INLINE bool3 operator< (Vec3 a, Vec3 b) { a.m = _mm_cmplt_ps(a.m, b.m); return a; }
//VM_INLINE bool3 operator> (Vec3 a, Vec3 b) { a.m = _mm_cmpgt_ps(a.m, b.m); return a; }
//VM_INLINE bool3 operator<=(Vec3 a, Vec3 b) { a.m = _mm_cmple_ps(a.m, b.m); return a; }
//VM_INLINE bool3 operator>=(Vec3 a, Vec3 b) { a.m = _mm_cmpge_ps(a.m, b.m); return a; }
//VM_INLINE Vec3 min(Vec3 a, Vec3 b) { a.m = _mm_min_ps(a.m, b.m); return a; }
//VM_INLINE Vec3 max(Vec3 a, Vec3 b) { a.m = _mm_max_ps(a.m, b.m); return a; }

VM_INLINE Vec3 operator+ (Vec3 a, Vec3 b) { return Vec3(_mm_add_ps(a.m, b.m)); }
VM_INLINE Vec3 operator- (Vec3 a, Vec3 b) { return Vec3(_mm_sub_ps(a.m, b.m)); }
VM_INLINE Vec3 operator* (Vec3 a, Vec3 b) { return Vec3(_mm_mul_ps(a.m, b.m)); }
VM_INLINE Vec3 operator/ (Vec3 a, Vec3 b) { return Vec3(_mm_div_ps(a.m, b.m)); }
VM_INLINE Vec3 operator* (Vec3 a, float b) { return Vec3(_mm_mul_ps(a.m, _mm_set1_ps(b))); }
VM_INLINE Vec3 operator/ (Vec3 a, float b) { return Vec3(_mm_div_ps(a.m, _mm_set1_ps(b))); }
VM_INLINE Vec3 operator* (float a, Vec3 b) { return Vec3(_mm_mul_ps(_mm_set1_ps(a), b.m)); }
VM_INLINE Vec3 operator/ (float a, Vec3 b) { return Vec3(_mm_div_ps(_mm_set1_ps(a), b.m)); }
VM_INLINE Vec3& operator+= (Vec3& a, Vec3 b) { a = a + b; return a; }
VM_INLINE Vec3& operator-= (Vec3& a, Vec3 b) { a = a - b; return a; }
VM_INLINE Vec3& operator*= (Vec3& a, Vec3 b) { a = a * b; return a; }
VM_INLINE Vec3& operator/= (Vec3& a, Vec3 b) { a = a / b; return a; }
VM_INLINE Vec3& operator*= (Vec3& a, float b) { a = a * b; return a; }
VM_INLINE Vec3& operator/= (Vec3& a, float b) { a = a / b; return a; }
VM_INLINE bool3 operator==(Vec3 a, Vec3 b) { return Vec3(_mm_cmpeq_ps(a.m, b.m)); }
VM_INLINE bool3 operator!=(Vec3 a, Vec3 b) { return Vec3(_mm_cmpneq_ps(a.m, b.m)); }
VM_INLINE bool3 operator< (Vec3 a, Vec3 b) { return Vec3(_mm_cmplt_ps(a.m, b.m)); }
VM_INLINE bool3 operator> (Vec3 a, Vec3 b) { return Vec3(_mm_cmpgt_ps(a.m, b.m)); }
VM_INLINE bool3 operator<=(Vec3 a, Vec3 b) { return Vec3(_mm_cmple_ps(a.m, b.m)); }
VM_INLINE bool3 operator>=(Vec3 a, Vec3 b) { return Vec3(_mm_cmpge_ps(a.m, b.m)); }
VM_INLINE Vec3 min(Vec3 a, Vec3 b) { return Vec3(_mm_min_ps(a.m, b.m)); }
VM_INLINE Vec3 max(Vec3 a, Vec3 b) { return Vec3(_mm_max_ps(a.m, b.m)); }

VM_INLINE Vec3 operator- (Vec3 a) { return Vec3(_mm_xor_ps(a.m, _mm_set1_ps(-0.0))); }

VM_INLINE float hmin(Vec3 v)
{
    v = min(v, SHUFFLE3(v, 1, 0, 2));
    return min(v, SHUFFLE3(v, 2, 0, 1)).getX();
}
VM_INLINE float hmax(Vec3 v)
{
    v = max(v, SHUFFLE3(v, 1, 0, 2));
    return max(v, SHUFFLE3(v, 2, 0, 1)).getX();
}

VM_INLINE Vec3 cross(const Vec3 &a, const Vec3 &b)
{
    // x  <-  a.y*b.z - a.z*b.y
    // y  <-  a.z*b.x - a.x*b.z
    // z  <-  a.x*b.y - a.y*b.x
    // We can save a shuffle by grouping it in this wacky order:
    return (a.zxy() * b - a * b.zxy()).zxy();
}

// Returns a 3-bit code where bit0..bit2 is X..Z
VM_INLINE unsigned mask(Vec3 v) { return _mm_movemask_ps(v.m) & 7; }
// Once we have a comparison, we can branch based on its results:
VM_INLINE bool any(bool3 v) { return mask(v) != 0; }
VM_INLINE bool all(bool3 v) { return mask(v) == 7; }

VM_INLINE Vec3 clamp(Vec3 t, Vec3 a, Vec3 b) { return min(max(t, a), b); }
VM_INLINE float sum(Vec3 v) { return v.getX() + v.getY() + v.getZ(); }
VM_INLINE float dot(Vec3 a, Vec3 b) { return sum(a * b); }

VM_INLINE float length(Vec3 v) { return sqrtf(dot(v, v)); }
VM_INLINE float sqLength(Vec3 v) { return dot(v, v); }
VM_INLINE Vec3 unit_vector(Vec3 v) { return v * (1.0f / length(v)); }
VM_INLINE Vec3 lerp(Vec3 a, Vec3 b, float t) 
{ 
    return (1 - t) * a + t * b;
    //return a + (b - a) * t; 
}


VM_INLINE void assert_unit(Vec3 v)
{
    assert(fabsf(sqLength(v) - 1.0f) < 0.01f);
}

VM_INLINE Vec3 random_in_unit_sphere(__m128i& state4)
{
    Vec3 p;

    do
    {
        p = Vec3(myrand02_x4(state4)) - Vec3(1.0f);

    } while (sqLength(p) >= 1);

    return p;
}


// float8


struct float8
{
    VM_INLINE float8() {}
    VM_INLINE explicit float8(const float* p) { m = _mm256_load_ps(p); }
    VM_INLINE explicit float8(float v) { m = _mm256_set1_ps(v); }
    VM_INLINE explicit float8(__m256 v) { m = v; }

    VM_INLINE void store_aligned(float* p)
    {
        assert((uintptr_t)p % alignof(__m256) == 0);
        _mm256_store_ps(p, m);
    }

    VM_INLINE void store_unaligned(float* p)
    {
        _mm256_storeu_ps(p, m);
    }

    __m256 m;
};


typedef float8 bool8;

VM_INLINE float8 operator+ (float8 a, float8 b) { a.m = _mm256_add_ps(a.m, b.m); return a; }
VM_INLINE float8 operator- (float8 a, float8 b) { a.m = _mm256_sub_ps(a.m, b.m); return a; }
VM_INLINE float8 operator* (float8 a, float8 b) { a.m = _mm256_mul_ps(a.m, b.m); return a; }
VM_INLINE bool8 operator==(float8 a, float8 b) { a.m = _mm256_cmp_ps(a.m, b.m, _CMP_EQ_OQ); return a; }
VM_INLINE bool8 operator!=(float8 a, float8 b) { a.m = _mm256_cmp_ps(a.m, b.m, _CMP_NEQ_UQ); return a; }
VM_INLINE bool8 operator< (float8 a, float8 b) { a.m = _mm256_cmp_ps(a.m, b.m, _CMP_LT_OS); return a; }
VM_INLINE bool8 operator> (float8 a, float8 b) { a.m = _mm256_cmp_ps(a.m, b.m, _CMP_GT_OS); return a; }
VM_INLINE bool8 operator<=(float8 a, float8 b) { a.m = _mm256_cmp_ps(a.m, b.m, _CMP_LE_OS); return a; }
VM_INLINE bool8 operator>=(float8 a, float8 b) { a.m = _mm256_cmp_ps(a.m, b.m, _CMP_GE_OS); return a; }
VM_INLINE bool8 operator&(bool8 a, bool8 b) { a.m = _mm256_and_ps(a.m, b.m); return a; }
VM_INLINE bool8 operator|(bool8 a, bool8 b) { a.m = _mm256_or_ps(a.m, b.m); return a; }
VM_INLINE float8 operator- (float8 a) { a.m = _mm256_xor_ps(a.m, _mm256_set1_ps(-0.0f)); return a; }
VM_INLINE float8 min(float8 a, float8 b) { a.m = _mm256_min_ps(a.m, b.m); return a; }
VM_INLINE float8 max(float8 a, float8 b) { a.m = _mm256_max_ps(a.m, b.m); return a; }

VM_INLINE float8 fma(float8 a, float8 b, float8 add)
{
    return float8(_mm256_fmadd_ps(a.m, b.m, add.m));
}

#else

struct Vec3
{
    float x, y, z;

    Vec3()
    {
    }

    Vec3(float x0, float y0, float z0)
    {
        x = x0;
        y = y0;
        z = z0;
    }

    VM_INLINE float getX() const { return x; }
    VM_INLINE float getY() const { return y; }
    VM_INLINE float getZ() const { return z; }

    VM_INLINE Vec3 & operator+=(const Vec3 & o)
    {
        x += o.x;
        y += o.y;
        z += o.z;
        return *this;
    }

    VM_INLINE Vec3& operator/=(float t)
    {
        float k = 1.0f / t;
        x *= k;
        y *= k;
        z *= k;
        return *this;
    }

    VM_INLINE Vec3& operator*=(float t)
    {
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }
};

VM_INLINE Vec3 operator-(const Vec3 & v)
{
    return Vec3(-v.x, -v.y, -v.z);
}


VM_INLINE Vec3 operator+(const Vec3 & a, const Vec3 & b)
{
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

VM_INLINE Vec3 operator*(const Vec3 & a, const Vec3 & b)
{
    return Vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}

VM_INLINE Vec3 operator*(const Vec3 & a, float t)
{
    return Vec3(a.x * t, a.y * t, a.z * t);
}

VM_INLINE Vec3 operator*(float t, const Vec3 & a)
{
    return Vec3(a.x * t, a.y * t, a.z * t);
}

VM_INLINE Vec3 operator-(const Vec3 & a, const Vec3 & b)
{
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

VM_INLINE float dot(const Vec3 &a, const Vec3 &b)
{
    return a.x* b.x + a.y * b.y + a.z * b.z;
}

VM_INLINE float length(const Vec3 &v)
{
    return sqrtf(dot(v, v));
}

VM_INLINE float sqLength(const Vec3 &v)
{
    return dot(v, v);
}

VM_INLINE Vec3 lerp(const Vec3 &a, const Vec3 &b, float t)
{
    //return (1 - t) * a + t * b;
    return a + (b - a) * t;
}

VM_INLINE Vec3 cross(const Vec3 & a, const Vec3 & b)
{
    return Vec3(a.y * b.z - a.z * b.y, -(a.x * b.z - a.z * b.x), a.x * b.y - a.y * b.x);
}

VM_INLINE Vec3 unit_vector(const Vec3 & v)
{
    return v * (1.0f / length(v));
}

VM_INLINE Vec3 random_in_unit_sphere()
{
    Vec3 p;

    do
    {
        __m128 r = myrand02_x4();
        float* scalars = (float*)& r;
        p = Vec3(scalars[0], scalars[1], scalars[2]) - Vec3(1, 1, 1);

    } while (sqLength(p) >= 1);

    return p;
}

#endif
