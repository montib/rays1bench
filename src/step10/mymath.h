#pragma once

#include <math.h>


static uint32_t xorshift_state = 1236787;

inline static uint32_t XorShift32(uint32_t& state)
{
    uint32_t x = state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 15;
    state = x;
    return x;
}

inline float myrand()
{
    return (XorShift32(xorshift_state) & 0xFFFFFF) / 16777216.0f;
}


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

    float getLength() const
    {
        return sqrtf(x * x + y * y + z * z);
    }

    float getSquaredLength() const
    {
        return x * x + y * y + z * z;
    }

    Vec3 getUnitVector() const
    {
        float k = 1.0f / sqrtf(x * x + y * y + z * z);
        return Vec3(x * k, y * k, z * k);
    }

    Vec3 &operator+=(const Vec3 &o)
    {
        x += o.x;
        y += o.y;
        z += o.z;
        return *this;
    }

    Vec3 &operator/=(float t)
    {
        float k = 1.0f / t;
        x *= k;
        y *= k;
        z *= k;
        return *this;
    }

    Vec3 &operator*=(float t)
    {
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }
};

inline Vec3 operator-(const Vec3 &v)
{
    return Vec3(-v.x, -v.y, -v.z);
}


inline Vec3 operator+(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vec3 operator*(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline Vec3 operator*(const Vec3 &a, float t)
{
    return Vec3(a.x * t, a.y * t, a.z * t);
}

inline Vec3 operator*(float t, const Vec3 &a)
{
    return Vec3(a.x * t, a.y * t, a.z * t);
}

inline Vec3 operator-(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline float dot(const Vec3 &a, const Vec3 &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3 &a, const Vec3 &b)
{
    return Vec3(a.y * b.z - a.z * b.y, -(a.x * b.z - a.z * b.x), a.x * b.y - a.y * b.x);
}

inline Vec3 unit_vector(const Vec3 &v)
{
    return v.getUnitVector();
}

inline Vec3 lerp(float t, const Vec3 &a, const Vec3 &b)
{
    return (1 - t) * a + t * b;
}


inline Vec3 random_in_unit_sphere()
{
    Vec3 p;

    do
    {
        p = 2.0f * Vec3((float)myrand(), (float)myrand(), (float)myrand()) - Vec3(1, 1, 1);

    } while (p.getSquaredLength() >= 1);

    return p;
}
