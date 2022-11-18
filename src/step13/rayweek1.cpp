#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif
#endif


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <ctime>
#include <mutex>
#include <thread>
#include <atomic>

#ifndef _countof
#define _countof(a)    (sizeof(a) / sizeof((a)[0]))
#endif


#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <immintrin.h>
#include <chrono>

#define MULTITHREADED   1

#include "../common/common.h"
#include "mymath.h"
#include "soa_sphere.h"

#define SIMD_WIDTH  8

// we can't use clock() in multi-threaded programs
// https://stackoverflow.com/questions/2962785/c-using-clock-to-measure-time-in-multi-threaded-programs
// https://stackoverflow.com/a/19471677

class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return std::chrono::duration_cast<second_>
        (clock_::now() - beg_).count(); }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};



static int tiles_required(int tile_w, int width)
{
    int n = width / tile_w;
    if (n * tile_w < width)
        n++;

    return n;
}


static inline char* align_ptr(char* ptr, size_t align)
{
    return (char*)(((uintptr_t)ptr + (align - 1)) & ~(uintptr_t)(align - 1));
}


class Scene;

struct ThreadData
{
    // in
    Scene* scene;
    Pix* image;
    int image_w;
    int image_h;
    int tile_w_in_pixels;
    int tile_h_in_pixels;
    int samples_per_pixel;

    uint32_t state; // xorshift random state
    __m128i state4; // simd xorshoft random state

    // out
    uint64_t out_num_rays;
};


struct Ray
{
    Ray()
    {
    }

    Ray(const Vec3 &origin, const Vec3 &dir)
    {
        _origin = origin;
        _dir = unit_vector(dir);
    }

    Vec3 point_at_parameter(float t) const
    {
        return _origin + t * _dir;
    }

    Vec3 _origin;
    Vec3 _dir;
};


class Material;

struct HitRecord
{
    float t;
    Vec3 p;
    Vec3 normal;
    Material *material;
};


class Material
{
public:
    virtual ~Material() {}
    virtual bool scatter(const Ray &ray_in, const HitRecord &rec, Vec3 *attenuation, Ray *scattered, ThreadData* td) const = 0;
};


struct Hitable
{
    ~Hitable()
    {
        for (uint32_t i = 0; i < _soa_spheres.getCount(); ++i)
            delete _soa_spheres.getData()->material[i];
        _soa_spheres.cleanUp();

        delete[] _hitables;
        _hitables = nullptr;
    }


    bool hit(const Ray& r, float t_min, float t_max, HitRecord* result_rec) const
    {
        HitRecord temp_rec;

        //
        // spheres
        //

        const float* __restrict cx = _soa_spheres.getData()->center_x;
        const float* __restrict cy = _soa_spheres.getData()->center_y;
        const float* __restrict cz = _soa_spheres.getData()->center_z;
        const float* __restrict inv_radii = _soa_spheres.getData()->inv_radius;
        const float* __restrict radii_sq = _soa_spheres.getData()->radius_sq;

        // We store temporary values on the stack.
        // Some alternatives if we do not want to remove any fix limit:
        // 1. best and fastest: split data into chunks of MAX_SPHERES
        // 2. use alloca() to allocate memory on stack (causes stackoverflow if we have millions or billions of spheres :))
        //    float* __restrict discriminants =  (float*)align_ptr((char*)alloca(sizeof(float) * _soa_spheres.getCount() + 31), 32);
        //    float* __restrict all_nb = (float*)align_ptr((char*)alloca(sizeof(float) * _soa_spheres.getCount() + 31), 32);
        // 3. allocate temporary memory for every thread in ThreadData and use them

        constexpr int MAX_SPHERES = 1024;
        uint32_t num_positives = 0;
        uint32_t positive_idx[MAX_SPHERES];
        alignas(32) float discriminants[MAX_SPHERES];
        alignas(32) float all_nb[MAX_SPHERES];


#if 1 //USE_FLOAT8
        float8 ray_dir_x(r._dir.getX());
        float8 ray_dir_y(r._dir.getY());
        float8 ray_dir_z(r._dir.getZ());
        float8 ray_origin_x(r._origin.getX());
        float8 ray_origin_y(r._origin.getY());
        float8 ray_origin_z(r._origin.getZ());


        for (uint32_t i = 0, N = _soa_spheres.getCount(); i < N; i += SIMD_WIDTH)
        {
            const float8 cox = float8(cx+i) - ray_origin_x;
            const float8 coy = float8(cy+i) - ray_origin_y;
            const float8 coz = float8(cz+i) - ray_origin_z;

            float8 nb = fma(coz, ray_dir_z, fma(coy, ray_dir_y, cox * ray_dir_x));
            nb.store_aligned(all_nb + i);

            float8 c = fma(coz, coz, fma(coy, coy, cox * cox)) - float8(radii_sq + i);
            float8 discr = nb * nb - c;

            discr.store_aligned(discriminants + i);

            int mask = (~_mm256_movemask_ps(discr.m)) & 255;
            //int mask = _mm256_movemask_ps((discr > float8(0.0f)).m);

            if (mask)
            {
                if (mask & (1<<0))
                    positive_idx[num_positives++] = i;
                if (mask & (1<<1))
                    positive_idx[num_positives++] = i + 1;
                if (mask & (1<<2))
                    positive_idx[num_positives++] = i + 2;
                if (mask & (1<<3))
                    positive_idx[num_positives++] = i + 3;
                if (mask & (1<<4))
                    positive_idx[num_positives++] = i + 4;
                if (mask & (1<<5))
                    positive_idx[num_positives++] = i + 5;
                if (mask & (1<<6))
                    positive_idx[num_positives++] = i + 6;
                if (mask & (1<<7))
                    positive_idx[num_positives++] = i + 7;
            }
        }

#else

        const __m256 ray_dir_x = _mm256_set1_ps(r._dir.getX());
        const __m256 ray_dir_y = _mm256_set1_ps(r._dir.getY());
        const __m256 ray_dir_z = _mm256_set1_ps(r._dir.getZ());
        const __m256 ray_origin_x = _mm256_set1_ps(r._origin.getX());
        const __m256 ray_origin_y = _mm256_set1_ps(r._origin.getY());
        const __m256 ray_origin_z = _mm256_set1_ps(r._origin.getZ());


        for (uint32_t i = 0, N = _soa_spheres.getCount(); i < N; i += SIMD_WIDTH)
        {
            const __m256 cox = _mm256_sub_ps(_mm256_load_ps(&cx[i]), ray_origin_x);
            const __m256 coy = _mm256_sub_ps(_mm256_load_ps(&cy[i]), ray_origin_y);
            const __m256 coz = _mm256_sub_ps(_mm256_load_ps(&cz[i]), ray_origin_z);

            __m256 nb = _mm256_mul_ps(cox, ray_dir_x);
            nb = _mm256_fmadd_ps(coy, ray_dir_y, nb);
            nb = _mm256_fmadd_ps(coz, ray_dir_z, nb);
            _mm256_store_ps(&all_nb[i], nb);

            __m256 c = _mm256_mul_ps(cox, cox);
            c = _mm256_fmadd_ps(coy, coy, c);
            c = _mm256_fmadd_ps(coz, coz, c);
            c = _mm256_sub_ps(c, _mm256_load_ps(&radii_sq[i]));

            __m256 discr = _mm256_sub_ps(_mm256_mul_ps(nb, nb), c);
            _mm256_store_ps(&discriminants[i], discr);

            int mask = (~_mm256_movemask_ps(discr)) & 255;

            if (mask)
            {
                if (mask & (1 << 0))
                    positive_idx[num_positives++] = i;
                if (mask & (1 << 1))
                    positive_idx[num_positives++] = i + 1;
                if (mask & (1 << 2))
                    positive_idx[num_positives++] = i + 2;
                if (mask & (1 << 3))
                    positive_idx[num_positives++] = i + 3;
                if (mask & (1 << 4))
                    positive_idx[num_positives++] = i + 4;
                if (mask & (1 << 5))
                    positive_idx[num_positives++] = i + 5;
                if (mask & (1 << 6))
                    positive_idx[num_positives++] = i + 6;
                if (mask & (1 << 7))
                    positive_idx[num_positives++] = i + 7;
            }
        }
#endif

        int hit_index = -1;
        float hit_t;

        for (int i = 0, N = num_positives; i < N; ++i)
        {
            int idx = positive_idx[i];

            float inv_radius = inv_radii[idx];

            // is this a place holder sphere?
            if (inv_radius == 0)
                continue;

            float discr_sq = sqrtf(discriminants[idx]);
            float nb = all_nb[idx];

            float temp = (nb - discr_sq);
            if (temp < t_max && temp > t_min)
            {
                t_max = temp;
                hit_t = temp;
                hit_index = idx;
                continue;
            }

            temp = (nb + discr_sq);
            if (temp < t_max && temp > t_min)
            {
                t_max = temp;
                hit_t = temp;
                hit_index = idx;
                continue;
            }
        }

        if (hit_index != -1)
        {
            result_rec->t = hit_t;
            result_rec->p = r.point_at_parameter(hit_t);
            result_rec->normal = (result_rec->p - Vec3(cx[hit_index], cy[hit_index], cz[hit_index])) * inv_radii[hit_index];
            result_rec->material = _soa_spheres.getData()->material[hit_index];
        }

        //
        // hitable lists
        //

        //for (uint32_t i = 0; i < _num_hitables; ++i)
        //{
        //    if (_hitables[i].hit(r, t_min, closest_so_far, &temp_rec))
        //    {
        //        hit_anything = true;
        //        closest_so_far = temp_rec.t;
        //        *result_rec = temp_rec;
        //    }
        //}

        return hit_index != -1;
    }

    // spheres
    SphereSOA _soa_spheres;

    // meshes, etc.
    // ...

    // list
    Hitable* _hitables = 0;
    uint32_t _num_hitables = 0;
};


Vec3 random_in_unit_disk(uint32_t& state)
{
    Vec3 p;
    do
    {
        p = Vec3((float)myrand02(state), (float)myrand02(state), 0) - Vec3(1, 1, 0);
    } while (sqLength(p) >= 1.0f);

    return p;
}

struct Camera
{
    void init(Vec3 lookfrom, Vec3 lookat, Vec3 vup, float vfov, float aspect, float aperture, float focus_dist)
    {
        _lensRadius = aperture / 2;
        float theta = vfov * (float)M_PI / 180;
        float half_height = tanf(theta / 2);
        float half_width = aspect * half_height;
        _origin = lookfrom;
        _w = unit_vector(lookfrom - lookat);
        _u = unit_vector(cross(vup, _w));
        _v = cross(_w, _u);
        _lowerLeftCorner = _origin - half_width * focus_dist * _u - half_height * focus_dist * _v - focus_dist * _w;
        _horizontal = 2 * half_width * focus_dist * _u;
        _vertical = 2 * half_height * focus_dist * _v;
    }

    VM_INLINE Ray getRay(float s, float t, uint32_t& state)
    {
        Vec3 rd = _lensRadius * random_in_unit_disk(state);
        Vec3 offset = _u * rd.getX() + _v * rd.getY();
        return Ray(_origin + offset, _lowerLeftCorner + s * _horizontal + t * _vertical - _origin - offset);
    }

    Vec3 _origin;
    Vec3 _lowerLeftCorner;
    Vec3 _horizontal;
    Vec3 _vertical;
    Vec3 _u, _v, _w;
    float _lensRadius;
};

class Lambertian : public Material
{
public:
    Lambertian(const Vec3 &a) : albedo(a)
    {
    }

    bool scatter(const Ray &ray_in, const HitRecord &rec, Vec3 *attenuation, Ray *scattered, ThreadData* td) const override
    {
        Vec3 target = rec.p + rec.normal + random_in_unit_sphere(td->state4);
        *scattered = Ray(rec.p, target - rec.p);
        *attenuation = albedo;
        return true;
    }

    Vec3 albedo;
};

Vec3 reflect(const Vec3 &v, const Vec3 &n)
{
    return v - 2 * dot(v, n) * n;
}

class Metal : public Material
{
public:
    Metal(const Vec3 &a, float f) : albedo(a)
    {
        fuzz = f < 1 ? f : 1;
    }

    bool scatter(const Ray &ray_in, const HitRecord &rec, Vec3 *attenuation, Ray *scattered, ThreadData* td) const override
    {
        Vec3 reflected = reflect(ray_in._dir, rec.normal);
        *scattered = Ray(rec.p, reflected + fuzz * random_in_unit_sphere(td->state4));
        *attenuation = albedo;
        return dot(scattered->_dir, rec.normal) > 0;
    }

    Vec3 albedo;
    float fuzz;
};

bool refract(const Vec3 &uv, const Vec3 &n, float ni_over_nt, Vec3 *refracted)
{
    assert_unit(uv);

    float dt = dot(uv, n);
    float discriminant = 1.0f - ni_over_nt * ni_over_nt * (1 - dt * dt);
    if (discriminant > 0)
    {
        *refracted = ni_over_nt * (uv - n * dt) - n * sqrtf(discriminant);
        return true;
    }
    else
        return false;
}

float schlick(float cosine, float ref_idx)
{
    float r0 = (1 - ref_idx) / (1 + ref_idx);
    r0 = r0 * r0;
    return r0 + (1 - r0) * (float)powf((1 - cosine), 5);
}

class Dielectric : public Material
{
    float _refIdx;

public:
    Dielectric(float ri) : _refIdx(ri)
    {
    }

    bool scatter(const Ray &ray_in, const HitRecord &rec, Vec3 *attenuation, Ray *scattered, ThreadData *td) const override
    {
        *attenuation = Vec3(1, 1, 1);

        Vec3 outward_normal;
        Vec3 reflected = reflect(ray_in._dir, rec.normal);
        float ni_over_nt;
        Vec3 refracted;
        float reflect_prob;
        float cosine;

        if (dot(ray_in._dir, rec.normal) > 0)
        {
            outward_normal = -rec.normal;
            ni_over_nt = _refIdx;
            cosine = _refIdx * dot(ray_in._dir, rec.normal);
        }
        else
        {
            outward_normal = rec.normal;
            ni_over_nt = 1.0f / _refIdx;
            cosine = -dot(ray_in._dir, rec.normal);
        }

        if (refract(ray_in._dir, outward_normal, ni_over_nt, &refracted))
        {
            reflect_prob = schlick(cosine, _refIdx);
        }
        else
        {
            reflect_prob = 1;
        }

        if (myrand01(td->state) < reflect_prob)
        {
            *scattered = Ray(rec.p, reflected);
        }
        else
            *scattered = Ray(rec.p, refracted);

        return true;
    }
};


Vec3 color(const Ray& r, Hitable *world, int depth, ThreadData *td)
{
    ++td->out_num_rays;
    HitRecord rec;
    if (world->hit(r, 0.001f, FLT_MAX, &rec))
    {
        Vec3 attenuation;
        Ray scattered;
        if (depth < MAX_BOUNCES && rec.material->scatter(r, rec, &attenuation, &scattered, td))
        {
            return attenuation * color(scattered, world, depth + 1, td);
        }
        else
            return Vec3(0, 0, 0);
    }
    else
    {
        Vec3 unit_direction = r._dir;   // -1..+1
        float t = 0.5f * (unit_direction.getY() + 1.0f);         // 0..1
        return lerp(Vec3(1.0f, 1.0f, 1.0f), Vec3(0.5f, 0.7f, 1.0f), t);
    }
}


class Scene
{
public:
    Camera camera;
    Hitable* hitables = nullptr;

    ~Scene()
    {
        delete hitables;
    }
};


Scene *create_small_scene()
{
    Scene* scene = new Scene;
    Hitable* world = new Hitable;

    scene->hitables = world;

    Vec3 lookfrom(2, 1, 2);
    Vec3 lookat(0, 0, 0);
    float dist_to_focus = 5.0f;
    float aperture = 0.1f;

    scene->camera.init(lookfrom, lookat, Vec3(0, 1, 0), 60, (float)SCREEN_W / (float)SCREEN_H, aperture, dist_to_focus);

    world->_soa_spheres.reserve(5 + SIMD_WIDTH);

    world->_soa_spheres.add(Vec3(0, 0, -1), 0.5f, new Lambertian(Vec3(0.1f, 0.2f, 0.5f)));
    world->_soa_spheres.add(Vec3(0, -100.5f, -1), 100.0f, new Lambertian(Vec3(0.8f, 0.8f, 0)));
    world->_soa_spheres.add(Vec3(1, 0, -1), 0.5f, new Metal(Vec3(0.8f, 0.6f, 0.2f), 0.3f));
    world->_soa_spheres.add(Vec3(-1, 0, -1), 0.5f, new Dielectric(1.5f));
    world->_soa_spheres.add(Vec3(-1, 0, -1), -0.45f, new Dielectric(1.5f));

    // make sure num spheres is multiple of SIMD width
    while (scene->hitables->_soa_spheres.getCount() % SIMD_WIDTH != 0)
        scene->hitables->_soa_spheres.add(Vec3(999999999, 999999999, 999999999), 0, nullptr);

    return scene;
}


Scene* create_medium_scene()
{
    // the aras_p scene :)
    Scene* scene = new Scene;
    Hitable* world = new Hitable;

    scene->hitables = world;

    Vec3 look_from(0, 2, 3);
    Vec3 look_at(0, 0, 0);
    float dist_to_focus = 3;
    float aperture = 0.1f * 0.2f;

    scene->camera.init(look_from, look_at, Vec3(0, 1, 0), 60, (float)SCREEN_W / (float)SCREEN_H, aperture, dist_to_focus);

    world->_soa_spheres.reserve(46 + SIMD_WIDTH);

    world->_soa_spheres.add(Vec3(0, -100.5, -1), 100, new Lambertian(Vec3(0.8f, 0.8f, 0.8f)));
    world->_soa_spheres.add(Vec3(2, 0, -1), 0.5f, new Lambertian(Vec3(0.8f, 0.4f, 0.4f)));
    world->_soa_spheres.add(Vec3(0, 0, -1), 0.5f, new Lambertian(Vec3(0.4f, 0.8f, 0.4f)));
    world->_soa_spheres.add(Vec3(-2, 0, -1), 0.5f, new Metal(Vec3(0.4f, 0.4f, 0.8f), 0));
    world->_soa_spheres.add(Vec3(2, 0, 1), 0.5f, new Metal(Vec3(0.4f, 0.8f, 0.4f), 0));
    world->_soa_spheres.add(Vec3(0, 0, 1), 0.5f, new Metal(Vec3(0.4f, 0.8f, 0.4f), 0.2f));
    world->_soa_spheres.add(Vec3(-2, 0, 1), 0.5f, new Metal(Vec3(0.4f, 0.8f, 0.4f), 0.6f));
    world->_soa_spheres.add(Vec3(0.5f, 1, 0.5f), 0.5f, new Dielectric(1.5f));
    world->_soa_spheres.add(Vec3(-1.5f, 1.5f, 0.f), 0.3f, new Lambertian(Vec3(0.8f, 0.6f, 0.2f)));
    world->_soa_spheres.add(Vec3(4, 0, -3), 0.5f, new Lambertian(Vec3(0.1f, 0.1f, 0.1f)));
    world->_soa_spheres.add(Vec3(3, 0, -3), 0.5f, new Lambertian(Vec3(0.2f, 0.2f, 0.2f)));
    world->_soa_spheres.add(Vec3(2, 0, -3), 0.5f, new Lambertian(Vec3(0.3f, 0.3f, 0.3f)));
    world->_soa_spheres.add(Vec3(1, 0, -3), 0.5f, new Lambertian(Vec3(0.4f, 0.4f, 0.4f)));
    world->_soa_spheres.add(Vec3(0, 0, -3), 0.5f, new Lambertian(Vec3(0.5f, 0.5f, 0.5f)));
    world->_soa_spheres.add(Vec3(-1, 0, -3), 0.5f, new Lambertian(Vec3(0.6f, 0.6f, 0.6f)));
    world->_soa_spheres.add(Vec3(-2, 0, -3), 0.5f, new Lambertian(Vec3(0.7f, 0.7f, 0.7f)));
    world->_soa_spheres.add(Vec3(-3, 0, -3), 0.5f, new Lambertian(Vec3(0.8f, 0.8f, 0.8f)));
    world->_soa_spheres.add(Vec3(-4, 0, -3), 0.5f, new Lambertian(Vec3(0.9f, 0.9f, 0.9f)));
    world->_soa_spheres.add(Vec3(4, 0, -4), 0.5f, new Metal(Vec3(0.1f, 0.1f, 0.1f), 0));
    world->_soa_spheres.add(Vec3(3, 0, -4), 0.5f, new Metal(Vec3(0.2f, 0.2f, 0.2f), 0));
    world->_soa_spheres.add(Vec3(2, 0, -4), 0.5f, new Metal(Vec3(0.3f, 0.3f, 0.3f), 0));
    world->_soa_spheres.add(Vec3(1, 0, -4), 0.5f, new Metal(Vec3(0.4f, 0.4f, 0.4f), 0));
    world->_soa_spheres.add(Vec3(0, 0, -4), 0.5f, new Metal(Vec3(0.5f, 0.5f, 0.5f), 0));
    world->_soa_spheres.add(Vec3(-1, 0, -4), 0.5f, new Metal(Vec3(0.6f, 0.6f, 0.6f), 0));
    world->_soa_spheres.add(Vec3(-2, 0, -4), 0.5f, new Metal(Vec3(0.7f, 0.7f, 0.7f), 0));
    world->_soa_spheres.add(Vec3(-3, 0, -4), 0.5f, new Metal(Vec3(0.8f, 0.8f, 0.8f), 0));
    world->_soa_spheres.add(Vec3(-4, 0, -4), 0.5f, new Metal(Vec3(0.9f, 0.9f, 0.9f), 0));
    world->_soa_spheres.add(Vec3(4, 0, -5), 0.5f, new Metal(Vec3(0.8f, 0.1f, 0.1f), 0));
    world->_soa_spheres.add(Vec3(3, 0, -5), 0.5f, new Metal(Vec3(0.8f, 0.5f, 0.1f), 0));
    world->_soa_spheres.add(Vec3(2, 0, -5), 0.5f, new Metal(Vec3(0.8f, 0.8f, 0.1f), 0));
    world->_soa_spheres.add(Vec3(1, 0, -5), 0.5f, new Metal(Vec3(0.4f, 0.8f, 0.1f), 0));
    world->_soa_spheres.add(Vec3(0, 0, -5), 0.5f, new Metal(Vec3(0.1f, 0.8f, 0.1f), 0));
    world->_soa_spheres.add(Vec3(-1, 0, -5), 0.5f, new Metal(Vec3(0.1f, 0.8f, 0.5f), 0));
    world->_soa_spheres.add(Vec3(-2, 0, -5), 0.5f, new Metal(Vec3(0.1f, 0.8f, 0.8f), 0));
    world->_soa_spheres.add(Vec3(-3, 0, -5), 0.5f, new Metal(Vec3(0.1f, 0.1f, 0.8f), 0));
    world->_soa_spheres.add(Vec3(-4, 0, -5), 0.5f, new Metal(Vec3(0.5f, 0.1f, 0.8f), 0));
    world->_soa_spheres.add(Vec3(4, 0, -6), 0.5f, new Lambertian(Vec3(0.8f, 0.1f, 0.1f)));
    world->_soa_spheres.add(Vec3(3, 0, -6), 0.5f, new Lambertian(Vec3(0.8f, 0.5f, 0.1f)));
    world->_soa_spheres.add(Vec3(2, 0, -6), 0.5f, new Lambertian(Vec3(0.8f, 0.8f, 0.1f)));
    world->_soa_spheres.add(Vec3(1, 0, -6), 0.5f, new Lambertian(Vec3(0.4f, 0.8f, 0.1f)));
    world->_soa_spheres.add(Vec3(0, 0, -6), 0.5f, new Lambertian(Vec3(0.1f, 0.8f, 0.1f)));
    world->_soa_spheres.add(Vec3(-1, 0, -6), 0.5f, new Lambertian(Vec3(0.1f, 0.8f, 0.5f)));
    world->_soa_spheres.add(Vec3(-2, 0, -6), 0.5f, new Lambertian(Vec3(0.1f, 0.8f, 0.8f)));
    world->_soa_spheres.add(Vec3(-3, 0, -6), 0.5f, new Lambertian(Vec3(0.1f, 0.1f, 0.8f)));
    world->_soa_spheres.add(Vec3(-4, 0, -6), 0.5f, new Metal(Vec3(0.5f, 0.1f, 0.8f), 0));
    world->_soa_spheres.add(Vec3(1.5f, 1.5f, -2), 0.3f, new Lambertian(Vec3(0.1f, 0.2f, 0.5f)));

    // make sure num spheres is multiple of SIMD width
    while (scene->hitables->_soa_spheres.getCount() % SIMD_WIDTH != 0)
        scene->hitables->_soa_spheres.add(Vec3(999999999, 999999999, 999999999), 0, nullptr);

    return scene;
}


Scene* create_large_scene()
{
    Scene* scene = new Scene;
    Hitable* world = new Hitable;

    scene->hitables = world;

    Vec3 lookfrom(3, 8, 15);
    Vec3 lookat(0, 0, 0);
    float dist_to_focus = 10.0f;
    float aperture = 0.1f;

    scene->camera.init(lookfrom, lookat, Vec3(0, 1, 0), 60, (float)SCREEN_W / (float)SCREEN_H, aperture, dist_to_focus);

    int W = 30;
    int H = 16;

    world->_soa_spheres.reserve(W * H + 4 + SIMD_WIDTH);

    srand(111);

    for (int y = 0; y < H; ++y)
    {
        for (int x = 0; x < W; ++x)
        {
            Vec3 pix_pos((x - W/2)*1.1f, 0, (y-H/2)*1.1f);

            // CRT random
            float r = (rand() & 0xff) / 255.0f;
            float g = (rand() & 0xff) / 255.0f;
            float b = (rand() & 0xff) / 255.0f;

            int i = x + y*W;
            float radius = 0.45f;
            Material* m;

            if (i % 20 == 0)
            {
                m = new Dielectric(1.2f + i * 0.05f);
            }
            else if (i % 10 == 0)
            {
                m = new Metal(Vec3(r, g, b), 0.01f + 0.5f * y / (float)(H));
                pix_pos += Vec3(0, 0.1f, 0);
            }
            else
            {
                m = new Lambertian(Vec3(r, g, b));
            }

            world->_soa_spheres.add(pix_pos, radius, m);
        }
    }
    //world->_soa_spheres.add(Vec3(0, -100.5f, -1), 100.0f, new Lambertian(Vec3(0.5f, 0.5f, 0.5f)));
    world->_soa_spheres.add(Vec3(0, -1000.5f, 0), 1000, new Lambertian(Vec3(0.5f, 0.5f, 0.5f)));

    world->_soa_spheres.add(Vec3(5, 3, 0), 2, new Metal(Vec3(0.5f, 0.5f, 0.8f), 0.65f));
    world->_soa_spheres.add(Vec3(0, 3, 0), 2, new Dielectric(1.5f));
    world->_soa_spheres.add(Vec3(-5, 3, 0), 2, new Metal(Vec3(0.8f, 0.2f, 0.2f), 0.05f));

    // make sure num spheres is multiple of SIMD width
    while (scene->hitables->_soa_spheres.getCount() % SIMD_WIDTH != 0)
        scene->hitables->_soa_spheres.add(Vec3(999999999, 999999999, 999999999), 0, nullptr);

    return scene;
}


void render_tile(int tile_index, ThreadData *thread_data)
{
    ThreadData td = *thread_data;

    assert(td.image_w >= td.tile_w_in_pixels);

    int num_tiles_x = tiles_required(td.tile_w_in_pixels, td.image_w);

    int tile_x = tile_index % num_tiles_x;
    int tile_y = tile_index / num_tiles_x;

    int y0 = tile_y * td.tile_h_in_pixels;
    int y1 = y0 + td.tile_h_in_pixels;
    int x0 = tile_x * td.tile_w_in_pixels;
    int x1 = x0 + td.tile_w_in_pixels;

    if (x1 > td.image_w)
        x1 = td.image_w;

    if (y1 > td.image_h)
        y1 = td.image_h;

    //printf("TILE %2d img: %p imw: %d imh: %d tw: %d th: %d x: %03d..%03d y: %03d..%03d td: %p\n", tile_index, td.image, td.image_w, td.image_h, td.tile_w_in_pixels, td.tile_h_in_pixels, x0, x1, y0, y1, thread_data);

    Vec3 inv_image_size(1.0f / td.image_w, 1.0f / td.image_h, 0);

    for (int y = y1 - 1; y >= y0; --y)
    {
        Pix* row = &td.image[y * td.image_w];

        for (int x = x0; x < x1; ++x)
        {
            Vec3 col(0, 0, 0);
            Vec3 xy((float)x, (float)y, 0);

            for (int s = 0; s < td.samples_per_pixel; ++s)
            {
                Vec3 uv = (Vec3(myrand01_x4(td.state4)) + xy) * inv_image_size;
                Ray r = td.scene->camera.getRay(uv.getX(), uv.getY(), td.state);

                col += color(r, td.scene->hitables, 0, &td);
            }

            col *= (float)(1.0f / td.samples_per_pixel);

            col = Vec3(sqrtf(col.getX()), sqrtf(col.getY()), sqrtf(col.getZ()));

            int ir = (int)(col.getX() * 255.99f);
            int ig = (int)(col.getY() * 255.99f);
            int ib = (int)(col.getZ() * 255.99f);

            row[x].r = (uint8_t)ir;
            row[x].g = (uint8_t)ig;
            row[x].b = (uint8_t)ib;
        }
    }

    thread_data->out_num_rays = td.out_num_rays;
    thread_data->state = td.state;
    thread_data->state4 = td.state4;
}


class TileRenderScheduler
{
public:
    uint64_t run(uint32_t num_tiles, uint32_t num_threads, ThreadData* thread_data)
    {
        std::thread* threads = new std::thread[num_threads];
        ThreadData* data = new ThreadData[num_threads];

        _num_tiles = num_tiles;
        _next_tile = 0;

        for (uint32_t i = 0; i < num_threads; ++i)
        {
            data[i] = *thread_data;

            // every thread has it's own random state
            data[i].state = 200*i + 10001;
            data[i].state4 = _mm_set_epi32(200 * i + 10001, 200 * i + 10003, 200 * i + 10005, 200 * i + 10007);

            threads[i] = std::thread(thread_function, this, i, &data[i]);
        }

        uint64_t num_rays = 0;

        for (uint32_t i = 0; i < num_threads; ++i)
        {
            threads[i].join();
            num_rays += data[i].out_num_rays;
        }

        delete[] data;
        delete[] threads;

        return num_rays;
    }

private:
    static void thread_function(TileRenderScheduler* scheduler, uint32_t thread_idx, ThreadData* td)
    {
        int tile;

        while ((tile = scheduler->get_tile()) != -1)
            render_tile(tile, td);
    }

    int get_tile()
    {
        int tile = _next_tile.fetch_add(1);

        if (tile < _num_tiles)
            return tile;
        
        return -1;
    }

    std::atomic<int> _next_tile = 0;
    int _num_tiles = 0;    
};


RESULT benchmark(Scene *scene, Pix *pixels, bool write_tga, const char *scene_name)
{
    RESULT result = { 0 };
    Timer timer;
    ThreadData td;

    td.scene = scene;
    td.image = pixels;
    td.image_w = SCREEN_W;
    td.image_h = SCREEN_H;
    td.tile_w_in_pixels = 32; // SCREEN_W;
    td.tile_h_in_pixels = 32; // SCREEN_H;
    td.samples_per_pixel = NUM_SAMPLES_PER_PIXEL;
    td.out_num_rays = 0;

    if (td.tile_w_in_pixels > td.image_w)
        td.tile_w_in_pixels = td.image_w;

    if (td.tile_h_in_pixels > td.image_h)
        td.tile_h_in_pixels = td.image_h;

    int num_tiles_x = tiles_required(td.tile_w_in_pixels, td.image_w);
    int num_tiles_y = tiles_required(td.tile_h_in_pixels, td.image_h);
    int num_tiles = num_tiles_x * num_tiles_y;
    int num_threads = std::thread::hardware_concurrency();

    if (num_tiles < 3 * num_threads && num_threads > 1)
        printf("warning: tile size is too big\n");


#if MULTITHREADED
    TileRenderScheduler scheduler;
    result.num_rays = scheduler.run(num_tiles, num_threads, &td);
#else
    // no threads
    td.state = 10001;
    td.state4 = _mm_set_epi32(1001, 1003, 1005, 1007);

    for (int i = 0; i < num_tiles; ++i)
    {
        render_tile(i, &td);
        result.num_rays += td.out_num_rays;
        td.out_num_rays = 0;
    }
#endif

    result.elapsed_seconds = timer.elapsed();

    uint64_t total_samples = SCREEN_W * SCREEN_H * NUM_SAMPLES_PER_PIXEL;

    printf("%s\n", scene_name);
    printf("elapsed time:   %.3fs\n", result.elapsed_seconds);
    printf("total samples:  %llu\n", (unsigned long long)total_samples);
    printf("total rays:     %llu\n", (unsigned long long)result.num_rays);
    printf("mrays/s:        %0.2f\n", result.get_mrays_per_sec());
    printf("threads:        %d/%u\n", num_threads, std::thread::hardware_concurrency());
    printf("tile size:      %dx%d\n", td.tile_w_in_pixels, td.tile_h_in_pixels);
    printf("\n");


    delete scene;

    if (write_tga)
    {
        char filename[128];
        sprintf(filename, "out_%s.tga", scene_name);

        tga_write_rgb24(filename, SCREEN_W, SCREEN_H, pixels);

//#ifdef _WIN32
//        char cmd[128];
//        sprintf(cmd, "start out_%s.tga", scene_name);
//#else
//        char cmd[128];
//        sprintf(cmd, "xdg-open out_%s.tga", scene_name);
//#endif
//
//        if (system(cmd) != 0)
//            printf("Meh");
    }

    return result;
}


int main(int argc, const char *argv[])
{
#if defined(_DEBUG) && defined(_WIN32)
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

    bool write_tga = false;
    int num_runs = 1;
    const static int MAX_NUMS = 32;
    RESULT results[MAX_NUMS];

    // arguments

    for (int i=1; i<argc; ++i)
    {
        if (strcmp(argv[i], "-w") == 0)
            write_tga = true;

        if (strcmp(argv[i], "-n") == 0 && i + 1 < argc)
        {
            ++i;
            int n = atoi(argv[i]);

            if (n >= 1 && n < MAX_NUMS)
                num_runs = n;
            else
                printf("Invalid num_runs parameter: %d\n", n);
        }
    }

    Pix *pixels = new Pix[SCREEN_W*SCREEN_H];

    memset(pixels, 0, sizeof(SCREEN_W * SCREEN_H * sizeof(pixels[0])));


    // run

    const char *version = "threads";

    for (int i=0; i<num_runs; ++i)
        results[i] = benchmark(create_small_scene(), pixels, write_tga, "small");

    log_results(version, "small", results, num_runs);


    for (int i = 0; i < num_runs; ++i)
        results[i] = benchmark(create_medium_scene(), pixels, write_tga, "medium");

    log_results(version, "medium", results, num_runs);


    for (int i = 0; i < num_runs; ++i)
        results[i] = benchmark(create_large_scene(), pixels, write_tga, "large");

    log_results(version, "large", results, num_runs);


    delete[] pixels;
}
