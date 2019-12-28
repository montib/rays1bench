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

#ifndef _countof
#define _countof(a)    (sizeof(a) / sizeof((a)[0]))
#endif


#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <immintrin.h>

#include "../common/common.h"
#include "math.h"
#include "soa_sphere.h"

#define SIMD_WIDTH  4



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
    virtual bool scatter(const Ray &ray_in, const HitRecord &rec, Vec3 *attenuation, Ray *scattered) const = 0;
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
        float* __restrict discriminants = _soa_spheres.getData()->discriminant;
        float* __restrict all_nb = _soa_spheres.getData()->nb;
        
        uint32_t positive_idx[1024];
        uint32_t num_positives = 0;

        const __m128 ray_dir_x = _mm_set1_ps(r._dir.getX());
        const __m128 ray_dir_y = _mm_set1_ps(r._dir.getY());
        const __m128 ray_dir_z = _mm_set1_ps(r._dir.getZ());
        const __m128 ray_origin_x = _mm_set1_ps(r._origin.getX());
        const __m128 ray_origin_y = _mm_set1_ps(r._origin.getY());
        const __m128 ray_origin_z = _mm_set1_ps(r._origin.getZ());


        for (uint32_t i = 0, N = _soa_spheres.getCount(); i < N; i += SIMD_WIDTH)
        {
            const __m128 cox = _mm_sub_ps(_mm_load_ps(&cx[i]), ray_origin_x);
            const __m128 coy = _mm_sub_ps(_mm_load_ps(&cy[i]), ray_origin_y);
            const __m128 coz = _mm_sub_ps(_mm_load_ps(&cz[i]), ray_origin_z);

            __m128 b = _mm_mul_ps(cox, ray_dir_x);
            b = _mm_add_ps(_mm_mul_ps(coy, ray_dir_y), b);
            b = _mm_add_ps(_mm_mul_ps(coz, ray_dir_z), b);
            _mm_storeu_ps(&all_nb[i], b);

            const __m128 radius_sq = _mm_load_ps(&radii_sq[i]);
            __m128 c = _mm_add_ps(_mm_add_ps(_mm_mul_ps(cox, cox), _mm_mul_ps(coy, coy)), _mm_mul_ps(coz, coz));
            c = _mm_sub_ps(c, radius_sq);

            __m128 discr = _mm_sub_ps(_mm_mul_ps(b, b), c);
            _mm_storeu_ps(&discriminants[i], discr);

            int mask = (~_mm_movemask_ps(discr)) & 15;

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
            }
        }

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
            float b = all_nb[idx];

            float temp = (b - discr_sq);
            if (temp < t_max && temp > t_min)
            {
                t_max = temp;
                hit_t = temp;
                hit_index = idx;                
                continue;
            }

            temp = (b + discr_sq);
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


Vec3 random_in_unit_disk()
{
    Vec3 p;
    do
    {
        p = Vec3((float)myrand02(), (float)myrand02(), 0) - Vec3(1, 1, 0);
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

    VM_INLINE Ray getRay(float s, float t)
    {
        Vec3 rd = _lensRadius * random_in_unit_disk();
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

    bool scatter(const Ray &ray_in, const HitRecord &rec, Vec3 *attenuation, Ray *scattered) const override
    {
        Vec3 target = rec.p + rec.normal + random_in_unit_sphere();
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

    bool scatter(const Ray &ray_in, const HitRecord &rec, Vec3 *attenuation, Ray *scattered) const override
    {
        Vec3 reflected = reflect(ray_in._dir, rec.normal);
        *scattered = Ray(rec.p, reflected + fuzz * random_in_unit_sphere());
        *attenuation = albedo;
        return dot(scattered->_dir, rec.normal) > 0;
    }

    Vec3 albedo;
    float fuzz;
};

bool refract(const Vec3 &uv, const Vec3 &n, float ni_over_nt, Vec3 *refracted)
{
    //assert_unit(uv);
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

    bool scatter(const Ray &ray_in, const HitRecord &rec, Vec3 *attenuation, Ray *scattered) const override
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

        if (myrand01() < reflect_prob)
        {
            *scattered = Ray(rec.p, reflected);
        }
        else
            *scattered = Ray(rec.p, refracted);

        return true;
    }
};


Vec3 color(const Ray& r, Hitable *world, int depth, uint32_t *ray_count)
{
    ++(*ray_count);
    HitRecord rec;
    if (world->hit(r, 0.001f, FLT_MAX, &rec))
    {
        Vec3 attenuation;
        Ray scattered;
        if (depth < MAX_BOUNCES && rec.material->scatter(r, rec, &attenuation, &scattered))
        {
            return attenuation * color(scattered, world, depth + 1, ray_count);
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


RESULT benchmark(Scene *scene, Pix *pixels, bool write_tga, const char *scene_name)
{
    RESULT result = { 0 };

    xorshift_state = DEFAULT_RANDOM_STATE;
    xorshift_state4 = DEFAULT_RANDOM4_STATE;


    clock_t t0 = std::clock();

    uint32_t num_rays = 0;

    for (int y = SCREEN_H - 1; y >= 0; --y)
    {
        for (int x = 0; x < SCREEN_W; ++x)
        {
            Vec3 col(0, 0, 0);
            Vec3 xy((float)x, (float)y, 0);

            for (int s = 0; s < NUM_SAMPLES_PER_PIXEL; ++s)
            {
                Vec3 uv = (Vec3(myrand01_x4()) + xy) * Vec3(1.0f / SCREEN_W, 1.0f / SCREEN_H, 0);
                Ray r = scene->camera.getRay(uv.getX(), uv.getY());

                col += color(r, scene->hitables, 0, &num_rays);
            }

            col /= (float)NUM_SAMPLES_PER_PIXEL;

            col = Vec3(sqrtf(col.getX()), sqrtf(col.getY()), sqrtf(col.getZ()));

            int ir = (int)(col.getX() * 255.99f);
            int ig = (int)(col.getY() * 255.99f);
            int ib = (int)(col.getZ() * 255.99f);
            assert(ir <= 255 && ig <= 255 && ib <= 255);
            pixels[y * SCREEN_W + x].r = (uint8_t)ir;
            pixels[y * SCREEN_W + x].g = (uint8_t)ig;
            pixels[y * SCREEN_W + x].b = (uint8_t)ib;
        }
    }

    clock_t t1 = std::clock();
    result.elapsed_seconds = (t1 - t0) / (double)CLOCKS_PER_SEC;
    result.num_rays = num_rays;

    uint64_t total_samples = SCREEN_W * SCREEN_H * NUM_SAMPLES_PER_PIXEL;

    printf("%s\n", scene_name);
    printf("elapsed time:   %.3fs\n", result.elapsed_seconds);
    printf("total samples:  %llu\n", (unsigned long long)total_samples);
    printf("total rays:     %llu\n", (unsigned long long)result.num_rays);
    printf("mrays/s:        %0.2f\n", result.get_mrays_per_sec());

    printf("\n");
    

    delete scene;
    
    if (write_tga)
    {
    char filename[128];
        sprintf(filename, "out_%s.tga", scene_name);

    tga_write_rgb24(filename, SCREEN_W, SCREEN_H, pixels);
    
//#ifdef _WIN32
//     char cmd[128];
//        sprintf(cmd, "start out_%s.tga", scene_name);
//#else
//     char cmd[128];
//        sprintf(cmd, "xdg-open out_%s.tga", scene_name);
//#endif
//
//     if (system(cmd) != 0)
//         printf("Meh");
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

    const char *version = "SIMD Vec3, randoms";

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
