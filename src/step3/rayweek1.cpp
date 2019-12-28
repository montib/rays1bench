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

#include "../common/common.h"

static uint32_t xorshift_state = 1236787;

static uint32_t XorShift32(uint32_t &state)
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

struct Ray
{
    Ray()
    {
    }

    Ray(const Vec3 &origin, const Vec3 &dir)
    {
        _origin = origin;
        _dir = dir;
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


class Sphere
{
    Vec3 _center;
    float _radius;
    Material *_material;

public:
    Sphere(const Vec3& center, float radius, Material* material)
    {
        _center = center;
        _radius = radius;
        _material = material;
    }
        
    ~Sphere()
    {
       delete _material;
    }


    bool hit(const Ray &r, float t_min, float t_max, HitRecord *rec) const
    {
        Vec3 oc = r._origin - _center;
        float a = dot(r._dir, r._dir);
        float b = dot(oc, r._dir);
        float c = dot(oc, oc) - _radius * _radius;
        float discriminant = b * b - a * c;
        if (discriminant > 0)
        {
            float temp = (-b - sqrtf(b * b - a * c)) / a;
            if (temp < t_max && temp > t_min)
            {
                rec->t = temp;
                rec->p = r.point_at_parameter(rec->t);
                rec->normal = (rec->p - _center) * (1.0f / _radius);   // TODO / radius
                rec->material = _material;
                return true;
            }

            temp = (-b + sqrtf(b * b - a * c)) / a;
            if (temp < t_max && temp > t_min)
            {
                rec->t = temp;
                rec->p = r.point_at_parameter(rec->t);
                rec->normal = (rec->p - _center) * (1.0f / _radius);   // TODO / radius
                rec->material = _material;
                return true;
            }
        }

        return false;
    }
};


struct Hitable
{
    ~Hitable()
    {
        for (uint32_t i = 0; i < _num_spheres; ++i)
            delete _spheres[i];

        delete[] _spheres;
        _spheres = nullptr;
        _num_spheres = 0;

        delete[] _hitables;
        _hitables = nullptr;
    }

    void alloc_spheres(uint32_t n)
    {
        assert(_spheres == nullptr);
        _spheres = new Sphere*[n];
        memset(_spheres, 0, sizeof(_spheres[0]));
        _num_spheres = n;
    }

    bool hit(const Ray& r, float t_min, float t_max, HitRecord* rec) const
    {
        HitRecord temp_rec;
        bool hit_anything = false;
        float closest_so_far = t_max;

        // spheres
        for (uint32_t i = 0; i < _num_spheres; ++i)
        {
            if (_spheres[i]->hit(r, t_min, closest_so_far, &temp_rec))
            {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                *rec = temp_rec;
            }
        }

        // hitable lists
        for (uint32_t i = 0; i < _num_hitables; ++i)
        {
            if (_hitables[i].hit(r, t_min, closest_so_far, &temp_rec))
            {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                *rec = temp_rec;
            }
        }

        return hit_anything;
    }

    // spheres
    Sphere** _spheres = 0;
    uint32_t _num_spheres = 0;

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
        p = 2.0 * Vec3((float)myrand(), (float)myrand(), 0) - Vec3(1, 1, 0);
    } while (dot(p, p) >= 1.0f);

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
    Ray getRay(float s, float t)
    {
        Vec3 rd = _lensRadius * random_in_unit_disk();
        Vec3 offset = _u * rd.x + _v * rd.y;
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
        Vec3 reflected = reflect(ray_in._dir.getUnitVector(), rec.normal);
        *scattered = Ray(rec.p, reflected + fuzz * random_in_unit_sphere());
        *attenuation = albedo;
        return dot(scattered->_dir, rec.normal) > 0;
    }

    Vec3 albedo;
    float fuzz;
};

bool refract(const Vec3 &v, const Vec3 &n, float ni_over_nt, Vec3 *refracted)
{
    Vec3 uv = unit_vector(v);
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
            cosine = _refIdx * dot(ray_in._dir, rec.normal) / ray_in._dir.getLength();
        }
        else
        {
            outward_normal = rec.normal;
            ni_over_nt = 1.0f / _refIdx;
            cosine = -dot(ray_in._dir, rec.normal) / ray_in._dir.getLength();
        }

        if (refract(ray_in._dir, outward_normal, ni_over_nt, &refracted))
        {
            reflect_prob = schlick(cosine, _refIdx);
        }
        else
        {
            reflect_prob = 1;
        }

        if (myrand() < reflect_prob)
        {
            *scattered = Ray(rec.p, reflected);
        }
        else
            *scattered = Ray(rec.p, refracted);

        return true;
    }
};



float hit_sphere(const Vec3 &center, float radius, const Ray &r) 
{
    Vec3 oc = r._origin - center;
    float a = dot(r._dir, r._dir);
    float b = 2 * dot(oc, r._dir);
    float c = dot(oc, oc) - radius * radius;
    float discriminant = b * b - 4 * a * c;

    if (discriminant < 0)
        return -1;

    return (-b - sqrtf(discriminant)) / (2 * a);
}


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
        Vec3 unit_direction = unit_vector(r._dir);   // -1..+1
        float t = 0.5f * (unit_direction.y + 1.0f);         // 0..1
        return lerp(t, Vec3(1.0f, 1.0f, 1.0f), Vec3(0.5f, 0.7f, 1.0f));
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

    world->alloc_spheres(5);

    world->_spheres[0] = new Sphere(Vec3(0, 0, -1), 0.5f, new Lambertian(Vec3(0.1f, 0.2f, 0.5f)));
    world->_spheres[1] = new Sphere(Vec3(0, -100.5f, -1), 100.0f, new Lambertian(Vec3(0.8f, 0.8f, 0)));
    world->_spheres[2] = new Sphere(Vec3(1, 0, -1), 0.5f, new Metal(Vec3(0.8f, 0.6f, 0.2f), 0.3f));
    world->_spheres[3] = new Sphere(Vec3(-1, 0, -1), 0.5f, new Dielectric(1.5f));
    world->_spheres[4] = new Sphere(Vec3(-1, 0, -1), -0.45f, new Dielectric(1.5f));

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

    world->alloc_spheres(46);

    world->_spheres[0] = new Sphere(Vec3(0, -100.5, -1), 100, new Lambertian(Vec3(0.8f, 0.8f, 0.8f)));
    world->_spheres[1] = new Sphere(Vec3(2, 0, -1), 0.5f, new Lambertian(Vec3(0.8f, 0.4f, 0.4f)));
    world->_spheres[2] = new Sphere(Vec3(0, 0, -1), 0.5f, new Lambertian(Vec3(0.4f, 0.8f, 0.4f)));
    world->_spheres[3] = new Sphere(Vec3(-2, 0, -1), 0.5f, new Metal(Vec3(0.4f, 0.4f, 0.8f), 0));
    world->_spheres[4] = new Sphere(Vec3(2, 0, 1), 0.5f, new Metal(Vec3(0.4f, 0.8f, 0.4f), 0));
    world->_spheres[5] = new Sphere(Vec3(0, 0, 1), 0.5f, new Metal(Vec3(0.4f, 0.8f, 0.4f), 0.2f));
    world->_spheres[6] = new Sphere(Vec3(-2, 0, 1), 0.5f, new Metal(Vec3(0.4f, 0.8f, 0.4f), 0.6f));
    world->_spheres[7] = new Sphere(Vec3(0.5f, 1, 0.5f), 0.5f, new Dielectric(1.5f));
    world->_spheres[8] = new Sphere(Vec3(-1.5f, 1.5f, 0.f), 0.3f, new Lambertian(Vec3(0.8f, 0.6f, 0.2f)));
    world->_spheres[9] = new Sphere(Vec3(4, 0, -3), 0.5f, new Lambertian(Vec3(0.1f, 0.1f, 0.1f)));
    world->_spheres[10] = new Sphere(Vec3(3, 0, -3), 0.5f, new Lambertian(Vec3(0.2f, 0.2f, 0.2f)));
    world->_spheres[11] = new Sphere(Vec3(2, 0, -3), 0.5f, new Lambertian(Vec3(0.3f, 0.3f, 0.3f)));
    world->_spheres[12] = new Sphere(Vec3(1, 0, -3), 0.5f, new Lambertian(Vec3(0.4f, 0.4f, 0.4f)));
    world->_spheres[13] = new Sphere(Vec3(0, 0, -3), 0.5f, new Lambertian(Vec3(0.5f, 0.5f, 0.5f)));
    world->_spheres[14] = new Sphere(Vec3(-1, 0, -3), 0.5f, new Lambertian(Vec3(0.6f, 0.6f, 0.6f)));
    world->_spheres[15] = new Sphere(Vec3(-2, 0, -3), 0.5f, new Lambertian(Vec3(0.7f, 0.7f, 0.7f)));
    world->_spheres[16] = new Sphere(Vec3(-3, 0, -3), 0.5f, new Lambertian(Vec3(0.8f, 0.8f, 0.8f)));
    world->_spheres[17] = new Sphere(Vec3(-4, 0, -3), 0.5f, new Lambertian(Vec3(0.9f, 0.9f, 0.9f)));
    world->_spheres[18] = new Sphere(Vec3(4, 0, -4), 0.5f, new Metal(Vec3(0.1f, 0.1f, 0.1f), 0));
    world->_spheres[19] = new Sphere(Vec3(3, 0, -4), 0.5f, new Metal(Vec3(0.2f, 0.2f, 0.2f), 0));
    world->_spheres[20] = new Sphere(Vec3(2, 0, -4), 0.5f, new Metal(Vec3(0.3f, 0.3f, 0.3f), 0));
    world->_spheres[21] = new Sphere(Vec3(1, 0, -4), 0.5f, new Metal(Vec3(0.4f, 0.4f, 0.4f), 0));
    world->_spheres[22] = new Sphere(Vec3(0, 0, -4), 0.5f, new Metal(Vec3(0.5f, 0.5f, 0.5f), 0));
    world->_spheres[23] = new Sphere(Vec3(-1, 0, -4), 0.5f, new Metal(Vec3(0.6f, 0.6f, 0.6f), 0));
    world->_spheres[24] = new Sphere(Vec3(-2, 0, -4), 0.5f, new Metal(Vec3(0.7f, 0.7f, 0.7f), 0));
    world->_spheres[25] = new Sphere(Vec3(-3, 0, -4), 0.5f, new Metal(Vec3(0.8f, 0.8f, 0.8f), 0));
    world->_spheres[26] = new Sphere(Vec3(-4, 0, -4), 0.5f, new Metal(Vec3(0.9f, 0.9f, 0.9f), 0));
    world->_spheres[27] = new Sphere(Vec3(4, 0, -5), 0.5f, new Metal(Vec3(0.8f, 0.1f, 0.1f), 0));
    world->_spheres[28] = new Sphere(Vec3(3, 0, -5), 0.5f, new Metal(Vec3(0.8f, 0.5f, 0.1f), 0));
    world->_spheres[29] = new Sphere(Vec3(2, 0, -5), 0.5f, new Metal(Vec3(0.8f, 0.8f, 0.1f), 0));
    world->_spheres[30] = new Sphere(Vec3(1, 0, -5), 0.5f, new Metal(Vec3(0.4f, 0.8f, 0.1f), 0));
    world->_spheres[31] = new Sphere(Vec3(0, 0, -5), 0.5f, new Metal(Vec3(0.1f, 0.8f, 0.1f), 0));
    world->_spheres[32] = new Sphere(Vec3(-1, 0, -5), 0.5f, new Metal(Vec3(0.1f, 0.8f, 0.5f), 0));
    world->_spheres[33] = new Sphere(Vec3(-2, 0, -5), 0.5f, new Metal(Vec3(0.1f, 0.8f, 0.8f), 0));
    world->_spheres[34] = new Sphere(Vec3(-3, 0, -5), 0.5f, new Metal(Vec3(0.1f, 0.1f, 0.8f), 0));
    world->_spheres[35] = new Sphere(Vec3(-4, 0, -5), 0.5f, new Metal(Vec3(0.5f, 0.1f, 0.8f), 0));
    world->_spheres[36] = new Sphere(Vec3(4, 0, -6), 0.5f, new Lambertian(Vec3(0.8f, 0.1f, 0.1f)));
    world->_spheres[37] = new Sphere(Vec3(3, 0, -6), 0.5f, new Lambertian(Vec3(0.8f, 0.5f, 0.1f)));
    world->_spheres[38] = new Sphere(Vec3(2, 0, -6), 0.5f, new Lambertian(Vec3(0.8f, 0.8f, 0.1f)));
    world->_spheres[39] = new Sphere(Vec3(1, 0, -6), 0.5f, new Lambertian(Vec3(0.4f, 0.8f, 0.1f)));
    world->_spheres[40] = new Sphere(Vec3(0, 0, -6), 0.5f, new Lambertian(Vec3(0.1f, 0.8f, 0.1f)));
    world->_spheres[41] = new Sphere(Vec3(-1, 0, -6), 0.5f, new Lambertian(Vec3(0.1f, 0.8f, 0.5f)));
    world->_spheres[42] = new Sphere(Vec3(-2, 0, -6), 0.5f, new Lambertian(Vec3(0.1f, 0.8f, 0.8f)));
    world->_spheres[43] = new Sphere(Vec3(-3, 0, -6), 0.5f, new Lambertian(Vec3(0.1f, 0.1f, 0.8f)));
    world->_spheres[44] = new Sphere(Vec3(-4, 0, -6), 0.5f, new Metal(Vec3(0.5f, 0.1f, 0.8f), 0));
    world->_spheres[45] = new Sphere(Vec3(1.5f, 1.5f, -2), 0.3f, new Lambertian(Vec3(0.1f, 0.2f, 0.5f)));
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
    
    world->alloc_spheres(W * H + 4);
    uint32_t idx = 0;

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
                pix_pos.y += 0.1f;
            }
            else
            {
                m = new Lambertian(Vec3(r, g, b));
            }

            world->_spheres[idx++] = new Sphere(pix_pos, radius, m);
        }
    }

    world->_spheres[idx++] = new Sphere(Vec3(0, -1000.5f, 0), 1000, new Lambertian(Vec3(0.5f, 0.5f, 0.5f)));
    world->_spheres[idx++] = new Sphere(Vec3(5, 3, 0), 2, new Metal(Vec3(0.5f, 0.5f, 0.8f), 0.65f));
    world->_spheres[idx++] = new Sphere(Vec3(0, 3, 0), 2, new Dielectric(1.5f));
    world->_spheres[idx++] = new Sphere(Vec3(-5, 3, 0), 2, new Metal(Vec3(0.8f, 0.2f, 0.2f), 0.05f));


    // debug info: how far are the Sphere objects from each other ?

    if ((0))
    {        
        int num_bad = 0;
        long long worst = 0;
        long long average = 0;
        uint8_t *lowest_ptr = (uint8_t*)world->_spheres[0];
        uint8_t *highest_ptr = lowest_ptr;

        for (uint32_t j = 1; j < world->_num_spheres; j++)
        {
            // distance from previous Sphere object
            long long distance = (uint8_t*)world->_spheres[j] - (uint8_t*)world->_spheres[j - 1];

            if (distance != sizeof(Sphere))
            {
                num_bad++;

                average += abs(distance);

                if (abs(distance) > worst)
                    worst = abs(distance);

                printf("%3d. 0x%p distance from previous: %lld\n", j, world->_spheres[j], distance);
            }
                                                   
            // minimum / maximum values of sphere pointers
            if ((uint8_t *)world->_spheres[j] < lowest_ptr)
                lowest_ptr = (uint8_t *)world->_spheres[j];

            if ((uint8_t *)world->_spheres[j] > highest_ptr)
                highest_ptr = (uint8_t *)world->_spheres[j];
        }

        average /= world->_num_spheres;

        printf("sizeof(Sphere): %d, num_bad: %d / %d, worst case: %lld, average: %lld, hi-lo: %lld\n", (int)sizeof(Sphere), num_bad, world->_num_spheres, worst, average, (long long)(highest_ptr - lowest_ptr));
    }

    return scene;
}


RESULT benchmark(Scene *scene, Pix *pixels, bool write_tga, const char *scene_name)
{
    xorshift_state = 1236787;

    RESULT result = { 0 };
    clock_t t0 = std::clock();

    int nx = SCREEN_W;
    int ny = SCREEN_H;

    uint32_t num_rays = 0;

    for (int y = ny - 1; y >= 0; --y)
    {
        for (int x = 0; x < nx; ++x)
        {
            Vec3 col(0, 0, 0);

            for (int s = 0; s < NUM_SAMPLES_PER_PIXEL; ++s)
            {
                float u = (float)((x + myrand()) / nx);
                float v = (float)((y + myrand()) / ny);
                Ray r = scene->camera.getRay(u, v);
                col += color(r, scene->hitables, 0, &num_rays);
            }

            col /= (float)NUM_SAMPLES_PER_PIXEL;

            col = Vec3(sqrtf(col.x), sqrtf(col.y), sqrtf(col.z));

            int ir = (int)(col.x * 255.99f);
            int ig = (int)(col.y * 255.99f);
            int ib = (int)(col.z * 255.99f);
            assert(ir <= 255 && ig <= 255 && ib <= 255);
            pixels[y * nx + x].r = (uint8_t)ir;
            pixels[y * nx + x].g = (uint8_t)ig;
            pixels[y * nx + x].b = (uint8_t)ib;
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
 //    char cmd[128];
//        sprintf(cmd, "start out_%s.tga", scene_name);
 //#else
 //    char cmd[128];
//        sprintf(cmd, "xdg-open out_%s.tga", scene_name);
 //#endif
//
 //    if (system(cmd) != 0)
 //        printf("Meh");
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

    const char *version = "non-virtual hit()";

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
