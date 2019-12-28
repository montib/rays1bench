#pragma once

#ifdef QUICKBENCH

#define MAX_BOUNCES 50
//#define SCREEN_W 640
//#define SCREEN_H 360
#define SCREEN_W 80
#define SCREEN_H 60

#if MULTITHREADED
#define NUM_SAMPLES_PER_PIXEL (4 * 25)
#else
#define NUM_SAMPLES_PER_PIXEL 4
#endif

#else

#define MAX_BOUNCES 50
#define SCREEN_W (1280)
#define SCREEN_H (720)

#if MULTITHREADED
// multi-threaded version is too fast, need more work :)
#define NUM_SAMPLES_PER_PIXEL (10 * 25)
#else
#define NUM_SAMPLES_PER_PIXEL 10
#endif


#endif


#include <stdio.h>

struct RESULT
{
    double elapsed_seconds;
    uint64_t num_rays;

    double get_mrays_per_sec() const
    {        
        return elapsed_seconds ? (num_rays / elapsed_seconds / 1000000.0) : 0;
    }
};

inline void log_results(const char *version, const char *scene, const RESULT *results, int num_runs)
{
    RESULT result = {0};

    // average of results[]

    for (int i=0; i<num_runs; i++)
    {
        result.elapsed_seconds += results[i].elapsed_seconds;
        result.num_rays += results[i].num_rays;
    }

    result.elapsed_seconds /= num_runs;
    result.num_rays /= num_runs;

    // store

    char filename[128];
    sprintf(filename, "out_%s.txt", scene);
    FILE *f = fopen(filename, "wt");

    if (f)
    {
        fprintf(f, "%s|", version);
        fprintf(f, "%.3fs|", result.elapsed_seconds);
        fprintf(f, "%llu|", (long long unsigned)result.num_rays);
        fprintf(f, "%0.3f mrays/s|", result.get_mrays_per_sec());

        fclose(f);
    }
}


struct Pix
{
    uint8_t r, g, b;
};


inline bool tga_write_rgb24(const char* filename, int width, int height, Pix* pixels)   // !!! modifies pixels (swaps R and B)
{
    // https://en.wikipedia.org/wiki/Truevision_TGA

    FILE* f = fopen(filename, "wb");

    if (!f)
        return false;

    uint8_t header[18] = {
        0, // ID length
        0, // no color map
        2, // uncompressed, true color
        0, 0, 0, 0, 0, // Color map specification
        0, 0, 0, 0, // x and y origin
        (uint8_t)(width & 0x00FF),
        (uint8_t)((width & 0xFF00) >> 8),
        (uint8_t)(height & 0x00FF),
        (uint8_t)((height & 0xFF00) >> 8),
        24, // bpp
        0 };

    // swap r, b
    for (int i = 0; i < width * height; ++i)
    {
        uint8_t tmp = pixels[i].r;
        pixels[i].r = pixels[i].b;
        pixels[i].b = tmp;
    }

    fwrite(header, 1, sizeof(header), f);
    fwrite(pixels, 3, width * height, f);

    fclose(f);

    return true;
}