#ifndef FCWT_WAVELET_H
#define FCWT_WAVELET_H


#ifdef _WIN32
#ifdef FCWT_LIBRARY_DLL_BUILDING
#define FCWT_LIBRARY_API __declspec(dllexport)
#else
#if FCWT_LIBRARY_DLL
#define FCWT_LIBRARY_API __declspec(dllimport)
#else /* static or header-only library on Windows */
#define FCWT_LIBRARY_API
#endif
#endif
#else /* Unix */
#define FCWT_LIBRARY_API
#endif

#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <stdbool.h>
#include <vector>
#include <chrono>
#include <cassert>
#include <math.h>
#include <iostream>
#include <sstream>
#include <math.h>
#include <complex>
#include "fftw3.h"


#ifndef SINGLE_THREAD
#include <omp.h>
#endif
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
#include "fftw3.h"
#include <memory>
// check if avx is supported and include the header
#if defined(__AVX__)
#include <immintrin.h>



#define AVX
union U256f
{
    __m256 v;
    float a[8];
};
#endif

class Wavelet
{
public:
    Wavelet(){};
    virtual void generate(float *real, float *imag, const int size, const float scale) { printf("ERROR [generate time complex]: Override this virtual class"); };
    virtual void generate(const int size) { printf("ERROR [generate freq]: Override this virtual class"); };

    constexpr int getSupport(const float scale) { return (int)(fb * scale * 3.0f); };

    virtual void getWavelet(float scale, std::complex<float> *pwav, int pn) { printf("ERROR [getsupport]: Override this virtual class"); };

    float fb;
    int width;
    float four_wavelen;
    bool imag_frequency, doublesided;
    float *mother;
};

template <size_t mother_size>
class MorletStatic : public Wavelet
{
public:
    FCWT_LIBRARY_API MorletStatic(const float bandwidth) : fb(bandwidth)
    { // frequency domain
        mother = mother_static;
        four_wavelen = 0.9876f;
        fb2 = 2.0f * fb * fb;
        ifb = 1.0f / fb;
        imag_frequency = false;
        doublesided = false;
    }

    ~MorletStatic(){};

    void generate(const int size)
    { // frequency domain
        // Frequency domain, because we only need size. Default scale is always 2;
        width = size;

        float tmp1;
        float toradians = (2 * PI) / (float)size;
        float norm = sqrt(2 * PI) * IPI4;

        // calculate array
        for (int w = 0; w < width; w++)
        {
            tmp1 = (2.0f * ((float)w * toradians) * fb - 2.0f * PI * fb);
            tmp1 = -(tmp1 * tmp1) / 2;
            mother[w] = (norm * exp(tmp1));
        }
    }

    void generate(float *real, float *imag, const int size, const float scale)
    {
        float tmp1, tmp2;
        width = getSupport(scale);
        float norm = (float)size * ifb * IPI4;

        // cout << scale << " [";
        for (int t = 0; t < width * 2 + 1; t++)
        {
            tmp1 = (float)(t - width) / scale;
            tmp2 = exp(-(tmp1 * tmp1) / (fb2));

            real[t] = norm * tmp2 * cos(tmp1 * 2.0f * PI) / scale;
            imag[t] = norm * tmp2 * sin(tmp1 * 2.0f * PI) / scale;
            // cout << real[t]*real[t]+imag[t]*imag[t] << ",";
        }
    }

    void getWavelet(const float scale, std::complex<float> *pwav, const int pn)
    {
        int w = getSupport(scale);

        float *real = (float *)malloc(sizeof(float) * std::max(w * 2 + 1, pn));
        float *imag = (float *)malloc(sizeof(float) * std::max(w * 2 + 1, pn));
        for (int t = 0; t < std::max(w * 2 + 1, pn); t++)
        {
            real[t] = 0;
            imag[t] = 0;
        }

        generate(real, imag, pn, scale);

        for (int t = 0; t < pn; t++)
        {
            pwav[t].real(real[t]);
            pwav[t].imag(imag[t]);
        }

        delete real;
        delete imag;
    }

    float fb;

private:
    float ifb, fb2;
    float mother_static[mother_size << 1];
};

class Morlet : public Wavelet
{
public:
    FCWT_LIBRARY_API Morlet(const float bandwidth); // frequency domain
    ~Morlet() { free(mother); };

    void generate(const int size);                                              // frequency domain
    void generate(float *real, float *imag, const int size, const float scale); // time domain
    void getWavelet(float scale, std::complex<float> *pwav, const int pn);
    float fb;

private:
    float ifb, fb2;
};

#endif  /* FCWT_WAVELET_H */