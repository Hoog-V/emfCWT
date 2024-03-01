//
//  fcwt.h
//  fCWT
//
//  Created by Lukas Arts on 21/12/2020.
//  Copyright © 2021 Lukas Arts.
/*Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
#ifndef FCWT_H
#define FCWT_H

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
#include <complex>

#include <iostream>
#include <sstream>

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

#define PI 3.14159265358979323846264338327950288419716939937510582097494459072381640628620899862803482534211706798f
#define sqrt2PI 2.50662827463100050241576528f
#define IPI4 0.75112554446f

using namespace std;

enum SCALETYPE
{
    FCWT_LINSCALES,
    FCWT_LOGSCALES,
    FCWT_LINFREQS
};

class Wavelet
{
public:
    Wavelet(){};
    virtual void generate(float *real, float *imag, int size, float scale) { printf("ERROR [generate time complex]: Override this virtual class"); };
    virtual void generate(int size) { printf("ERROR [generate freq]: Override this virtual class"); };
    virtual int getSupport(float scale)
    {
        printf("ERROR [getsupport]: Override this virtual class");
        return 0;
    };
    virtual void getWavelet(float scale, complex<float> *pwav, int pn) { printf("ERROR [getsupport]: Override this virtual class"); };

    int width;
    float four_wavelen;
    bool imag_frequency, doublesided;
    float *mother;
};

template <size_t mother_size>
class MorletStatic : public Wavelet
{
public:
    FCWT_LIBRARY_API MorletStatic(float bandwidth)
    { // frequency domain
        mother = mother_static;
        four_wavelen = 0.9876f;
        fb = bandwidth;
        fb2 = 2.0f * fb * fb;
        ifb = 1.0f / fb;
        imag_frequency = false;
        doublesided = false;
    }

    ~MorletStatic(){};

    void generate(int size)
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

    void generate(float *real, float *imag, int size, float scale)
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

    int getSupport(float scale) { return (int)(fb * scale * 3.0f); };
    void getWavelet(float scale, complex<float> *pwav, int pn)
    {
        int w = getSupport(scale);

        float *real = (float *)malloc(sizeof(float) * max(w * 2 + 1, pn));
        float *imag = (float *)malloc(sizeof(float) * max(w * 2 + 1, pn));
        for (int t = 0; t < max(w * 2 + 1, pn); t++)
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
    float mother_static[mother_size];
};

class Morlet : public Wavelet
{
public:
    FCWT_LIBRARY_API Morlet(float bandwidth); // frequency domain
    ~Morlet() { free(mother); };

    void generate(int size);                                        // frequency domain
    void generate(float *real, float *imag, int size, float scale); // time domain
    int getSupport(float scale) { return (int)(fb * scale * 3.0f); };
    void getWavelet(float scale, complex<float> *pwav, int pn);
    float fb;

private:
    float ifb, fb2;
};

class Scales
{
public:
    Scales(){};

    virtual void FCWT_LIBRARY_API getScales(float *pfreqs, int pnf) = 0;
    virtual void FCWT_LIBRARY_API getFrequencies(float *pfreqs, int pnf) = 0;

    float *scales;
    int fs;
    float fourwavl;
    int nscales;

private:
};

class ScalesDynamic : public Scales
{
public:
    FCWT_LIBRARY_API ScalesDynamic(Wavelet *pwav, SCALETYPE st, int fs, float f0, float f1, int fn);

    void FCWT_LIBRARY_API getScales(float *pfreqs, int pnf);
    void FCWT_LIBRARY_API getFrequencies(float *pfreqs, int pnf);

private:
    void calculate_logscale_array(float base, float four_wavl, int fs, float f0, float f1, int fn);
    void calculate_linscale_array(float four_wavl, int fs, float f0, float f1, int fn);
    void calculate_linfreq_array(float four_wavl, int fs, float f0, float f1, int fn);
};

template <size_t scale_size>
class ScalesStatic : public Scales
{
public:
    FCWT_LIBRARY_API ScalesStatic(Wavelet *wav, SCALETYPE st, int afs, float af0, float af1, int afn)
    {

        fs = afs;
        scales = static_scales;
        fourwavl = wav->four_wavelen;
        nscales = afn;

        if (st == SCALETYPE::FCWT_LOGSCALES)
            calculate_logscale_array(2.0f, wav->four_wavelen, afs, af0, af1, afn);
        else if (st == SCALETYPE::FCWT_LINSCALES)
            calculate_linscale_array(wav->four_wavelen, afs, af0, af1, afn);
        else
            calculate_linfreq_array(wav->four_wavelen, afs, af0, af1, afn);
    }
    void FCWT_LIBRARY_API getScales(float *pfreqs, int pnf)
    {
        memcpy(pfreqs, scales, (pnf * sizeof(float)));
    }
    void FCWT_LIBRARY_API getFrequencies(float *pfreqs, int pnf)
    {
        for (int i = 0; i < pnf; i++)
        {
            pfreqs[i] = ((float)fs) / scales[i];
        };
    }

private:
    constexpr void calculate_logscale_array(float base, float four_wavl, int fs, float f0, float f1, int fn)
    {
        // If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;
        float nf0 = f0;
        float nf1 = f1;
        float s0 = (fs / nf1);
        float s1 = (fs / nf0);

        // Cannot pass the nyquist frequency
        assert(("Max frequency cannot be higher than the Nyquist frequency (fs/2)", f1 <= fs / 2));

        float power0 = log(s0) / log(base);
        float power1 = log(s1) / log(base);
        float dpower = power1 - power0;

        for (int i = 0; i < fn; i++)
        {
            float power = power0 + (dpower / (fn - 1)) * i;
            scales[i] = pow(base, power);
        }
    }
    constexpr void calculate_linscale_array(float four_wavl, int fs, float f0, float f1, int fn)
    {
        float nf0 = f0;
        float nf1 = f1;
        // If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;
        float s0 = fs / nf1;
        float s1 = fs / nf0;

        // Cannot pass the nyquist frequency
        assert(("Max frequency cannot be higher than the Nyquist frequency (fs/2)", f1 <= fs / 2));
        float ds = s1 - s0;

        for (int i = 0; i < fn; i++)
        {
            scales[i] = (s0 + (ds / fn) * i);
        }
    }
    constexpr void calculate_linfreq_array(float four_wavl, int fs, float f0, float f1, int fn)
    {
        float nf0 = f0;
        float nf1 = f1;
        // If a signal has fs=100hz and you want to measure [0.1-50]Hz, you need scales 2 to 1000;

        // Cannot pass the nyquist frequency
        assert(("Max frequency cannot be higher than the Nyquist frequency (fs/2)", f1 <= fs / 2));
        float df = nf1 - nf0;

        for (int i = 0; i < fn; i++)
        {
            scales[fn - i - 1] = (((float)fs) / (nf0 + (df / fn) * (float)i));
        }
    }
    float static_scales[scale_size];
};

class FCWT
{
public:
    FCWT_LIBRARY_API FCWT(Wavelet *pwav, int pthreads, bool puse_optimalization_schemes, bool puse_normalization) : wavelet(pwav),
                                                                                                                    threads(pthreads),
                                                                                                                    use_optimalization_schemes(puse_optimalization_schemes),
                                                                                                                    use_normalization(puse_normalization){};

    void FCWT_LIBRARY_API create_FFT_optimization_plan(int pmaxsize, int poptimizationflags);
    void FCWT_LIBRARY_API create_FFT_optimization_plan(int pmaxsize, string poptimizationflags);
    void FCWT_LIBRARY_API cwt(float *pinput, int psize, complex<float> *poutput, Scales *scales);
    void FCWT_LIBRARY_API cwt(complex<float> *pinput, int psize, complex<float> *poutput, Scales *scales);
    void FCWT_LIBRARY_API cwt(float *pinput, int psize, Scales *scales, complex<float> *poutput, int pn1, int pn2);
    void FCWT_LIBRARY_API cwt(complex<float> *pinput, int psize, Scales *scales, complex<float> *poutput, int pn1, int pn2);

    Wavelet *wavelet;

private:
    void cwt(float *pinput, int psize, complex<float> *poutput, Scales *scales, bool complexinput);
    void cwt_static(float *pinput, int psize, float *poutput, float *scales);
    void cwt_dynamic(float *pinput, int psize, float *poutput, float *scales);
    void convolve(fftwf_plan p, fftwf_complex *Ihat, fftwf_complex *O1, complex<float> *out, Wavelet *wav, int size, int newsize, float scale, bool lastscale);
    void convolve(float *in, complex<float> *out, Wavelet *wav, float scale);
    void fftbased(fftwf_plan p, fftwf_complex *Ihat, fftwf_complex *O1, float *out, float *mother, int size, float scale, bool imaginary, bool doublesided);
    void firbased(float *in, float *out, Wavelet *wav, float scale);
    void fft_normalize(complex<float> *out, int size);
    void main(float *Rinput, float *Routput);
    void load_FFT_optimization_plan();
    void daughter_wavelet_multiplication(fftwf_complex *input, fftwf_complex *output, float const *mother, float scale, int isize, bool imaginary, bool doublesided);

    void calculate_logscale_array(float base, float four_wavl, float *scales);
    void calculate_linscale_array(float four_wavl, float *scales);

    int threads;
    int size;
    float fs, f0, f1, fn;
    bool use_optimalization_schemes;
    bool use_normalization;
};

inline int find2power(int n)
{
    int m, m2;

    m = 0;
    m2 = 1 << m; /* 2 to the power of m */
    while (m2 - n < 0)
    {
        m++;
        m2 <<= 1; /* m2 = m2*2 */
    }
    return (m);
}

#endif