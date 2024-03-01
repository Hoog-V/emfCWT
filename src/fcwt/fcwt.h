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


#define PI 3.14159265358979323846264338327950288419716939937510582097494459072381640628620899862803482534211706798f
#define sqrt2PI 2.50662827463100050241576528f
#define IPI4 0.75112554446f

#include "fcwt_wavelet.h"
#include "fcwt_scales.h"

class FCWT
{
public:
    FCWT_LIBRARY_API FCWT(Wavelet *pwav, const int pthreads, const bool puse_optimalization_schemes, const bool puse_normalization) : wavelet(pwav),
                                                                                                                                      threads(pthreads),
                                                                                                                                      use_optimalization_schemes(puse_optimalization_schemes),
                                                                                                                                      use_normalization(puse_normalization){};

    void FCWT_LIBRARY_API create_FFT_optimization_plan(const int pmaxsize, const int poptimizationflags);
    void FCWT_LIBRARY_API create_FFT_optimization_plan(const int pmaxsize, const std::string poptimizationflags);
    void FCWT_LIBRARY_API cwt(const float *pinput, const int psize, std::complex<float> *poutput, const Scales *scales);
    void FCWT_LIBRARY_API cwt(const std::complex<float> *pinput, const int psize, std::complex<float> *poutput, const Scales *scales);
    void FCWT_LIBRARY_API cwt(const float *pinput, const int psize, const Scales *scales, std::complex<float> *poutput, const int pn1, const int pn2);
    void FCWT_LIBRARY_API cwt(const std::complex<float> *pinput, const int psize, const Scales *scales, std::complex<float> *poutput, const int pn1, const int pn2);

    Wavelet *wavelet;

private:
    void cwt(const float *pinput, const int psize, std::complex<float> *poutput, const Scales *scales, const bool complexinput);
    void convolve(const fftwf_plan p, fftwf_complex *Ihat, fftwf_complex *O1, std::complex<float> *out, const Wavelet *wav, const int size, const int newsize, const float scale, const bool lastscale);
    void fftbased(const fftwf_plan p, fftwf_complex *Ihat, fftwf_complex *O1, float *out, float const *mother, const int size, const float scale, const bool imaginary, const bool doublesided);
    void fft_normalize(std::complex<float> *out, int size);
    void load_FFT_optimization_plan();
    void daughter_wavelet_multiplication(fftwf_complex *input, fftwf_complex *output, float const *mother, const float scale, const int isize, const bool imaginary, const bool doublesided);

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