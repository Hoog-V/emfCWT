#ifndef FCWT_SCALES_H
#define FCWT_SCALES_H

enum SCALETYPE
{
    FCWT_LINSCALES,
    FCWT_LOGSCALES,
    FCWT_LINFREQS
};



class Scales
{
public:
    Scales(){};

    virtual void FCWT_LIBRARY_API getScales(float *pfreqs, const int pnf) = 0;
    virtual void FCWT_LIBRARY_API getFrequencies(float *pfreqs, const int pnf) = 0;

    float *scales;
    int fs;
    float fourwavl;
    int nscales;
};

class ScalesDynamic : public Scales
{
public:
    FCWT_LIBRARY_API ScalesDynamic(Wavelet *pwav, SCALETYPE st, const int fs, const float f0, const float f1, const int fn);

    void FCWT_LIBRARY_API getScales(float *pfreqs, const int pnf);
    void FCWT_LIBRARY_API getFrequencies(float *pfreqs, const int pnf);

private:
    void calculate_logscale_array(const float base, const float four_wavl, const int fs, const float f0, const float f1, const int fn);
    void calculate_linscale_array(const float four_wavl, const int fs, const float f0, const float f1, const int fn);
    void calculate_linfreq_array(const float four_wavl, const int fs, const float f0, const float f1, const int fn);
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
        const float fac = (s1 - s0) / fn;

        for (int i = 0; i < fn; i++)
        {
            scales[i] = (s0 + fac * i);
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

#endif  /* FCWT_SCALES_H */