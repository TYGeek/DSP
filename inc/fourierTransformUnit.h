#ifndef DSP_FOURIERTRANSFORMUNIT_H
#define DSP_FOURIERTRANSFORMUNIT_H

#include <vector>

struct TFrequencyDomainSignal;
struct TDiscreteTimeDomainSignal;
struct TComplexNumber;

using scalarVec = std::vector<double>;
using complexVec  = std::vector<TComplexNumber>;
using func = std::function<double(double)>; // x=f(t)

struct IFourierTransform {
    virtual TFrequencyDomainSignal forwardTransform(TDiscreteTimeDomainSignal const& tds) = 0;
    virtual TDiscreteTimeDomainSignal inverseTransform(TFrequencyDomainSignal const& fds) = 0;
    virtual ~IFourierTransform() = default;
};

struct TDiscreteFourierTransform: IFourierTransform
{
    TFrequencyDomainSignal forwardTransform(TDiscreteTimeDomainSignal const& tds) override;
    TDiscreteTimeDomainSignal inverseTransform(TFrequencyDomainSignal const& fds) override;
};

struct TFastFourierTransform: IFourierTransform
{
    TFrequencyDomainSignal forwardTransform(TDiscreteTimeDomainSignal const& tds) override = 0;
    TDiscreteTimeDomainSignal inverseTransform(TFrequencyDomainSignal const& fds) override = 0;
};


struct TDiscreteTimeDomainSignal
{
    TDiscreteTimeDomainSignal() = default;
    // get discrete signal from continuous x=f(t)
    TDiscreteTimeDomainSignal(size_t N, double fs, func const& fn);

    size_t N;                   // number of samples
    double fs;                  // sampling rate [Hz] at which original signal was sampled
    double ts;                  // time sampling [s]

    scalarVec t;
    scalarVec x;
};

struct TFrequencyDomainSignal
{
    TFrequencyDomainSignal() = default;

    scalarVec get_Mag();
    scalarVec get_Phase();

    size_t N;           // number of samples [1]
    double fs;          // sampling rate at which original signal was sampled [Hz]
    double ts;          // 1/fs time sampling [s]

    complexVec X;       // complex number (Re + jIm)
    scalarVec Freq;     // sampling frequency [Hz]
};

#endif //DSP_FOURIERTRANSFORMUNIT_H
