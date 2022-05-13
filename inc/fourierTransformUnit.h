#ifndef DSP_FOURIERTRANSFORMUNIT_H
#define DSP_FOURIERTRANSFORMUNIT_H

#include <vector>


struct TDiscreteTimeDomainSignal;
struct TComplexNumber;
struct TComplexDomainSignal;

using scalarVec = std::vector<double>;
using complexVec  = std::vector<TComplexNumber>;
using func = std::function<double(double)>; // x=f(t)

enum class ETransform
{
    forward = -1,
    inverse = 1
};

struct TFourierTransform {

    // function to calculate real fourier transform
    static TComplexDomainSignal forwardTransform(TDiscreteTimeDomainSignal const& tds);
    static TDiscreteTimeDomainSignal inverseTransform(TComplexDomainSignal const& fds);
    /**
    * The standard DFT is designed to accept complex input sequences,
    * that is, real inputs have nonzero Re sample values,
    * and the Im sample values are assumed to be zero.
    **/
    static TComplexDomainSignal transform(TComplexDomainSignal const& cds, ETransform direction);

};


struct TDiscreteTimeDomainSignal
{
    TDiscreteTimeDomainSignal();
    TDiscreteTimeDomainSignal(TDiscreteTimeDomainSignal&&) = default;
    // get discrete signal from continuous x=f(t)
    TDiscreteTimeDomainSignal(size_t N, double fs, func const& fn);

    size_t N;                   // number of samples
    double fs;                  // sampling rate [Hz] at which original signal was sampled

    scalarVec t;
    scalarVec x;
};

struct TComplexDomainSignal
{
    TComplexDomainSignal();
    TComplexDomainSignal(TComplexDomainSignal const& cds);
    explicit TComplexDomainSignal(TDiscreteTimeDomainSignal const& tds);
    TComplexDomainSignal(TComplexDomainSignal&& cds) noexcept ;

    scalarVec get_Mag();
    scalarVec get_Phase();
    scalarVec get_Re();
    scalarVec get_Im();

    size_t N;           // number of samples [1]
    double fs;          // sampling rate at which original signal was sampled [Hz]

    complexVec X;       // complex number (Re + jIm)
    scalarVec Freq;     // sampling frequency [Hz]
};

#endif //DSP_FOURIERTRANSFORMUNIT_H
