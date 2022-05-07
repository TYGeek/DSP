#ifndef DSP_UTILSUNIT_H
#define DSP_UTILSUNIT_H

#include <vector>
#include <cmath>

const double Pi = 3.1415926;
const double RAD = Pi / 180.0;

using Vec = std::vector<double>;
using func = std::function<double(double)>; // x=f(t)

struct TFrequencyDomainSignalPolar;
struct TFrequencyDomainSignalRectangular;

struct TDiscreteTimeDomainSignal
{
    // get discrete signal from continuous x=f(t)
    TDiscreteTimeDomainSignal(size_t N, double fs, func const& fn);
    // inverse discrete fourier transform
    explicit TDiscreteTimeDomainSignal(TFrequencyDomainSignalRectangular const& fds_rect);

    size_t N;                   // number of samples
    double fs;                  // sampling rate [Hz] at which original signal was sampled
    double ts;                  // time sampling [s]

    Vec t;
    Vec x;
};

struct TFrequencyDomainSignalRectangular
{
    // forward discrete fourier transform
    explicit TFrequencyDomainSignalRectangular(TDiscreteTimeDomainSignal const& tds);
    // convert from polar to rectangular format
    explicit TFrequencyDomainSignalRectangular(TFrequencyDomainSignalPolar const& fds_polar);

    size_t N;       // number of samples [1]
    double fs;      // sampling rate at which original signal was sampled [Hz]
    double ts;      // 1/fs time sampling [s]

    Vec Re;     // [1]
    Vec Im;     // [1]
    Vec Freq;   // sampling frequency [Hz]
private:

};

struct TFrequencyDomainSignalPolar
{
    // convert to polar format
    explicit TFrequencyDomainSignalPolar(TFrequencyDomainSignalRectangular const& fds_rect);

    size_t N;       // number of samples [1]
    double fs;      // sampling rate at which original signal was sampled [Hz]
    double ts;      // time sampling [s]

    Vec Mag;        // [1]
    Vec Phase;      // [rad]
    Vec UWPhase;    // unwrapped phase [rad]
    Vec Freq;       // sampling frequency [Hz]
    Vec Power;      // power spectrum
};

TDiscreteTimeDomainSignal inverseDiscreteFourierTransform(TFrequencyDomainSignalRectangular const& FDS);
TFrequencyDomainSignalRectangular forwardDiscreteFourierTransform(TDiscreteTimeDomainSignal const& TDS);

// n - time index integer sequence [1];
// ts - constant time period between samples [s];
// A - sine amplitude [1];
// f0 - frequency [Hz];
// Q - start phase [rad];
static double generateDiscreteSignal(double n, double ts, double A, double f0, double Q = 0.0)
{
    return A * sin( 2*Pi * f0 * n * ts + Q );
}


// ts - constant time period between samples [s];
// A - sine amplitude [1];
// f0 - frequency [Hz];
// Q - start phase [rad];
static double generateContinuousSignal(double t, double f0, double A = 1.0,  double Q = 0.0)
{
    return A * sin( 2*Pi * f0 * t + Q );
}

static double generateTrigonometricFunctionSignal(double t, double f0, double N)
{
    double signal = 0.0;
    for (int k = 1; k <= N; ++k) {
        signal += pow(-1.0, k-1) * 2/k * sin( 2*Pi * f0 * k * t );
    }
    return signal;
}

#endif //DSP_UTILSUNIT_H
