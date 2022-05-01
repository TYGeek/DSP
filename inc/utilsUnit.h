#ifndef DSP_UTILSUNIT_H
#define DSP_UTILSUNIT_H

#include <vector>
#include <cmath>

const double Pi = 3.1415926;
const double RAD = Pi / 180.0;

using Vec = std::vector<double>;

struct TFrequencyDomainSignalPolar;
struct TFrequencyDomainSignalRectangular;

struct TTimeDomainSignal
{
    explicit TTimeDomainSignal(size_t N);

    size_t N; // number of samples in the time domain
    Vec t;
    Vec x;
};

struct TFrequencyDomainSignalRectangular
{
    explicit TFrequencyDomainSignalRectangular(size_t N);
    explicit TFrequencyDomainSignalRectangular(TFrequencyDomainSignalPolar const& polar);

    size_t N;
    Vec Re;
    Vec Im;
};

struct TFrequencyDomainSignalPolar
{
    explicit TFrequencyDomainSignalPolar(size_t N);
    explicit TFrequencyDomainSignalPolar(TFrequencyDomainSignalRectangular const& rectangular);

    size_t N;
    Vec Mag;
    Vec Phase;
};

TTimeDomainSignal inverseDiscreteFourierTransform(TFrequencyDomainSignalRectangular const& FDS);
TFrequencyDomainSignalRectangular forwardDiscreteFourierTransform(TTimeDomainSignal const& TDS);


static double generateSignal(double t, double A, double w, double Q)
{
    return A * sin(w*t + Q);
}

static double generateTrigonometricFunctionSignal(double t, double w0, double N)
{
    double signal = 0.0;
    for (int k = 1; k <= N; ++k) {
        signal += pow(-1.0, k-1) * 2/k * sin(w0 * k*t);
    }
    return signal;
}

static double ak_function(int k, double w0, Vec const& t_vec, Vec const& f_vec)
{
    double sum = 0.0;
    double T = 2*Pi/w0;

    for(size_t i = 0; i < t_vec.size(); i++)
        sum += f_vec.at(i) * cos(w0 * t_vec.at(i) * k);

    sum *= 2 / T;

    return sum;
}

static double bk_function(double k, double w0, Vec const& t_vec, Vec const& f_vec)
{
    double sum = 0.0;
    double T = 2*Pi/w0;

    for(size_t i = 0; i < t_vec.size(); i++)
        sum += f_vec.at(i) * sin(w0 * t_vec.at(i) * k) ;

    sum *= 2 / T;

    return sum;
}



#endif //DSP_UTILSUNIT_H
