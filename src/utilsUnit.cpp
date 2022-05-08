#include "utilsUnit.h"
#include <cmath>

double generateDiscreteSignal(double n, double ts, double A, double f0, double Q)
{
    // n - time index integer sequence [1];
    // ts - constant time period between samples [s];
    // A - sine amplitude [1];
    // f0 - frequency [Hz];
    // Q - start phase [rad];

    return A * sin( 2*Pi * f0 * n * ts + Q );
}

double generateContinuousSignal(double t, double f0, double A,  double Q)
{
    // ts - constant time period between samples [s];
    // A - sine amplitude [1];
    // f0 - frequency [Hz];
    // Q - start phase [rad];

    return A * sin( 2*Pi * f0 * t + Q );
}


double generateTrigonometricFunctionSignal(double t, double f0, double N)
{
    double signal = 0.0;
    for (int k = 1; k <= N; ++k) {
        signal += pow(-1.0, k-1) * 2/k * sin( 2*Pi * f0 * k * t );
    }
    return signal;
}
