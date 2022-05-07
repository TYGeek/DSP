#include "matplotlibcpp.h"
#include "utilsUnit.h"
#include <vector>
#include "complexNumberUnit.h"


namespace plt = matplotlibcpp;

double continuousSignal(double t, double f0, double f1)
{
    // x(t) = sin( 2*Pi * f0 * t ) + 0.5*sin( 2*Pi * f1 * t + 3*Pi/4 );
    return generateContinuousSignal(t, f0) + generateContinuousSignal(t,f1,0.5,3*Pi/4);
}

double continuousTimeSignal(double t)
{
    // x(t) = sin( 2*Pi * f0 * t ) + 0.5*sin( 2*Pi * f1 * t + 3*Pi/4 );
    return sin( 2*Pi * 1000 * t ) + 0.5*sin( 2*Pi * 2000 * t + 3*Pi/4 );
}

double DFTRectangular(int m, Vec const& x)
{
    double X_DFT = 0.0;
    size_t N = x.size();
    for(int n = 0; n < N; ++n)
    {
        X_DFT += x.at(n) * cos( 2*Pi * n * m / N );
        X_DFT -= x.at(n) * sin( 2*Pi * n * m / N );
    }
    return X_DFT;
}

int main() {
    // set parameter for continuous input signal
    double f0 = 1000.0;         // [Hz]
    double f1 = 2000.0;         // [Hz]

    // calculate discrete input signal
    double fs = 8000.0;         // [samples/sec]
    double ts = 1/fs;           // [sec]
    int N = 8;                  // samples of discrete signal

    // The N separate DFT analysis frequencies [Hz]
    std::vector<double> f_analyses(N);
    for (int m = 0; m < N; ++m) {
        f_analyses.at(m) = m * fs / N;
    }

    // The N input samples on which we need to perform DFT
    std::vector<double> x(N);
    for (int n = 0; n < N; ++n) {
        x.at(n) = continuousSignal(n*ts, f0 ,f1);
    }

    std::vector<double>X_DFT(N);
    for (int m = 0; m < N; ++m) {
        X_DFT.at(m) = DFTRectangular(m, x);
    }

    TDiscreteTimeDomainSignal TDS(N, fs, continuousTimeSignal);
    TFrequencyDomainSignalRectangular FDS_REC = forwardDiscreteFourierTransform(TDS);
//    plt::plot(FDS_REC.Freq, FDS_REC.Re, "sb");

    TFrequencyDomainSignalPolar FDS_POL(FDS_REC);
//    plt::plot(FDS_POL.Freq, FDS_POL.Mag,  "xb:");
//    plt::plot(FDS_POL.Freq, FDS_POL.Phase, ".r:");

    TDiscreteTimeDomainSignal INV(FDS_REC);
//    plt::plot( INV.x);

    double Q = 90*RAD;
    TComplexNumber p(4, 2); // complex number
    plt::plot({ p.Re }, { p.Im }, "xb");
    TComplexNumber q(Q);            // rotor
    TComplexNumber p1 = p*q;
    plt::plot({ p1.Re }, { p1.Im }, "xr");

    plt::grid(true);
    plt::show();

    // s= "<marker><color><line>"
    //    ".r" - red dots, no connecting line
}