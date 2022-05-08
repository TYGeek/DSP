#include "matplotlibcpp.h"
#include "complexNumberUnit.h"
#include "utilsUnit.h"
#include "fourierTransformUnit.h"
#include <vector>


namespace plt = matplotlibcpp;

double continuousTimeSignal(double t)
{
    // x(t) = sin( 2*Pi * f0 * t ) + 0.5*sin( 2*Pi * f1 * t + 3*Pi/4 );
    return sin( 2*Pi * 1000 * t ) + 0.5*sin( 2*Pi * 2000 * t + 3*Pi/4 );
}

int main() {
    // calculate discrete input signal
    double fs = 8000.0;         // [samples/sec]
    double ts = 1/fs;           // [sec]
    int N = 8;                  // samples of discrete signal

    std::unique_ptr<IFourierTransform> transform = std::make_unique<TDiscreteFourierTransform>();
    TDiscreteTimeDomainSignal TDS(N, fs, continuousTimeSignal);
//    plt::plot(TDS.x);
    TFrequencyDomainSignal FDS = transform->forwardTransform(TDS);
    plt::plot(FDS.Freq, FDS.get_Mag(), "sb");
//    plt::plot(FDS.Freq, FDS.get_Phase(), ".r");

    TDiscreteTimeDomainSignal INV = transform->inverseTransform(FDS);
//    plt::plot( INV.x);

    double Q = 90*RAD;
    TComplexNumber p(4, 2); // complex number
//    plt::plot({ p.Re }, { p.Im }, "xb");
    TComplexNumber q(Q);            // rotor
    TComplexNumber p1 = p*q;
//    plt::plot({ p1.Re }, { p1.Im }, "xr");

    plt::grid(true);
    plt::show();

    // s= "<marker><color><line>"
    //    ".r" - red dots, no connecting line
}