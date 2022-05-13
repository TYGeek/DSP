#include "matplotlibcpp.h"
#include "complexNumberUnit.h"
#include "utilsUnit.h"
#include "fourierTransformUnit.h"
#include <vector>
#include "vectorUnit.h"


namespace plt = matplotlibcpp;

double continuousTimeSignal(double t)
{
    // x(t) = sin( 2*Pi * f0 * t ) + 0.5*sin( 2*Pi * f1 * t + 3*Pi/4 );
    //     A=1          f=1000hz      A2=0.5         f=2000
    return 1*sin( 2*Pi * 1000 * t ) + 0.5*sin( 2*Pi * 2000 * t + 3*Pi/4 );

}


int main() {
    // calculate discrete input signal
    double fs = 8000.0;         // [samples/sec]
    int N = fs/1000;             // samples of discrete signal [1]
    if ( floor( log2(N)) != ceil( log2(N) ) )
        N = pow(2, ceil( log2(N) ));


    TDiscreteTimeDomainSignal TDS(N, fs, continuousTimeSignal);
//    plt::plot(TDS.x);

    TComplexDomainSignal FDS = TFourierTransform::forwardTransform(TDS);
//    plt::plot(FDS.Freq, FDS.get_Mag(), "sb");
//    plt::plot(FDS.Freq, FDS.get_Phase(), ".r");

    TDiscreteTimeDomainSignal FDS_INV = TFourierTransform::inverseTransform(FDS);
//    plt::plot( FDS_INV.x);

    TComplexDomainSignal CDS = TFourierTransform::transform(TComplexDomainSignal(TDS), ETransform::forward);
//    plt::plot(CDS.Freq, CDS.get_Mag(), "sb");

    TComplexDomainSignal CDS_INV = TFourierTransform::transform(CDS, ETransform::inverse);
//    plt::plot( CDS_INV.get_Re());

    plt::grid(true);
    plt::show();
//    s = "<marker><color><line>" ( s = ".r" - red dots, no connecting line )


}