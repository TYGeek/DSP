#include "fourierTransformUnit.h"
#include "complexNumberUnit.h"
#include "utilsUnit.h"
#include <cmath>

// -------------------------------------- TDiscreteTimeDomainSignal ------------------------------------------------
TDiscreteTimeDomainSignal::TDiscreteTimeDomainSignal(size_t N, double fs, const func &fn):
        N{N}, fs{fs}, ts{1/fs}, t(N), x(N)
{
    // generate the N input samples discrete time signal on which we need to perform DFT
    for (int n = 0; n < N; ++n) {
        x.at(n) = fn(n*ts);
    }
}

// ---------------------------------------- TFrequencyDomainSignal -------------------------------------------------
scalarVec TFrequencyDomainSignal::get_Mag() {
    scalarVec vec_Mag;
    vec_Mag.reserve(X.size());

    for (auto& comp : X) {
        vec_Mag.push_back(calc_Mag(comp));
    }

    return vec_Mag;
}

scalarVec TFrequencyDomainSignal::get_Phase() {
    scalarVec vec_Phase;
    vec_Phase.reserve(X.size());

    for (auto& comp : X) {
        vec_Phase.push_back(calc_Phase(comp));
    }

    return vec_Phase;
}

TFrequencyDomainSignal TDiscreteFourierTransform::forwardTransform(const TDiscreteTimeDomainSignal &tds)
{
    TFrequencyDomainSignal fds;

    fds.N =  ( tds.N % 2 == 0 ) ? ( tds.N/2 ) : ( (tds.N+1)/2 ) ; // if N even number N/2, if odd (N+1)/2
    fds.fs = tds.fs;
    fds.ts = tds.ts;
    fds.X.resize(fds.N);
    fds.Freq.resize(fds.N, 0.0);

    // The N separate DFT analysis frequencies [Hz]
    for (int m = 0; m < fds.N; ++m) {
        fds.Freq.at(m) = m * fds.fs / tds.N;
    }


    // determine the DFT of our x(n) input
    for (int k = 0; k < fds.N; ++k) {
        for (int n = 0; n < tds.N; ++n) {
            TComplexNumber rotor(2*Pi * k * n / tds.N );
            fds.X.at(k).Re += tds.x.at(n) * cos( 2*Pi * k * n / tds.N ) ;
            fds.X.at(k).Im -= tds.x.at(n) * sin( 2*Pi * k * n / tds.N ) ;
        }
        // restore magnitude
//        fds.X.at(k).Re *= 1 / sqrt(tds.N);
//        fds.X.at(k).Im *= 1 / sqrt(tds.N);
    }

    return fds;

}

TDiscreteTimeDomainSignal TDiscreteFourierTransform::inverseTransform(const TFrequencyDomainSignal &fds)
{
    TDiscreteTimeDomainSignal tds;
    tds.N =  fds.N*2;
    tds.fs = fds.fs;
    tds.ts = fds.ts;
    tds.t.resize(tds.N);
    tds.x.resize(tds.N);

    // find the cosine and sine wave amplitude
    scalarVec Re(fds.N);
    scalarVec Im(fds.N);
    for (int k = 0; k < fds.N; ++k) {
        Re.at(k) = fds.X.at(k).Re / fds.N;
        Im.at(k) = -fds.X.at(k).Im / fds.N;
    }
    Re.at(0) = Re.at(0) / 2.0;
    Re.at(fds.N-1) = Re.at(fds.N-1) / 2.0;

    // restore time domain signal
    for (int i = 0; i < tds.N; ++i) {
        for (int k = 0; k < fds.N; ++k) {
            tds.x.at(i) += Re.at(k) * cos( 2*Pi * k * i / tds.N ) ;
            tds.x.at(i) += Im.at(k) * sin( 2*Pi * k * i / tds.N ) ;
        }
//        x.at(i) *= 1/ (N);
    }

    return tds;
}
