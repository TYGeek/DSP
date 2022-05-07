#include "utilsUnit.h"
#include <cmath>

// ---------------------------------------- Useful functions ----------------------------------------------------------
TDiscreteTimeDomainSignal inverseDiscreteFourierTransform(TFrequencyDomainSignalRectangular const& fds)
{
    return TDiscreteTimeDomainSignal(fds);
}

TFrequencyDomainSignalRectangular forwardDiscreteFourierTransform(TDiscreteTimeDomainSignal const& tds)
{
    return TFrequencyDomainSignalRectangular(tds);
}

// -------------------------------------- TDiscreteTimeDomainSignal ------------------------------------------------
TDiscreteTimeDomainSignal::TDiscreteTimeDomainSignal(size_t N, double fs, const func &fn):
N{N}, fs{fs}, ts{1/fs}, t(N), x(N)
{
    // generate the N input samples discrete time signal on which we need to perform DFT
    for (int n = 0; n < N; ++n) {
        x.at(n) = fn(n*ts);
    }
}

TDiscreteTimeDomainSignal::TDiscreteTimeDomainSignal(const TFrequencyDomainSignalRectangular &fds_rect):
        N{fds_rect.N}, fs{fds_rect.fs}, ts{fds_rect.ts}, t(N), x(N)
{
    // find the cosine and sine wave amplitude
    Vec Re(fds_rect.N);
    Vec Im(fds_rect.N);
    for (int k = 0; k < fds_rect.N; ++k) {
        Re.at(k) = fds_rect.Re.at(k) / fds_rect.N;
        Im.at(k) = -fds_rect.Im.at(k) / fds_rect.N;
    }
    Re.at(0) = Re.at(0) / 2.0;
    Re.at(fds_rect.N-1) = Re.at(fds_rect.N-1) / 2.0;


    // restore time domain signal
    for (int k = 0; k < fds_rect.N; ++k) {
        for (int i = 0; i < N; ++i) {
            x.at(i) += Re.at(k) * cos( 2*Pi * k * i / N ) ;
            x.at(i) += Im.at(k) * sin( 2*Pi * k * i / N ) ;
        }

    }
}

TFrequencyDomainSignalRectangular::TFrequencyDomainSignalRectangular(TDiscreteTimeDomainSignal const& tds):
        N{tds.N}, fs{tds.fs}, ts{tds.ts}, Re(N), Im(N), Freq(N)
{

    // The N separate DFT analysis frequencies [Hz]
    for (int m = 0; m < N; ++m) {
        Freq.at(m) = m * fs / N;
    }

    // determine the DFT of our x(n) input
    for (int k = 0; k < N; ++k) {
        for (int n = 0; n < tds.N; ++n) {
            Re.at(k) += tds.x.at(n) * cos(2 * Pi * k * n / tds.N ) ;
            Im.at(k) -= tds.x.at(n) * sin(2 * Pi * k * n / tds.N ) ;
        }

    }

}

TFrequencyDomainSignalRectangular::TFrequencyDomainSignalRectangular(const TFrequencyDomainSignalPolar &fds_polar):
        N{fds_polar.N}, fs{fds_polar.fs}, ts{fds_polar.ts}, Re(N), Im(N), Freq(N)
{
    // restore Re and Im part
    for (int k = 0; k < N; ++k) {
        Re.at(k) = fds_polar.Mag.at(k) * cos(fds_polar.Phase.at(k));
        Im.at(k) = fds_polar.Mag.at(k) * sin(fds_polar.Phase.at(k));
    }

    // copy frequency
    std::copy(fds_polar.Freq.begin(), fds_polar.Freq.end(), Freq.begin());

}

TFrequencyDomainSignalPolar::TFrequencyDomainSignalPolar(const TFrequencyDomainSignalRectangular &fds_rect):
        N{fds_rect.N}, fs{fds_rect.fs}, ts{fds_rect.ts}, Mag(N), Phase(N), UWPhase(N), Freq(N), Power(N)
{
    // copy frequency
    std::copy(fds_rect.Freq.begin(), fds_rect.Freq.end(), Freq.begin());

    for (int k = 0; k < N; ++k) {
        Mag.at(k) = sqrt( pow(fds_rect.Re.at(k), 2) +
                             pow(fds_rect.Im.at(k), 2) );

        Phase.at(k) = atan2(fds_rect.Im.at(k) , fds_rect.Re.at(k));

        // correct phase
//        if (fds_rect.Re.at(k) < 0 and fds_rect.Im.at(k) < 0 )
//            Phase[k] -= Pi;
//        if (fds_rect.Re.at(k) < 0 and fds_rect.Im.at(k) >= 0 )
//            Phase[k] += Pi;
    }

    // phase unwrapping
    UWPhase.at(0) = 0;
    for (int k = 1; k < N; ++k) {
        int c = static_cast<int>( (UWPhase.at(k - 1) - Phase.at(k)) / (2 * Pi) );
        UWPhase.at(k) = Phase.at(k) + c*2*Pi;
    }

    // power spectrum
    for (int k = 0; k < N; ++k) {
        Power.at(k) = pow(Mag.at(k), 2);
    }
}
