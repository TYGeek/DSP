#include "utilsUnit.h"
#include <cmath>

// convert between (rectangular -> polar) Frequency Domain Signal representation
TFrequencyDomainSignalPolar rectangularToPolar(TFrequencyDomainSignalRectangular const& rectangular)
{
    TFrequencyDomainSignalPolar polar(rectangular.N);

    for (int k = 0; k < polar.N; ++k) {
        polar.Mag.at(k) = sqrt( pow(rectangular.Re.at(k), 2) +
                                pow(rectangular.Im.at(k), 2) );

        polar.Phase.at(k) = atan2(rectangular.Im.at(k), rectangular.Re.at(k));

        // correct phase
        if ( rectangular.Re.at(k) < 0 and rectangular.Im.at(k) < 0 )
            polar.Phase[k] -= Pi;
        if ( rectangular.Re.at(k) < 0 and rectangular.Im.at(k) >= 0 )
            polar.Phase[k] += Pi;
    }

    return polar;
}

// convert between polar -> rectangular Frequency Domain Signal representation
TFrequencyDomainSignalRectangular polarToRectangular(TFrequencyDomainSignalPolar const& polar)
{
    TFrequencyDomainSignalRectangular rectangular(polar.N);

    for (int k = 0; k < rectangular.N; ++k) {
        rectangular.Re.at(k) = polar.Mag.at(k) * cos(polar.Phase.at(k));
        rectangular.Im.at(k) = polar.Mag.at(k) * sin(polar.Phase.at(k));
    }

    return rectangular;
}

// -------------------------------- TFrequencyDomainSignalRectangular -------------------------------------------------
TFrequencyDomainSignalRectangular::TFrequencyDomainSignalRectangular(size_t N): Re(N), Im(N) { }

TFrequencyDomainSignalRectangular::TFrequencyDomainSignalRectangular(const TFrequencyDomainSignalPolar &polar) {
    *this = polarToRectangular(polar);
}

// ----------------------------------- TFrequencyDomainSignalPolar ----------------------------------------------------
TFrequencyDomainSignalPolar::TFrequencyDomainSignalPolar(size_t N): Mag(N), Phase(N) { }

TFrequencyDomainSignalPolar::TFrequencyDomainSignalPolar(const TFrequencyDomainSignalRectangular &rectangular) {
    *this = rectangularToPolar(rectangular);
}


// -------------------------------------- TTimeDomainSignal -----------------------------------------------------------
TTimeDomainSignal::TTimeDomainSignal(size_t N):N(N), t(N), x(N) { }



// ---------------------------------------- Useful functions ----------------------------------------------------------
TTimeDomainSignal inverseDiscreteFourierTransform(TFrequencyDomainSignalRectangular const& FDS) {

    // find the cosine and sine wave amplitude
    Vec Re(FDS.N);
    Vec Im(FDS.N);

    for (int k = 0; k < FDS.N; ++k) {
        Re.at(k) = FDS.Re.at(k) / FDS.N;
        Im.at(k) = -FDS.Im.at(k) / FDS.N;
    }

    TTimeDomainSignal TDS(FDS.N * 2);

    // restore time domain signal
    for (int k = 0; k < FDS.N; ++k) {
        for (int i = 0; i < TDS.N; ++i) {
            TDS.x.at(i) += Re.at(k) * cos(2*Pi * k * i / TDS.N);
            TDS.x.at(i) += Im.at(k) * sin(2*Pi * k * i / TDS.N);
        }
    }

    return TDS;
}


TFrequencyDomainSignalRectangular forwardDiscreteFourierTransform(TTimeDomainSignal const& TDS)
{
    TFrequencyDomainSignalRectangular FDS(TDS.N/2);

    for (int k = 0; k < FDS.N; ++k) {
        for (int i = 0; i < TDS.N; ++i) {
            FDS.Re.at(k) += TDS.x.at(i) * cos(2*Pi*k*i / TDS.N);
            FDS.Im.at(k) -= TDS.x.at(i) * sin(2*Pi*k*i / TDS.N);
        }
    }

    return FDS;
}