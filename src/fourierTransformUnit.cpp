#include "fourierTransformUnit.h"
#include "complexNumberUnit.h"
#include "utilsUnit.h"
#include <cmath>


// -------------------------------------- TDiscreteTimeDomainSignal ------------------------------------------------
TDiscreteTimeDomainSignal::TDiscreteTimeDomainSignal():N{0}, fs{0.0}, t(N), x(N) { }

TDiscreteTimeDomainSignal::TDiscreteTimeDomainSignal(size_t N, double fs, const func &fn):
    N{N}, fs{fs}, t(N), x(N)
{
    // generate the N input samples discrete time signal on which we need to perform DFT
    double ts = 1/fs; // time sampling
    for (int n = 0; n < N; ++n) {
        t.at(n) = n*ts;
        x.at(n) = fn(n*ts);
    }
}

// ---------------------------------------- TComplexDomainSignal -----------------------------------------------------
TComplexDomainSignal::TComplexDomainSignal():N{0}, fs{0.0}, X(N), Freq(N) {}

TComplexDomainSignal::TComplexDomainSignal(const TDiscreteTimeDomainSignal &tds):
    N{tds.N}, fs{tds.fs}, X(N), Freq(N)
{
    for (int i = 0; i < N; ++i) {
        X.at(i).Re = tds.x.at(i);
        X.at(i).Im = 0.0;
    }

    // The N separate analysis frequencies [Hz]
    for (int m = 0; m < N; ++m) {
        Freq.at(m) = m * fs / N;
    }
}

TComplexDomainSignal::TComplexDomainSignal(TComplexDomainSignal const& cds):
    N{cds.N}, fs{cds.fs}, X(cds.X), Freq(cds.Freq)
{

}

TComplexDomainSignal::TComplexDomainSignal(TComplexDomainSignal &&cds) noexcept:
        N{cds.N}, fs{cds.fs}, X(std::move(cds.X)), Freq(std::move(cds.Freq))
{

}

scalarVec TComplexDomainSignal::get_Mag() {
    scalarVec vec_Mag;
    vec_Mag.reserve(X.size());

    for (auto& comp : X) {
        vec_Mag.push_back(calc_Mag(comp));
    }

    return vec_Mag;
}

scalarVec TComplexDomainSignal::get_Phase() {
    scalarVec vec_Phase;
    vec_Phase.reserve(X.size());

    for (auto& comp : X) {
        vec_Phase.push_back(calc_Phase(comp));
    }

    return vec_Phase;
}

scalarVec TComplexDomainSignal::get_Re() {
    scalarVec Re;
    Re.reserve(X.size());

    for (auto& item : X) {
        Re.push_back(item.Re);
    }

    return Re;
}

scalarVec TComplexDomainSignal::get_Im() {
    scalarVec Im;
    Im.reserve(X.size());

    for (auto& item : X) {
        Im.push_back(item.Im);
    }

    return Im;
}

// ---------------------------------------- TFourierTransform --------------------------------------------------------
// TODO: to fix args by const reference instead copy!
complexVec DFT(complexVec x, ETransform trans) // 1) [1, 3]
{
    int N = static_cast<int>( x.size() ); // 1) 2;
    int direction = static_cast<int>(trans);
    complexVec X(N);

    if (trans == ETransform::inverse)
    {
        for (auto& item: x) {
            item /= N;
        }
    }

    // determine the DFT of our x(n) input
    for (int k = 0; k < N; ++k) {
        for (int n = 0; n < N; ++n) {
            TComplexNumber rotor(direction*2*Pi * k*n / N );
            X.at(k) += x.at(n) * rotor;
        }
    }

    return X;
}


complexVec FFT(complexVec const& x, ETransform trans) // 1) [1, 2, 3, 4];  2) [1, 3]
{
    // amount of X(m) complex vector
    int N = static_cast<int>( x.size() ); // 1) [4]; 2) [2];
    int N2 = N/2;
    int direction = static_cast<int>(trans);

    // if FFT_2
    if (N2 == 1) { // 1) False 2) true;
        return DFT(x, trans); // rvalue
    }
    else {
        // separate on even and odd
        complexVec x_even(N2); // 1) [1, 3]
        complexVec x_odd(N2); // 1) [2, 4]
        for (int n = 0; n < N2; ++n) {
            x_even.at(n) = x.at(2 * n);
            x_odd.at(n) = x.at(2 * n + 1);
        }

        // recursion FFT from FFT_N to FFT_2
        complexVec X_even = FFT(x_even, trans); // 1) [C0, C1]
        complexVec X_odd = FFT(x_odd, trans); // 1) [C0 ,C1]

        // calculate W_N_m
        complexVec W_N_m(N2); // 1) N = 4; N2 = 2;
        for (int f0 = 0; f0 < N2; ++f0) {
            W_N_m.at(f0) = TComplexNumber(direction*2*Pi * f0 / N);
        }

        complexVec X(N); // 1) N = 4
        for (int n = 0; n < N2; ++n) {
            // C0 =  {f0 + w0*f2} + w0 * {f1 + f3*w0} = f0 + w0*f2 + w0f1 + w0*w0*f3
            X.at(n) = X_even.at(n) + W_N_m.at(n) * X_odd.at(n);
            // C2 =  {f0 + w0*f2} - w0 * {f1 + f3*w0} = f0 + w0*f2 - w0f1 - w0*w0*f3
            X.at(N2 + n) = X_even.at(n) - W_N_m.at(n) * X_odd.at(n);
        }

        return std::move(X);
    }
}

TComplexDomainSignal TFourierTransform::forwardTransform(const TDiscreteTimeDomainSignal &tds)
{
    TComplexDomainSignal cds;
    cds.N =  ( tds.N % 2 == 0 ) ? ( tds.N/2 ) : ( (tds.N+1)/2 ) ; // if N even number N/2, if odd (N+1)/2
    cds.fs = tds.fs;
    cds.X.resize(cds.N);
    cds.Freq.resize(cds.N, 0.0);

    // The N separate DFT analysis frequencies [Hz]
    for (int m = 0; m < cds.N; ++m) {
        cds.Freq.at(m) = m * cds.fs / tds.N;
    }

    // determine the DFT of our x(n) input
    for (int k = 0; k < cds.N; ++k) {
        for (int n = 0; n < tds.N; ++n) {
            cds.X.at(k).Re += tds.x.at(n) * cos( -2*Pi * k * n / tds.N );
            cds.X.at(k).Im += tds.x.at(n) * sin( -2*Pi * k * n / tds.N );
        }
        // included the normalization factor, 2/N
        cds.X.at(k).Re *= ( 2.0 / tds.N );
        cds.X.at(k).Im *= ( 2.0 / tds.N );
    }
    return cds;
}

TDiscreteTimeDomainSignal TFourierTransform::inverseTransform(const TComplexDomainSignal &cds)
{
    TDiscreteTimeDomainSignal tds;
    tds.N =  cds.N * 2;
    tds.fs = cds.fs;
    tds.t.resize(tds.N);
    tds.x.resize(tds.N);

    // Before using the synthesis equation, the values for Re X[0] and Re X[N/2] must be divided by two.
    scalarVec Re(cds.N);
    scalarVec Im(cds.N);
    for (int k = 0; k < cds.N; ++k) {
        Re.at(k) = cds.X.at(k).Re;
        Im.at(k) = cds.X.at(k).Im;
    }
    Re.at(0) = Re.at(0) / 2;
    Re.at(cds.N - 1) = Re.at(cds.N - 1) / 2;

    // restore time domain signal
    for (int i = 0; i < tds.N; ++i) {
        for (int k = 0; k < cds.N; ++k) {
            tds.x.at(i) += Re.at(k) * cos( 2*Pi * k * i / tds.N ) ;
            tds.x.at(i) -= Im.at(k) * sin( 2*Pi * k * i / tds.N ) ;
        }
    }

    return tds;
}


TComplexDomainSignal TFourierTransform::transform(TComplexDomainSignal const& input, ETransform trans) {
    TComplexDomainSignal output;
    output.N = input.N;
    output.fs = input.fs;
    output.X.resize(output.N);
    output.Freq.resize(output.N, 0.0);

    // The N separate DFT analysis frequencies [Hz]
    for (int m = 0; m < output.N; ++m) {
        output.Freq.at(m) = m * output.fs / input.N;
    }

    output.X = FFT(input.X, trans);

    // included the normalization factor, 2/N
    if( trans == ETransform::forward )
    {
        for (auto& item: output.X)
        {
            item *= ( 2.0f / static_cast<float>(output.N) ) ;
        }
    }


    return output;
}

// ----------------------------------------- Windows -----------------------------------------------------
double rectangularWindow(double n, double N)
{
    return 1;
}

double hanningWindow(double t, double N)
{
    return 0.5 - 0.5 * cos( 2*Pi * t / N );
}

double hammingWindow(double t, double N)
{
    return 0.54 - 0.46 * cos( 2*Pi * t / N );
}