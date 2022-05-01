#include "matplotlibcpp.h"
#include "utilsUnit.h"
#include <vector>


namespace plt = matplotlibcpp;


int main() {
    double w0 = 2*Pi;
    double t_start = 0.0;
    double t_end = 2*Pi/w0 + t_start;
    double t_step = 0.001;
    int N = static_cast<int>((t_end - t_start) / t_step);

    std::vector<double> t(N);
    std::vector<double> f1(N);
    std::vector<double> f2(N);
    std::vector<double> f3(N);
    std::vector<double> f4(N);

    // create signals with different form
    for(size_t i = 0; i < N; ++i)
    {
        t.at(i) = ( i == 0 ? t_start : t.at(i-1) + t_step );
        f1.at(i) = generateTrigonometricFunctionSignal(t.at(i), w0, 1);
        f2.at(i) = generateTrigonometricFunctionSignal(t.at(i), w0, 3);
        f3.at(i) = generateTrigonometricFunctionSignal(t.at(i), w0, 5);
        f4.at(i) = generateTrigonometricFunctionSignal(t.at(i), w0, 100);
    }

    // Fourier transform
    int K_garmonic = 10;
    std::vector<double> ak_vector(K_garmonic);
    std::vector<double> bk_vector(K_garmonic);

    for (int k = 0; k < K_garmonic; ++k) {
        ak_vector.at(k) = ak_function(k, w0, t, f4);
        bk_vector.at(k) = bk_function(k, w0, t, f4);
    }

    std::vector<double> f_vector(N);

    for (size_t i = 0; i < N; ++i) {
        for (int k = 0; k < K_garmonic; ++k) {
            f_vector.at(i) += ak_vector.at(k) * cos(w0 * k * t.at(i)) +
                                 bk_vector.at(k) * sin(w0 * k * t.at(i));
        }
    }


//    plt::plot(t, f4);
//    plt::plot(t, f2);
//    plt::plot(t, f3);
//    plt::plot(t, f4);
//    plt::plot(t, f_vector);

//    plt::plot(ak_vector);
//    plt::plot(bk_vector);

    plt::grid(true);
    plt::legend();
    plt::show();
}