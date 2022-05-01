#include "signalUnit.h"
#include <cmath>
#include "utilsUnit.h"


TSignal::TSignal(double A, double fs, double Q): A{A},
                                                 w{fs*2*Pi},
                                                 Q{Q},
                                                 T{1/fs},
                                                 fs{fs}
{

}

double TSignal::signal(double t)
{
    return A * sin(w*t + Q);
}
