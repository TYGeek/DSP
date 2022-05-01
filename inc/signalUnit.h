#ifndef DSP_SIGNALUNIT_H
#define DSP_SIGNALUNIT_H

class TSignal
{
public:
    TSignal(double A, double fs, double Q);
    double signal(double t);
private:
    double A;   // amplitude [1]
    double w;   // angular rate [rad/s] -> w = 2Pi/T
    double Q;   // phase t=0 [rad]
    double T;   // period [s]
    double fs;  // frequency [Hz] -> fs = 1/T
};

#endif //DSP_SIGNALUNIT_H