#ifndef DSP_UTILSUNIT_H
#define DSP_UTILSUNIT_H

const double Pi = 3.1415926;
const double RAD = Pi / 180.0;

double generateDiscreteSignal(double n, double ts, double A, double f0, double Q = 0.0);
double generateContinuousSignal(double t, double f0, double A = 1.0,  double Q = 0.0);
double generateTrigonometricFunctionSignal(double t, double f0, double N);

#endif //DSP_UTILSUNIT_H
