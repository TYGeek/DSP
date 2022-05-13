#include "complexNumberUnit.h"
#include <algorithm>
#include <cmath>
#include <iostream>

TComplexNumber::TComplexNumber():Re{0.0}, Im{0.0} { };
TComplexNumber::TComplexNumber(double Q): Re{cos(Q)}, Im{sin(Q)} { }
TComplexNumber::TComplexNumber(double Re, double Im):Re{Re}, Im{Im} { }
TComplexNumber::TComplexNumber(const TComplexNumber &comp) = default;
TComplexNumber::TComplexNumber(TComplexNumber&& comp) noexcept = default;

TComplexNumber &TComplexNumber::operator=(TComplexNumber comp) {
    std::swap(Re, comp.Re);
    std::swap(Im, comp.Im);
    return *this;
}

TComplexNumber TComplexNumber::operator+(TComplexNumber const& comp) const {
    return {Re + comp.Re, Im + comp.Im};
}

TComplexNumber& TComplexNumber::operator+=(const TComplexNumber &comp) {
    Re += comp.Re;
    Im += comp.Im;
    return *this;
}

TComplexNumber TComplexNumber::operator-(const TComplexNumber &comp) const {
    return {Re - comp.Re, Im - comp.Im};
}

TComplexNumber &TComplexNumber::operator-=(const TComplexNumber &comp) {
    Re -= comp.Re;
    Im -= comp.Im;
    return *this;
}

TComplexNumber TComplexNumber::operator*(float scalar) const {
    return {Re*scalar, Im*scalar };
}

TComplexNumber operator*(double scalar, TComplexNumber const& comp) {
    return {scalar * comp.Re, scalar * comp.Im };
}

TComplexNumber &TComplexNumber::operator*=(float scalar) {
    Re *= scalar;
    Im *= scalar;
    return *this;
}

TComplexNumber& TComplexNumber::operator/=(float scalar) {
    Re /= scalar;
    Im /= scalar;
    return *this;
}

TComplexNumber TComplexNumber::operator/(float scalar) const {
    return {Re/scalar, Im/scalar};
}

TComplexNumber TComplexNumber::operator*(const TComplexNumber &comp) const {
    TComplexNumber res(0.0, 0.0);
    res.Re = Re * comp.Re - Im * comp.Im;
    res.Im = Re * comp.Im + Im * comp.Re;
    return res;
}

double calc_Mag(const TComplexNumber &comp) {
    return sqrt(comp.Re*comp.Re + comp.Im*comp.Im);
}

void TComplexNumber::print() const {
    std::cout << Re << " + i" << Im ;
}

TComplexNumber calc_Conjugate(TComplexNumber comp) {
    comp.Im -= comp.Im;
    return std::move(comp);
}

double calc_Phase(const TComplexNumber &comp) {
    // TODO: check restriction below
    // if (Re = 0) and (Im > 0) then return Pi/2;
    // if (Re = 0) and (Im = 0) then return 0;
    // if (Re = 0) and (Im < 0) then return -Pi/2;
    return atan2(comp.Im , comp.Re);
}

double calc_Power(const TComplexNumber &comp) {
    return comp.Re*comp.Re + comp.Im*comp.Im;
}

double calc_Power_dB(const TComplexNumber &comp) {
    // normalized should be
    // return 20 * log10( calc_Power(m) / calc_Power(m)_MAX )
    return 10 * log10(calc_Power(comp));
}




