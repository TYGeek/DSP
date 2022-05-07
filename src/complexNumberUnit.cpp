#include "complexNumberUnit.h"
#include <algorithm>
#include <cmath>
#include <iostream>


TComplexNumber::TComplexNumber(double Q):Re{cos(Q)}, Im{sin(Q)} { }
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
    return {Re*scalar, Im*scalar};
}

TComplexNumber &TComplexNumber::operator*=(float scalar) {
    Re *= scalar;
    Im *= scalar;
    return *this;
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
    return atan2(comp.Im , comp.Re);
}
