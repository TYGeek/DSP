#ifndef DSP_COMPLEXNUMBERUNIT_H
#define DSP_COMPLEXNUMBERUNIT_H


struct TComplexNumber {
    TComplexNumber();
    explicit TComplexNumber(double Q);      // complex rotor
    TComplexNumber(double Re, double Im);   // complex number
    TComplexNumber(TComplexNumber const& comp);
    TComplexNumber(TComplexNumber&& comp) noexcept ;

    // Assigment operation (copy-and-swap idiom)
    TComplexNumber& operator=(TComplexNumber comp);
    // Binary + add operation
    TComplexNumber operator+(TComplexNumber const& comp) const;
    TComplexNumber& operator+=(TComplexNumber const & comp);
    // Binary - subtract operation
    TComplexNumber operator-(TComplexNumber const& comp) const;
    TComplexNumber& operator-=(TComplexNumber const& comp);
    // Multiplication by scalar
    TComplexNumber operator*(float scalar) const;
    TComplexNumber& operator*=(float scalar);
    // Divide by scalar
    TComplexNumber operator/(float scalar) const;
    TComplexNumber& operator/=(float scalar);
    // Multiplication by complex number (cross product)
    TComplexNumber operator*(TComplexNumber const& comp) const;
    // print complex number
    void print() const;

    double Re, Im;

private:
    // Friend function calculate magnitude of a complex
    friend double calc_Mag(TComplexNumber const& comp);
    // Friend function calculate phase of a complex
    friend double calc_Phase(TComplexNumber const& comp);
    // Friend function calculate power spectrum
    friend double calc_Power(TComplexNumber const& comp);
    // Friend function calculate power spectrum in dB
    friend double calc_Power_dB(TComplexNumber const& comp);
    // Friend function calculate conjugate of a complex
    friend TComplexNumber calc_Conjugate(TComplexNumber comp);
    // Friend function calculate multiplication by scalar
    friend TComplexNumber operator*(double scalar, TComplexNumber const& comp);
};

TComplexNumber operator+(double scalar, TComplexNumber const& comp);
double calc_Mag(TComplexNumber const& comp);
double calc_Phase(TComplexNumber const& comp);
double calc_Power(TComplexNumber const& comp);
double calc_Power_dB(TComplexNumber const& comp);
TComplexNumber calc_Conjugate(TComplexNumber comp);

#endif //DSP_COMPLEXNUMBERUNIT_H
