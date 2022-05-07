#ifndef DSP_COMPLEXNUMBERUNIT_H
#define DSP_COMPLEXNUMBERUNIT_H


struct TComplexNumber {
    TComplexNumber()=default;
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
    // // Friend function calculate conjugate of a complex
    friend TComplexNumber calc_Conjugate(TComplexNumber comp);
};


double calc_Mag(TComplexNumber const& comp);
double calc_Phase(TComplexNumber const& comp);
TComplexNumber calc_Conjugate(TComplexNumber comp);

#endif //DSP_COMPLEXNUMBERUNIT_H
