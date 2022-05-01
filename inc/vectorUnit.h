#ifndef DSP_VECTORUNIT_H
#define DSP_VECTORUNIT_H

#include <vector>

class TVector
{
public:
    // Constructors
    TVector() = default;
    TVector(TVector const& vec);
    TVector(TVector&& vec) noexcept ;
    explicit TVector(std::vector<float> vec);

    // Assigment operation
    TVector& operator=(TVector vec);
    TVector& operator=(TVector&& vec) noexcept ;

    // Logical operation
    bool operator==(TVector const& vec) const;
    bool operator!=(TVector const& vec) const;

    // Set the vector to zero
    void zero();
    // Unary minus (return negative of the Vector)
    TVector operator-() const;
    // Binary + add operation
    TVector operator+(TVector vec) const;
    TVector& operator+=(TVector const & vec);
    // Binary - subtract operation
    TVector operator-(TVector const& vec) const;
    TVector& operator-=(TVector const& vec);
    // Multiplication by scalar
    TVector operator*(float scalar) const;
    TVector& operator*=(float scalar);
    // Vector dot product
    float operator*(TVector const& vec) const;
    // Normalize the vector
    void normalize();
    // Variables
    std::vector<float> v;
private:

    // Friend function calculate multiplication scalar by vector
    friend TVector operator*(float scalar, TVector vec);
    // Friend function calculate magnitude of a vector
    friend float vectorMag(TVector const& vec);
    // Friend function compute distance between two vectors
    friend float distance(TVector const& a, TVector const& b);
    // Friend function compute correlation between two vectors
    friend float correlation(TVector const& a, TVector const& b);
};

TVector operator*(float scalar, TVector vec);
float vectorMag(TVector const& vec);
float distance(TVector const& a, TVector const& b);
float correlation(TVector const& a, TVector const& b);


#endif //DSP_VECTORUNIT_H
