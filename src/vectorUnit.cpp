#include "vectorUnit.h"
#include <algorithm>
#include <cmath>

TVector::TVector(const TVector &vec) {
    v.reserve(vec.v.size()-1);
    for (auto i : vec.v) {
        v.push_back(i);
    }
}

TVector::TVector(TVector &&vec) noexcept: v{std::move(vec.v)} { }
TVector::TVector(std::vector<float> vec): v{std::move(vec)}{ }

TVector& TVector::operator=(TVector vec) {
    std::swap(v, vec.v);
    return *this;
}

bool TVector::operator==(const TVector &vec) const {
    return std::equal(v.begin(), v.end(), vec.v.begin());
}

bool TVector::operator!=(const TVector &vec) const {
    return not std::equal(v.begin(), v.end(), vec.v.begin());
}

void TVector::zero() {
    std::fill(v.begin(), v.end(), 0.0);
}

TVector TVector::operator-() const {
    TVector vec_negate(*this);
    for (auto& item : vec_negate.v) {
        item *= -1;
    }

    return vec_negate;
}

TVector TVector::operator+(TVector vec) const {

    for(size_t i = 0; i < vec.v.size(); i++)
        vec.v.at(i) += v.at(i);

    return std::move(vec);
}

TVector &TVector::operator+=(TVector const & vec) {
    for (int i = 0; i < vec.v.size(); ++i) {
        v.at(i) += vec.v.at(i);
    }
    return *this;
}

TVector TVector::operator-(TVector const& vec) const {
    TVector res(*this);
    for (int i = 0; i < vec.v.size(); ++i) {
        res.v.at(i) -= vec.v.at(i);
    }
    return res;
}

TVector &TVector::operator-=(const TVector &vec) {
    for (int i = 0; i < vec.v.size(); ++i) {
        v.at(i) -= vec.v.at(i);
    }
    return *this;
}

TVector TVector::operator*(float scalar) const {
    TVector res(*this);
    std::for_each(res.v.begin(), res.v.end(), [&scalar](float& n){ n*=scalar; });
    return res;
}

TVector &TVector::operator*=(float scalar) {
    std::for_each(v.begin(), v.end(), [&scalar](float& n){ n*=scalar; });
    return *this;
}

float TVector::operator*(const TVector &vec) const {
    float sum = 0.0;
    for(size_t i = 0; i < vec.v.size(); i++)
        sum += v.at(i) * vec.v.at(i);

    return sum;
}

void TVector::normalize() {
    float magSq = 0.0;

    for (auto item: v) {
        magSq += item*item;
    }

    if(magSq > 0) {
        float oneOverMag = 1.0f / sqrt(magSq);
        for (auto& item: v) {
            item *= oneOverMag;
        }
    }
}

float vectorMag(const TVector &vec) {
    float magSq = 0.0;

    for (auto item: vec.v) {
        magSq += item*item;
    }

    if(magSq > 0)
    {
        // magSq = magSq / vec.v.size(); // normalize vector to it size
        return sqrt(magSq);
    }
    else
        return 0;
}

float distance(const TVector &a, const TVector &b) {
    TVector res = b - a;
    return vectorMag(res);
}

TVector operator*(float scalar, TVector vec) {
    for(auto& item : vec.v)
        item *= scalar;

    return std::move(vec);
}

float correlation(const TVector &a, const TVector &b) {
    if (a.v.size() != b.v.size())
        return 0.0;

    float dotProduct = a * b;
    float magA = vectorMag(a);
    float magB = vectorMag(b);
    float r = dotProduct / (magA * magB);
    return r;
}
