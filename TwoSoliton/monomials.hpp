#include <stdexcept>
#include "complexchain.hpp"
#include <iostream>
static const double equality_precision = 1e-19;
class Monomial{
    public:
        int x;
        int y;
        int t;
        ComplexChain value;
    
        Monomial(int x_, int y_, int t_, ComplexChain val_);
        Monomial(ComplexChain val_);
        ComplexChain Evaluate(const ComplexNumber& x, const ComplexNumber& y, const ComplexNumber& t) const;
        Monomial Simplify() const;
        bool isZero() const;
        friend bool operator<(const Monomial& m1, const Monomial& m2);
        friend Monomial operator+(const Monomial& m1, const Monomial& m2);
        friend Monomial operator-(const Monomial& m1, const Monomial& m2);
        friend Monomial operator*(const Monomial& m1, const Monomial& m2);
        friend bool operator==(const Monomial& m1, const Monomial& m2);
        friend Monomial Derive(const Monomial& m, int idx);
        friend std::ostream& operator<<(std::ostream& os, const Monomial& m);
};
