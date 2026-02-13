#include <stdexcept>
#include <complex> 
#include <iostream>
static const long double equality_precision = 1e-8;
class Monomial{
    private:
        int x;
        int y;
        int t;
        std::complex<long double> value;
    public:
        Monomial(int x_, int y_, int t_, std::complex<long double> val_);
        Monomial(std::complex<long double> val_);
        bool isZero() const;
        friend bool operator<(const Monomial& m1, const Monomial& m2);
        friend Monomial operator+(const Monomial& m1, const Monomial& m2);
        friend Monomial operator-(const Monomial& m1, const Monomial& m2);
        friend Monomial operator*(const Monomial& m1, const Monomial& m2);
        friend bool operator==(const Monomial& m1, const Monomial& m2);
        friend Monomial Derive(const Monomial& m, int idx);
        friend std::ostream& operator<<(std::ostream& os, const Monomial& m);
};
