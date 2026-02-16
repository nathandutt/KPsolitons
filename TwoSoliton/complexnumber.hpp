#include <iostream>
#include <complex>
#include <cmath>

static const double precision=1e-10;
static const double zero_cutoff = -200;
struct ComplexNumber{
	double abs;
	double arg;
    ComplexNumber();
	ComplexNumber(const std::complex<double>& c);
	ComplexNumber(const double abs_, const double arg_);
    bool isZero() const;
	std::complex<double> ToComplex() const;
    double ToDouble() const;
};
ComplexNumber Pow(const ComplexNumber& c, int pow);
ComplexNumber operator*(const ComplexNumber& c1, const ComplexNumber& c2);
ComplexNumber operator+(const ComplexNumber& c1, const ComplexNumber& c2);
ComplexNumber operator-(const ComplexNumber& c1, const ComplexNumber& c2);
ComplexNumber operator/(const ComplexNumber& c1, const ComplexNumber& c2);
bool operator==(const ComplexNumber& c1, const ComplexNumber& c2);
bool operator<(const ComplexNumber& c1, const ComplexNumber& c2);
std::ostream& operator<<(std::ostream& os, const ComplexNumber& c);

