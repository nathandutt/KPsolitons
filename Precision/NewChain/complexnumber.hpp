#include <iostream>
#include <stdio>
#include <complex>
#include <numbers>
static const double significant_digits = 4.;
static const double equality_precision = 1e-8;
static const double PI = std::numbers::pi;

class ComplexNumber{
    private:
        double abs;
        double arg;
    public:
        ComplexNumber();
        ComplexNumber(std::complex<double> c);
        ComplexNumber(double abs, double arg);

        friend std::complex<double> ToComplex(const ComplexNumber& c);
        friend ComplexNumber operator*(const ComplexNumber& c1, const ComplexNumber& c2);
        friend ComplexNumber operator/(const ComplexNumber& c1, const ComplexNumber& c2);
        friend ComplexNumber operator+(const ComplexNumber& c1, const ComplexNumber& c2);
        friend ComplexNumber operator-(const ComplexNumber& c1, const ComplexNumber& c2);
        friend ostream& operator<<(ostream& os, const ComplexNumber& c);

        friend bool operator==(const ComplexNumber& c1, const ComplexNumber& c2);
        friend bool operator<(const ComplexNumber& c1, const ComplexNumber& c2);
        friend bool operator>(const ComplexNumber& c1, const ComplexNumber& c2);
        friend bool Close(const ComplexNumber& c1, const ComplexNumber& c2);
};
