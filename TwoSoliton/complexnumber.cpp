#include "complexnumber.hpp"

double wrap(const double arg){
    auto a = arg + M_PI;
    a = std::fmod(a, 2*M_PI);
    if(a<0.) a = a+2*M_PI;
    return a - M_PI;
}
ComplexNumber::ComplexNumber(const std::complex<double>& c){
    if(std::fabs(c) < 1e-20){
        abs = zero_cutoff - 10.;
        arg = 0.;
        return;
    }
	abs = std::log(std::fabs(c));
	arg = wrap(std::arg(c));
}

ComplexNumber::ComplexNumber(){
    abs = 0.;
    arg = 0.;
}
ComplexNumber::ComplexNumber(const double abs_, const double arg_){
	abs = abs_; arg = wrap(arg_);
}

bool ComplexNumber::isZero() const{
    return abs < zero_cutoff;
}

std::complex<double> ComplexNumber::ToComplex() const{
	return std::exp(std::complex<double>(abs, arg));
}

//Supposes that our number is real..
double ComplexNumber::ToDouble() const{
    return ToComplex().real();
}

ComplexNumber Pow(const ComplexNumber& c, int pow){
    return ComplexNumber(pow*c.abs, pow*c.arg);
}

ComplexNumber operator*(const ComplexNumber& c1, const ComplexNumber& c2){
	return ComplexNumber(c1.abs+c2.abs, c1.arg+c2.arg);
}

ComplexNumber operator/(const ComplexNumber& c1, const ComplexNumber& c2){
	return ComplexNumber(c1.abs-c2.abs, c1.arg-c2.arg);
}

ComplexNumber operator+(const ComplexNumber& c1, const ComplexNumber& c2){
	if(c1.abs>c2.abs){
		return c1*ComplexNumber(std::complex<double>(1., 0.) + (c2/c1).ToComplex());
	}
	return c2*ComplexNumber(std::complex<double>(1., 0.) + (c1/c2).ToComplex());
}

ComplexNumber operator-(const ComplexNumber& c1, const ComplexNumber& c2){
	if(c1.abs>c2.abs){
		return c1*ComplexNumber(std::complex<double>(1., 0.) - (c2/c1).ToComplex());
	}
	return c2*ComplexNumber(std::complex<double>(1., 0.) - (c1/c2).ToComplex());
}

bool operator==(const ComplexNumber& c1, const ComplexNumber& c2){
	bool abs_val = (std::fabs(c1.abs-c2.abs) < precision);
	bool arg_val = (std::fabs(c1.arg - c2.arg) < precision) || (std::fabs(std::fabs(c1.arg-c2.arg)- 2*M_PI) < precision);
	return abs_val && arg_val;
}

bool operator<(const std::complex<double>& c1, const std::complex<double>& c2){
	auto p_1 = std::pair<double, double>(c1.real(), c1.imag());
	auto p_2 = std::pair<double, double>(c2.real()-precision, c2.imag()-precision);
	return p_1 < p_2;
}

//First comparison of abs_value, and if equal lexicographical on real and imaginary
bool operator<(const ComplexNumber& c1, const ComplexNumber& c2){
	if(c1.abs < (c2.abs-precision)) return true;
	if(c1.abs > (c2.abs+precision)) return false;
	return (c1.ToComplex()) < (c2.ToComplex());
}

std::ostream& operator<<(std::ostream& os, const ComplexNumber& c){
    os << "exp(" << c.abs << "+"<< c.arg <<"i)";
    return os;
}
