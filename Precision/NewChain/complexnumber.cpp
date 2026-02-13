#include "complexnumber.hpp"

//Helper function to put angle into [-PI, PI]
double wrap(const double angle){ 
    auto a = std::fmod(angle-PI,2*PI);
    if(a<0) a += 2*PI;
    return a-PI;
}

//Constructors
ComplexNumber::ComplexNumber(const std::complex<double> c){
    abs = std::log(std::abs(c));
    arg = std::arg(c);
}

ComplexNumber::ComplexNumber(const double abs_, const double arg_){
    abs = abs_; arg = wrap(arg_);
}

std::complex<double> ToComplex(const ComplexNumber& c){
    return std::exp(std::complex<double>(abs, arg));
}

//Basic Operations
ComplexNumber operator*(const ComplexNumber& c1, const ComplexNumber& c2){
    return ComplexNumber(c1.abs + c2.abs, c1.arg + c2.arg);
}

ComplexNumber operator/(const ComplexNumber& c1, const ComplexNumber& c2){
    return ComplexNumber(c1.abs - c2.abs, c1.arg - c2.arg);
}

ComplexNumber operator+(const ComplexNumber& c1, const ComplexNumber& c2){
    if(c1.abs < c2.abs) return c2+c1;
    auto c_div = c2/c1;
    auto sum = std::complex<double>(1., 0.) +  ToComplex(c_div);
    return c1*ComplexNumber(sum);
}
ComplexNumber operator-(const ComplexNumber& c1, const ComplexNumber& c2){
    if(c1.abs < c2.abs) return c2-c1;
    auto c_div = c2/c1;
    auto sum = std::complex<double>(1., 0.) -  ToComplex(c_div);
    return c1*ComplexNumber(sum);
}

//Comparisons

bool operator==(const ComplexNumber& c1, const ComplexNumber& c2){
    auto val = (std::fabs(c1.abs - c2.abs) < equality_precision) &&
               ((std::fabs(c1.arg - c2.arg) < equality_precision) 
                    ||(std::fabs(2*PI-std::fabs(c1.arg - c2.arg)) < equality_precision));
    return val;
}

//< comparator: first compare abs field, if equal: check that both aren't equal, if not compare complex values lexicographically on imaginary first then real..
bool operator<(const ComplexNumber& c1, const ComplexNumber& c2){
    if(c1.abs < c2.abs-equality_precision) return true;
    if(c1.abs > c2.abs+equality_precision) return false;
    if(std::fabs(c1.arg - c2.arg) < equality_precision || 
            (std::fabs(c1.arg - c2.arg) - 2*PI) < equality_precision) return false;

    auto c1_complex = ToComplex(c1);
    auto c2_complex = ToComplex(c2);
    auto c1_pair = std::pair<double, double>(c1_complex.imag, c1_complex.real);
    auto c2_pair = std::pair<double, double>(c2_complex.imag, c2_complex.real);

    return c1_pair < c2_pair;
}

bool operator>(const ComplexNumber& c1, const ComplexNumber& c2){
    return !(c1<c2) && !(c1==c2);
}

bool Close(const ComplexNumber& c1, const ComplexNumber& c2){
    //TODO
}






