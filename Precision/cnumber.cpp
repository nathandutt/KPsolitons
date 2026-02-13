#include <complex>
#include <cmath>
#include "cnumber.hpp"
const comp I = comp(0., 1.);
const comp one = comp(1., 0.);

/*
 * Helper Functions for complex number wrapping
 * We need make sure arg is in [-pi, pi]
 */
long double wrap_pi(long double a) {
    a = std::fmod(a + PI, 2*PI);
    if (a < 0) a += 2*PI;
    a = a - PI; //Now to handle exp(0) + exp(i*PI) not being evaluated to zero we need to add
    if(std::fabs(a-PI) < equality_precision) return a-2*PI;
    return a;
}

int max(int a, int b){
    if(a>b) return a;
    return b;
}

/*
 *Constructors for cnumber
 *Either initialize fields directly or with a complex number.
 */
cnumber::cnumber(){
    arg = 0.;
    abs = 0.;
}
cnumber::cnumber(const comp c){
    arg = std::arg(c); //guaranteeed here
    abs = std::log(std::abs(c));
    
}
cnumber::cnumber(const real arg_, const real abs_){

    arg = wrap_pi(arg_); abs = abs_; 
    //Keep in mind here abs_ is actually log of absolute value of number

}

/*
 *Comparison operations on cnumber
 *There are multiple. Because we need different types of equality
 */


//There can be some "dangerous" numerical behaviour here: because comparison is lexicographical
//Sometimes wont simplify exp(0) + exp(i*pi) to 0 because they aren't close enough in lexicographical order.
//We try to mitigate that by only wrapping to pi if we're far enough away from zero negatively.
bool cnumber::inLHPlane() const{
    return arg < -1.*equality_precision;
}

cnumber cnumber::opposite() const{
    return cnumber(arg + PI, abs);
}

//This is used later for division purposes (knowing when division is not perfect) operator< doesn't work because it compares complex numbers too
bool cnumber::lessthan(const cnumber& c2) const{
    if(abs < (c2.abs-equality_precision)) return true;
    return false;
}
bool cnumber::operator==(const cnumber& c2) const{
   bool im_equal = (std::fabs(c2.arg - arg) < equality_precision) || (std::fabs(std::fabs(c2.arg - arg) < equality_precision));
   bool r_equal = (std::fabs(c2.abs - abs) < equality_precision) ;
   return im_equal && r_equal;
}

//Supposes we've clipped arg into [0, pi] or maybe [-pi, pi]
bool cnumber::operator<(const cnumber& c2) const{
    if(abs < (c2.abs-equality_precision)) return true;
    else if(fabs(abs - c2.abs) < equality_precision) return (arg < (c2.arg-equality_precision));
    return false;
}
cnumber cnumber::operator-(const cnumber& c2) const{
    auto c_dom = cnumber(arg, abs); auto c_weak = c2;
    if(abs < c2.abs){c_weak = c_dom; c_dom = c2;}
    auto delta_c = cnumber(c_weak.arg - c_dom.arg, c_weak.abs - c_dom.abs);
    auto num = one - std::exp(std::complex<long double>(delta_c.abs, delta_c.arg));
    auto c_small = cnumber(num);
    return c_dom*c_small;
}    

cnumber cnumber::operator+(const cnumber& c2) const{
    auto c_dom = cnumber(arg, abs); auto c_weak = c2;
    if(abs < c2.abs){c_weak = c_dom; c_dom = c2;}
    auto delta_c = cnumber(c_weak.arg - c_dom.arg, c_weak.abs - c_dom.abs);
    auto num = one + std::exp(std::complex<long double>(delta_c.abs, delta_c.arg));
    auto c_small = cnumber(num);
    return c_dom*c_small;
}    

cnumber cnumber::operator/(const cnumber& c2) const{
   return cnumber(arg-c2.arg, abs-c2.abs); 
}    

cnumber cnumber::operator*(const cnumber& c2) const{
   return cnumber(arg+c2.arg, abs+c2.abs); 
}    

cnumber cnumber::operator*(const double d) const{
    auto d_compl = std::complex<long double>(d, 0.);
    auto d_cnum = cnumber(d_compl);
    return cnumber(arg+d_cnum.arg, abs + d_cnum.abs);
}


std::ostream& operator<<(std::ostream& os, const cnumber& c){
    os << "exp(" << c.abs <<"+"<<c.arg<<"i)";
    return os;
}
