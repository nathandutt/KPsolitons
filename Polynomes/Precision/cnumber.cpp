#include <complex>
#include <cmath>
#include "cnumber.hpp"
#include <algorithm>
typdef std::complex<real> comp;
const comp I = comp(0., 1.);
const comp one = comp(1., 0.);

//We need make sure arg is in [-pi, pi]
long double wrap_pi(long double a) {
    a = std::fmod(a + PI, 2*PI);
    if (a < 0) a += 2*PI;
    return a - PI;
}
cnumber::cnumber(const comp& c){
    arg = std::arg(c); //guaranteeed here
    abs = std::log(std::abs(c));
    
}
cnumber::cnumber(const real arg_, const real abs_){

    arg = wrap_pi(arg_); abs = abs_;

}

bool cnumber::inLHPlane() const{
    return arg < 0.;
}
cnumber cnumber::opposite() const{
    return cnumber(arg + PI, abs);
}

cnumber::operator==(const cnumber& c2) const{
   bool val = (std::fabs(c2.arg - arg) < equality_precision) && (std::fabs(c2.abs - abs) < equality_precision) ;
   return val;
}

//Supposes we've clipped arg into [0, pi] or maybe [-pi, pi]
cnumber::operator<(const cnumber& c2) const{
    if(abs < (c2.abs-equality_precision)) return true;
    else if(fabs(abs - c2.abs) < equality_precision) return (arg < (c2.arg-equality_precision));
    return false;
}
cnumber::operator-(const cnumber& c2) const{
    auto c_dom = cnumber(arg, abs); auto c_weak = c2;
    if(abs < c2.abs){c_weak = c_dom; c_dom = c2;}
    auto delta_c = cnumber(c_weak.arg - c_dom.arg, c_weak.abs - c_dom.abs);
    auto num = one - std::exp(std::complex<long double>(delta_c.abs, delta_c.arg));
    auto c_small = cnumber(num);
    return c_dom*c_small;
}    

cnumber::operator+(const cnumber& c2) const{
    auto c_dom = cnumber(arg, abs); auto c_weak = c2;
    if(abs < c2.abs){c_weak = c_dom; c_dom = c2;}
    auto delta_c = cnumber(c_weak.arg - c_dom.arg, c_weak.abs - c_dom.abs);
    auto num = one + std::exp(std::complex<long double>(delta_c.abs, delta_c.arg));
    auto c_small = cnumber(num);
    return c_dom*c_small;
}    

cnumber::operator/(const cnumber& c2) const{
   return cnumber(arg-c2.arg, abs-c2.abs); 
}    

cnumber::operator*(const cnumber& c2) const{
   return cnumber(arg+c2.arg, abs+c2.abs); 
}    


