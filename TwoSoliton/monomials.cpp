#include "monomials.hpp" 

Monomial::Monomial(int x_, int y_, int t_, std::complex<long double> val_){
    if(x_ < 0 || y_ < 0 || t_<0){
        x=0; y=0; t=0; value = std::complex<long double>(0., 0.);
    }
    else{
        x=x_; y=y_; t=t_; value = val_;
    }
}
Monomial::Monomial(std::complex<long double> val_){
    x=0; y=0; t=0; value = val_;
}
bool Monomial::isZero() const{
    return fabs(value) < equality_precision;
}

bool operator<(const Monomial& m1, const Monomial& m2){
    if(m1.x < m2.x) return true;
    if(m1.y < m2.y) return true;
    return m1.t < m2.t;
}

Monomial operator*(const Monomial& m1, const Monomial& m2){
    if(m1.isZero() || m2.isZero()){
        return Monomial(std::complex<long double>(0., 0.));
    }
    return Monomial(m1.x+m2.x, m1.y+m2.y, m1.t + m2.t, m1.value*m2.value);
}

Monomial operator+(const Monomial& m1, const Monomial& m2){
    if(m1.isZero()) return m2;
    if(m2.isZero()) return m1;
    if(m1.x != m2.x || m1.y != m2.y || m1.t != m2.t) throw std::runtime_error("Non equal degree monomials being added");
    return Monomial(m1.x, m1.y, m1.t, m1.value + m2.value);
}

Monomial operator-(const Monomial& m1, const Monomial& m2){
    if(m1.isZero()) return std::complex<long double>(-1., 0.)*m2;
    if(m2.isZero()) return m1;
    if(m1.x != m2.x || m1.y != m2.y || m1.t != m2.t) throw std::runtime_error("Non equal degree monomials being subtracted");
    return Monomial(m1.x, m1.y, m1.t, m1.value - m2.value);
}

bool operator==(const Monomial& m1, const Monomial& m2){
    return (m1.x==m2.x) && (m1.y==m2.y) && (m1.t == m2.t);
}

std::ostream& operator<<(std::ostream& os, const Monomial& m){
   if(m.isZero()) return os;
   os << m.value;
   if(m.x>0) os << "X^" << m.x;
   if(m.y>0) os << "Y^" << m.y;
   if(m.t>0) os << "T^" << m.t;
   return os;
}
Monomial Derive(const Monomial& m, int idx){
    switch (idx)
    {
        case 0:
            return Monomial(m.x-1, m.y, m.t, m.value*std::complex<long double>(m.x));     
        case 1:
            return Monomial(m.x, m.y-1, m.t, m.value*std::complex<long double>(m.y));     
        case 2:
            return Monomial(m.x, m.y, m.t-1, m.value*std::complex<long double>(m.t));     
        default:
            throw std::runtime_error("Unknown derivation idx");
    }
}
