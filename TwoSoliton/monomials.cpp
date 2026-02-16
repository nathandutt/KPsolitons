#include "monomials.hpp" 

Monomial::Monomial(int x_, int y_, int t_, ComplexChain val_){
    if(x_ < 0 || y_ < 0 || t_<0){
        x=0; y=0; t=0; value = ComplexChain();
    }
    else{
        x=x_; y=y_; t=t_; value = val_;
    }
}
Monomial::Monomial(ComplexChain val_){
    x=0; y=0; t=0; value = val_;
}
bool Monomial::isZero() const{
    return value.isZero();
}

ComplexChain Monomial::Evaluate(const ComplexNumber& X,const ComplexNumber& Y, const ComplexNumber& T) const{
    return Pow(X, x)*Pow(Y, y)*Pow(T, t)*value;
}
Monomial Monomial::Simplify() const{
    return Monomial(x, y, t, value.Simplify());
}
bool operator<(const Monomial& m1, const Monomial& m2){
    if(m1.x < m2.x) return true;
    if(m1.x > m2.x) return false;
    if(m1.y < m2.y) return true;
    if(m1. y > m2.y) return false;
    return m1.t < m2.t;
}

Monomial operator*(const Monomial& m1, const Monomial& m2){
    if(m1.isZero() || m2.isZero()){
        return Monomial(ComplexChain());
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
    if(m1.isZero()) return ComplexChain(ComplexNumber(std::complex<double>(-1., 0.)))*m2;
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
            return Monomial(m.x-1, m.y, m.t, m.value*ComplexChain(m.x));     
        case 1:
            return Monomial(m.x, m.y-1, m.t, m.value*ComplexChain(m.y));     
        case 2:
            return Monomial(m.x, m.y, m.t-1, m.value*ComplexChain(m.t));     
        default:
            throw std::runtime_error("Unknown derivation idx");
    }
}
