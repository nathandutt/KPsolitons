#include "basicpolynomial.hpp"

BasicMonomial::BasicMonomial(){
    x=0; y=0; t=0; value =std::complex<double>(0., 0.);
}

BasicMonomial::BasicMonomial(const Monomial& m){
    x=m.x; y = m.y; t=m.t;
    value = ((m.value).OneTerm()).ToComplex();
}

bool operator<(const BasicMonomial& m1, const BasicMonomial& m2){
   if(m1.x < m2.x) return true;
   if(m1.x > m2.x) return false;
   if(m1.y < m2.y) return true;
   if(m1.y > m2.y) return false;
   return m1.t < m2.t;
}

BasicPolynomial::BasicPolynomial(){
    terms = std::set<BasicMonomial>{};
}

BasicPolynomial::BasicPolynomial(const Polynomial& P){
    terms = std::set<BasicMonomial>{};
    for(const auto& mon : P.terms){
        terms.insert(BasicMonomial(mon));
    }
}

double BasicPolynomial::Evaluate(const double x, const double y, const double t) const{
    std::complex<double> res = 0.;
    for(const auto& m : terms){
        auto n_t = std::pow(x, m.x)*std::pow(y, m.y)*std::pow(t, m.t)*m.value;
        res+= n_t;
    }
    return res.real();
}

std::ostream& operator<<(std::ostream& os, const BasicMonomial& M){
    os << M.value;
    if(M.x > 0){
        os << "X^" << M.x;
    }
    if(M.y > 0){
        os << "Y^" << M.y;
    }
    if(M.t > 0){
        os << "t^" << M.t;
    }
    return os;
}
std::ostream& operator<<(std::ostream& os, const BasicPolynomial& P){
    for(const auto& m : P.terms){
        os << m << " + ";
    }
    return os;
}
