#include "polynomial.hpp"
#include <algorithm>

//Monomial Class

Monomial::Monomial(int d, int x_d){
    degree = d; x_degree = x_d; coefficient = cchain();
}
Monomial::Monomial(const cnumber& c, int d, int x_d){
    degree = d; x_degree = x_d; coefficient = cchain(c);
}

Monomial::Monomial(const cchain& chain, int d, int x_d){
    degree = d; x_degree = x_d; coefficient = chain;
}

bool Monomial::isZero() const{
    return coefficient.isZero();
}
Monomial Monomial::operator+(const Monomial& m2) const{
    if(degree!=m2.degree || x_degree != m2.x_degree) 
        throw std::runtime_error("Non equal monomials being added");
    auto n_chain = cchain::Add(coefficient, m2.coefficient);
    return Monomial(n_chain, degree, x_degree);
}

Monomial Monomial::operator-(const Monomial& m2) const{
    if(degree!=m2.degree || x_degree != m2.x_degree) 
        throw std::runtime_error("Non equal monomials being added");
    auto one = cnumber(0., 1.);
    auto minus_one = cchain::cchain(one, -1);
    auto negative_m2 = cchain::Multiply(minus_one, m2.coefficient);
    auto n_chain = cchain::Add(coefficient, negative_m2);
    return Monomial(n_chain, degree, x_degree);
}

Monomial Monomial::operator*(const Monomial& m2) const{
    auto n_d = degree + m2.degree;
    auto n_d_x = x_degree + m2.x_degree;
    auto n_coeff = cchain::Multiply(coefficient, m2.coefficient);
    return Monomial(n_coeff, n_d, n_d_x);
}

//Next is also only perfect division
Monomial Monomial::operator/(const Monomial& m2) const{
    auto n_d = degree-m2.degree;
    auto n_d_x = x_degree - m2.x_degree;
    if(n_d < 0 || n_d_x < 0 || (n_d_x > n_d)) throw std::runtime_error("Non perfect monomial division");
    auto new_c = cchain::Divide(cofficient, m2.coefficient);
    return Monomial(new_c, n_d, n_d_x);
}


//Polynomial implementation
std::pair<int, int> IdxToPair(int idx){
    float discriminant = 1. + static_cast<float>(8*idx);
    int d = std::floor(-1. + std::sqrt(discriminant));
    int n = idx - d*(d+1)/2;
    return std::pair<int, int>(d, n);
}

int PairToIdx(std::pair<int, int> degree){
    const auto& [d, n] = degree;
    return d*(d+1)/2 + n;
}

Polynomial::Polynomial(int d){
    auto coeffs = std::vector<Monomial>{}
    for(int idx = 0; idx < (d+1)*(d+2)/2; idx++){
        auto [d_, n] =IdxToPair(idx);
        auto new_M = Monomial(d_, n);
        coeffs.emplace_back.(new_M);
    }
    coefficients = coeffs;
    max_degree = d;
    leading_idx = -1;
}
Polynomial::Polynomial(const cchain& c, int d, int d_x){
    auto coeffs = std::vector<Monomial>{}
    for(int idx = 0; idx < (d+1)*(d+2)/2; idx++){
        auto [d_, n] =IdxToPair(idx);
        auto new_M = Monomial(d, n);
        if(n == d_x) new_M = Monomial(c, d, d_x);
        coeffs.emplace_back.(new_M);
    }
    coefficients = coeffs;
    max_degree = d;
    leading_idx = PairToIdx(std::pair<int, int>(d, d_x));
}

bool Polynomial::isZero() const{
    return (leading_idx == -1);
}
void Polynomial::SetLeading(){
    leading_idx = -1;
    for(int i=0; i<(degree+1)*(degree+2)/2; i++){
        if(coefficients[i].isZero()) continue;
        leading_idx = i;
    }
}

//Next supposes leading_idx is well set;
Polynomial Polynomial::Reduce(){
    auto [d, n] = IdxToPair(leading_idx);
    auto new_P = Polynomial(d);
    for(int i = 0; i<(d+1)*(d+2); i++){
        new_P.coefficients[i] = coefficients[i];
    }
}

Polynomial Polynomial::operator+(const Polynomial& P2) const{
    auto d = std::max(P2.degree, degree);
    auto new_P = Polynomial(d);
    for(int i=0; i<(degree+1)*(degree+2)/2; i++){
        new_P.coefficients[i] = new_P.coefficients[i] + coefficients[i];
    }
    for(int i=0; i<(P2.degree+1)*(P2.degree+2)/2; i++){
        new_P.coefficients[i] = new_P.coefficients[i] + P2.coefficients[i];
    }
    return new_P;
    new_P.SetLeading();
    new_P.Reduce();
}

Polynomial Polynomial::operator-(const Polynomial& P2) const{
    auto d = std::max(P2.degree, degree);
    auto new_P = Polynomial(d);
    for(int i=0; i<(degree+1)*(degree+2)/2; i++){
        new_P.coefficients[i] = new_P.coefficients[i] - coefficients[i];
    }
    for(int i=0; i<(P2.degree+1)*(P2.degree+2)/2; i++){
        new_P.coefficients[i] = new_P.coefficients[i] - P2.coefficients[i];
    }
    new_P.SetLeading();
    new_P.Reduce();
    return new_P;

}

Polynomial Polynomial::operator*(const Monomial& M) const{
    auto new_P = Polynomial(degree);
    for(int i = 0; i<(degree+1)*(degree+2)/2; i++){
        new_P.coefficients[i] = M*coefficients[i];
    }
    new_P.SetLeading();
    new_P.Reduce();
    return new_P;

}

Polynomial Polynomial::operator*(const Polynomial& P2) const {
    auto new_P = Polynomial(degree+P2.degree);
    for(int i = 0; i<i(degree+1)*(degree+2)/2; i++){
        for(int j=0; j<(P2.degree+1)*(P2.degree+2)/2; j++){
            auto [d_1, n_1] = IdxToPair(i);
            auto [d_2, n_2] = IdxToPair(j);
            auto l = PairToIdx(std::pair<int, int>(d_1+d_2, n_1+n_2));
            new_P.coefficients[l] = 
                coefficients[i]*P2.coefficients[j] + new_P.coefficients[l];
        }
    }
    new_P.SetLeading();
    new_P.Reduce();
    return new_P;

}


//is a perfect division
friend Polynomial Divide(const Polynomial& P1, const Polynomial& P2){
    auto reste = P1;
    if(P2.degree > P1.degree) throw std::runtime_error("Non perfect polynomial division");
    if(P2.isZero()) throw std::runtime_error("Division by zero polynomial");
    auto quotient = Polynomial(P1.degree - P2.degree);
    auto leading_2 = P2.coefficients[P2.leading_idx];
    while(!reste.isZero()){
       auto l_r = reste.coefficients[reste.leading_idx];
       auto next_term = l_r/leading_2;
       auto poly = Polynomial(next_term.coefficient,next_term.degree, next_term.x_degree);
       quotient = quotient + poly;
       reste = reste - (P2*next_term);
    }
    return quotient;
}






