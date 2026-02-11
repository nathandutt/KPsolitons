#include "polynomial.hpp"

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
    auto n_chain = Add(coefficient, m2.coefficient);
    return Monomial(n_chain, degree, x_degree);
}

Monomial Monomial::operator-(const Monomial& m2) const{
    if((degree!=m2.degree) || (x_degree != m2.x_degree)) 
        throw std::runtime_error("Non equal monomials being added");
    auto minus_one = cchain(cnumber(PI, 0.));
    auto negative_m2 = Multiply(minus_one, m2.coefficient);
    auto n_chain = Add(coefficient, negative_m2);
    return Monomial(n_chain, degree, x_degree);
}

Monomial Monomial::operator*(const Monomial& m2) const{
    auto n_d = degree + m2.degree;
    auto n_d_x = x_degree + m2.x_degree;
    auto n_coeff = Multiply(coefficient, m2.coefficient);
    return Monomial(n_coeff, n_d, n_d_x);
}

//Next is also only perfect division
Monomial Monomial::operator/(const Monomial& m2) const{
    auto n_d = degree-m2.degree;
    auto n_d_x = x_degree - m2.x_degree;
    if(n_d < 0 || n_d_x < 0 || (n_d_x > n_d)) {
        std::cout << "Monomial 2: " << m2 << std::endl;
        std::cout << "Monomial 1: " << *this << std::endl;
        throw std::runtime_error("Non perfect monomial division");}
    auto new_c = Divide(coefficient, m2.coefficient);
    return Monomial(new_c, n_d, n_d_x);
}

std::ostream& operator<<(std::ostream& os, const Monomial& M){
    os << M.coefficient;
    if(M.x_degree!=0) os << "X^" << M.x_degree;
    if((M.degree-M.x_degree) > 0) os << "Y^" <<(M.degree-M.x_degree);
    return os;
}

//Polynomial implementation
std::pair<int, int> IdxToPair(int idx){
    float discriminant = 1. + static_cast<float>(8*idx);
    int d = std::floor(-0.5 + 0.5*std::sqrt(discriminant));
    int n = idx - d*(d+1)/2;
    return std::pair<int, int>(d, n);
}

int PairToIdx(std::pair<int, int> degree){
    const auto& [d, n] = degree;
    return d*(d+1)/2 + n;
}

Polynomial::Polynomial(int d){
    auto coeffs = std::vector<Monomial>{};
    for(int idx = 0; idx < (d+1)*(d+2)/2; idx++){
        auto [d_, n] =IdxToPair(idx);
        auto new_M = Monomial(d_, n);
        coeffs.emplace_back(new_M);
    }
    coefficients = coeffs;
    degree = d;
    leading_idx = -1;
}
Polynomial::Polynomial(const cchain& c, int d, int d_x){
    auto coeffs = std::vector<Monomial>{};
    for(int idx = 0; idx < (d+1)*(d+2)/2; idx++){
        auto [d_, n] =IdxToPair(idx);
        auto new_M = Monomial(d_, n);
        if((n == d_x) && (d_ == d)) new_M = Monomial(c, d, d_x);
        coeffs.emplace_back(new_M);
    }
    coefficients = coeffs;
    degree = d;
    leading_idx = PairToIdx(std::pair<int, int>(d, d_x));
}

bool Polynomial::isZero() const{
    return (leading_idx == -1);
}
void Polynomial::SetLeading(){
    int idx = (degree+1)*(degree+2)/2 -1;
    while(idx>-1 && coefficients[idx].isZero()){
        idx--;
    }
    leading_idx = idx;
}

//Next supposes leading_idx is well set;
Polynomial Polynomial::Reduce(){
    auto [d, n] = IdxToPair(leading_idx);
    auto new_P = Polynomial(d);
    for(int i = 0; i<=leading_idx; i++){
        new_P.coefficients[i] = coefficients[i];
    }
    new_P.SetLeading();
    return new_P;
}

Polynomial Polynomial::operator+(const Polynomial& P2) const{
    auto d = max(P2.degree, degree);
    
    auto new_P = Polynomial(d);
    for(int i=0; i<=leading_idx; i++){
        new_P.coefficients[i] = new_P.coefficients[i] + coefficients[i];
    }
    for(int i=0; i<=P2.leading_idx; i++){
        new_P.coefficients[i] = new_P.coefficients[i] + P2.coefficients[i];
    }
    new_P.SetLeading();
    new_P.Reduce();
    return new_P;
}

Polynomial Polynomial::operator-(const Polynomial& P2) const{
    auto d = std::max(P2.degree, degree);
    auto new_P = Polynomial(d);
    for(int i=0; i<=leading_idx; i++){
        new_P.coefficients[i] = new_P.coefficients[i] + coefficients[i];
    }
    for(int i=0; i<=P2.leading_idx; i++){
        new_P.coefficients[i] = new_P.coefficients[i] - P2.coefficients[i];
    }
    
    new_P.SetLeading();
    new_P.Reduce();
    return new_P;

}

Polynomial Polynomial::operator*(const Monomial& M) const{
    auto d = degree + M.degree;
    auto new_P = Polynomial(d);
    for(int i = 0; i<=leading_idx; i++){
        auto [d_p, n_p] = IdxToPair(i);
        auto n_d = d_p + M.degree; auto n_d_x = n_p + M.x_degree;
        auto n_idx = PairToIdx(std::pair<int, int>(n_d, n_d_x));
        new_P.coefficients[n_idx] = coefficients[i]*M;
    }
    new_P.SetLeading();
    return new_P;

}

Polynomial Polynomial::operator*(const Polynomial& P2) const {
    auto new_P = Polynomial(degree+P2.degree);
    for(int i = 0; i<=leading_idx; i++){
        for(int j=0; j<=P2.leading_idx; j++){
            auto [d_1, n_1] = IdxToPair(i);
            auto [d_2, n_2] = IdxToPair(j);
            auto l = PairToIdx(std::pair<int, int>(d_1+d_2, n_1+n_2));
            new_P.coefficients[l] = 
                coefficients[i]*P2.coefficients[j] + new_P.coefficients[l];
        }
    }
    new_P.SetLeading();
    return new_P;

}

std::ostream& operator<<(std::ostream& os, const Polynomial& P){
    os << "-----" << std::endl << "Leading term = " << P.leading_idx << std::endl;
    if (P.leading_idx<0) {os << "0"; return os;}
    for(int i =0; i < P.leading_idx; i++){
         
        if(P.coefficients[i].isZero()) continue;
        os <<P.coefficients[i] << " + ";
    }
    os << P.coefficients[P.leading_idx];
    return os;
}

//is a perfect division
Polynomial Divide(const Polynomial& P1, const Polynomial& P2){
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
       reste = reste - P2*next_term;
    }
    return quotient;
}





