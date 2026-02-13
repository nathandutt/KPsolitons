#include "polynomials.hpp"

//Constructors
Polynomial::Polynomial(){
    terms = std::set<Monomial>{};
}

Polynomial::Polynomial(const std::complex<long double>& c){
    terms = std::set<Monomial>{Monomial(c)};
}
Polynomial::Polynomial(Monomial M){
    terms = std::set<Monomial>{M};
}

//Operations

void Polynomial::operator+=(const Monomial& M){
    if(auto found = terms.find(M); found!=terms.end()){
        auto f_M = *found;
        auto n_M = f_M + M;
        terms.erase(found);
        if(!n_M.isZero()) terms.insert(n_M);
        return;
    }
    terms.insert(M);
}

void Polynomial::operator-=(const Monomial& M){
    if(auto found = terms.find(M); found!=terms.end()){
        auto f_M = *found;
        auto n_M = f_M - M;
        terms.erase(found);
        if(!n_M.isZero()) terms.insert(n_M);
        return;
    }
    terms.insert(M);
}

Polynomial operator+(const Polynomial& p1, const Polynomial& p2){
    auto n_p = Polynomial();
    for( const auto& M : p1.terms) n_p+=M;
    for( const auto& M : p2.terms) n_p+=M;
    return n_p;
}

Polynomial operator-(const Polynomial& p1, const Polynomial& p2){
    auto n_p = Polynomial();
    for( const auto& M : p1.terms) n_p+=M;
    for( const auto& M : p2.terms) n_p-=M;
    return n_p;
}

Polynomial operator*(const Polynomial& p1, const Polynomial& p2){
    auto n_p = Polynomial();
    for(const auto& M1 : p1.terms){
        for(const auto& M2: p2.terms){
            n_p += M1*M2;
        }
    }
    return n_p;
}

std::ostream& operator<<(std::ostream& os, const Polynomial& p){
    auto first = true;
    for(const auto& monomial : p.terms){
        if(!first) {
            os << " + ";
        }
        first = false;
        os << monomial;
    }
    return os;
}

Polynomial Derive(const Polynomial& p, int variable){
    auto n_p = Polynomial();
    for(const auto& M : p.terms){
        n_p+=Derive(M, variable);        
    }
    return n_p;
}
