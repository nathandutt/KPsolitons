#include "monomials.hpp"
#include <set>
class Polynomial{
    std::set<Monomial> terms;
    
    public:
        Polynomial();
        Polynomial(const std::complex<long double>& c);
        Polynomial(Monomial M);
        bool isZero() const;
         
        void operator+=(const Monomial& M);
        void operator -=(const Monomial& M);
        friend Polynomial operator+(const Polynomial& p1, const Polynomial& p2);
        friend Polynomial operator-(const Polynomial& p1, const Polynomial& p2);
        friend Polynomial operator*(const Polynomial& p1, const Polynomial& p2);
        friend Polynomial Derive(const Polynomial& p, int variable);
        friend std::ostream& operator<<(std::ostream& os, const Polynomial& p);
};
