#include "monomials.hpp"
class Polynomial{
    
    public:
        std::set<Monomial> terms;
        Polynomial();
        Polynomial(const ComplexChain& c);
        Polynomial(Monomial M);
        bool isZero() const;
         
        void operator+=(const Monomial& M);
        void operator -=(const Monomial& M);
        Polynomial Simplify() const;
        ComplexNumber Evaluate(const double x, const double y, const double t) const;
        friend Polynomial operator+(const Polynomial& p1, const Polynomial& p2);
        friend Polynomial operator-(const Polynomial& p1, const Polynomial& p2);
        friend Polynomial operator*(const Polynomial& p1, const Polynomial& p2);
        friend Polynomial Derive(const Polynomial& p, int variable);
        friend std::ostream& operator<<(std::ostream& os, const Polynomial& p);
};
