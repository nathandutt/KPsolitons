#include "polynomials.hpp"
struct BasicMonomial{
    int x; int y; int t; std::complex<double> value;

    BasicMonomial();
    BasicMonomial(const Monomial& m);
};
bool operator<(const BasicMonomial& m1, const BasicMonomial& m2);
class BasicPolynomial{
    public:
    std::set<BasicMonomial> terms;
    BasicPolynomial();
    BasicPolynomial(const Polynomial& P);

    double Evaluate(const double x, const double y, const double t) const;
};

std::ostream& operator<<(std::ostream& os, const BasicMonomial& M);
std::ostream& operator<<(std::ostream& os, const BasicPolynomial& P);
