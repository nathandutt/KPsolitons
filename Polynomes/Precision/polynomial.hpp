#include "cchain.hpp"
#include <vector>

struct Monomial{
    int degree;
    int x_degree;
    cchain coefficient;
    
    Monomial(int d, int x_d);
    Monomial(const cnumber& c, int  d, int x_d);
    Monomial(const cchain& chain, int d, int x_d);
    bool isZero() const;
    bool operator<(const Monomial& m2) const; //Compare degree (lexicographical)
    Monomial operator+(const Monomial& m2) const; //Add only of same degree
    Monomial operator-(const Monomial& m2) const;
    Monomial operator*(const Monomial& m2) const;
    Monomial operator/(const Monomial& m2) const;

};

//Polynomial is vector of Monomials
//
std::pair<int, int> IdxToPair(int idx);
int PairToIdx(std::pair<int, int> degrees);

class Polynomial{
    private:
        std::vector<Monomial> coefficients;
        int degree;
        int leading_idx;

    public:
        Polynomial(int d);
        Polynomial(const cchain& c, int d, int d_x);
        bool isZero() const;
        void SetLeading();
        Polynomial Reduce();
        Polynomial operator+(const Polynomial& P2) const;
        Polynomial operator-(const Polynomial& P2) const;
        Polynomial operator*(const Polynomial& P2) const;
        Polynomial operator*(const Monomial& M) const;
        friend Polynomial Divide(const Polynomial& P1, const Polynomial& P2);
        
};
