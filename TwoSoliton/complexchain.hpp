#include "complexnumber.hpp"
#include <set>
struct ComplexNode{
	ComplexNumber c;
	int multiplicity;
    ComplexNode();
	ComplexNode(const ComplexNumber& c);
	ComplexNode(const ComplexNumber& c, int mult);
};
bool operator<(const ComplexNode& c1, const ComplexNode& c2);

class ComplexChain{	
	std::set<ComplexNode> terms;

	public:
	ComplexChain();
	ComplexChain(const ComplexNode& c);
    ComplexChain(const ComplexNumber& c);
    ComplexChain(const std::complex<double>& c);
	void operator+=(const ComplexNode& node);
	void operator+=(const ComplexChain& chain);
	void operator-=(const ComplexNode& node);
	double ToDouble() const;
    ComplexChain Simplify() const;
    ComplexNumber OneTerm() const;
    bool isZero() const;
	friend ComplexChain operator+(const ComplexChain& c1, const ComplexChain& c2);
	friend ComplexChain operator-(const ComplexChain& c1, const ComplexChain& c2);
	friend ComplexChain operator*(const ComplexChain& c1, const ComplexChain& c2);
    friend std::ostream& operator<<(std::ostream& os, const ComplexChain& chain);
};

