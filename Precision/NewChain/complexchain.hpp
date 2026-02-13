#include "complexnumber.hpp"
#include <memory>
#include <set>

class ComplexChain{
    private:
        std::set<ComplexNumber, operator<> numbers;
    public:
        ComplexChain();
        ComplexChain(const ComplexNumber& c);
        
        bool isZero() const;
        void operator+(const ComplexNumber& c);

        friend ComplexChain operator*(const ComplexChain& c1, const ComplexChain& c2);
        friend ComplexChain operator/(const ComplexChain& c1, const ComplexChain& c2);
        friend ComplexChain operator+(const ComplexChain& c1, const ComplexChain& c2);
        friend ComplexChain operator-(const ComplexChain& c1, const ComplexChain& c2);

        friend std::ostream& operator<<(std::ostream& os, const ComplexChain& c);
}
