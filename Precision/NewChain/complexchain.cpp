#include "complexchain.hpp"

ComplexChain::ComplexChain(const ComplexNumber& c){
    numbers = std::set<ComplexNumber>(c);
}

bool ComplexChain::isZero() const{
    return numbers.empty();
}

void ComplexChain::operator+(const ComplexNumber& c){
    if(numbers.empty()){
        numbers.insert(c);
        return;
    }

}
