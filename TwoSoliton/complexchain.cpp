#include "complexchain.hpp"

ComplexNode::ComplexNode(const ComplexNumber& c_){
	if(c_.arg < -1.*precision){
		multiplicity = -1;
		c = ComplexNumber(c_.abs, c_.arg+M_PI);	
	}
	else{
		c=c_;
		multiplicity = 1;
	}

}

ComplexNode::ComplexNode(const ComplexNumber& c_, int mult){
	if(c_.arg < -1.*precision){
		multiplicity = -1*mult;
		c = ComplexNumber(c_.abs, c_.arg+M_PI);	
	}
	else{
		c=c_;
		multiplicity = mult;
	}

}
ComplexNode::ComplexNode(){
    c = ComplexNumber();
    multiplicity = 0.;
}
bool operator<(const ComplexNode& c1, const ComplexNode& c2){
	return c1.c < c2.c;
}

ComplexChain::ComplexChain(const ComplexNode& c){
	terms = std::set<ComplexNode>{c};
}

ComplexChain::ComplexChain(){
	terms = std::set<ComplexNode>{};
}
ComplexChain::ComplexChain(const ComplexNumber& c){
    if(c.isZero()){
        terms = std::set<ComplexNode>{};
        return;
    }
    terms = std::set<ComplexNode>{ComplexNode(c)};
}
ComplexChain::ComplexChain(const std::complex<double>& c){
    auto com = ComplexNode(ComplexNumber(c));
    terms = std::set<ComplexNode>{com};
}

bool ComplexChain::isZero() const{
    return terms.empty();
}
double ComplexChain::ToDouble() const{
	auto comp = ComplexNumber();
	for(const auto& node : terms){
		auto c_multi = ComplexNumber(std::complex<double>(node.multiplicity, 0.));
		auto new_c = c_multi * node.c;
		comp = comp + new_c;
	}
    auto ret = comp.ToComplex();
    return ret.real();
}
ComplexChain ComplexChain::Simplify() const{
    auto comp = ComplexNumber();
    auto first = true;
    for(const auto& node : terms){
        if(first){
            auto multi = ComplexNumber(std::complex<double>(node.multiplicity, 0.));
            comp = multi * node.c;
            first = false; continue;
        }
		auto c_multi = ComplexNumber(std::complex<double>(node.multiplicity, 0.));
        auto new_c = c_multi * node.c;
        comp = comp + new_c;
    }
    return ComplexChain(comp);

}

ComplexNumber ComplexChain::OneTerm() const{
    auto node = *((Simplify().terms).begin());
    return ComplexNumber(node.multiplicity)*node.c;
}
void ComplexChain::operator+=(const ComplexNode& node){
	if(auto search = terms.find(node); search!=terms.end()){
		auto found_node = *(search);
		terms.erase(search);
		if((found_node.multiplicity + node.multiplicity)==0) return;
		auto new_node = ComplexNode(found_node.c, found_node.multiplicity + node.multiplicity);
		terms.insert(new_node);
		return;
	}

	terms.insert(node);
}
void ComplexChain::operator+=(const ComplexChain& chain){
    for(const auto& node : chain.terms){
        *this += node;
    }
}
void ComplexChain::operator-=(const ComplexNode& node){
	if(auto search = terms.find(node); search!=terms.end()){
		auto found_node = *(search);
		terms.erase(search);
		if((found_node.multiplicity - node.multiplicity==0)) return;
		auto new_node = ComplexNode(found_node.c, found_node.multiplicity - node.multiplicity);
		terms.insert(new_node);
		return;
	}

	terms.insert(ComplexNode(node.c, -1*node.multiplicity));
}

ComplexChain operator+(const ComplexChain& c1, const ComplexChain& c2){
	auto new_c = ComplexChain();
	for(const auto& node : c1.terms){
		new_c+=node;
	}
	for(const auto& node : c2.terms){
		new_c+=node;
	}
	return new_c;
}

ComplexChain operator-(const ComplexChain& c1, const ComplexChain& c2){
	auto new_c = ComplexChain();
	for(const auto& node : c1.terms){
		new_c+=node;
	}
	for(const auto& node : c2.terms){
		new_c-=node;
	}
	return new_c;
}

ComplexChain operator*(const ComplexChain& c1, const ComplexChain& c2){
	auto new_c = ComplexChain();
	for(const auto& node1 : c1.terms){
		for(const auto& node2 : c2.terms){
			auto n_node = ComplexNode(node1.c*node2.c, node1.multiplicity*node2.multiplicity);
			new_c+=n_node;
		}
	}
	return new_c;
}

std::ostream& operator<<(std::ostream& os, const ComplexChain& chain){
    if((chain.terms).empty()){
        os << "0"; return os;
    }
    os << "[";
    bool first = true;
    for(const auto& node : chain.terms){
        if(!first) os << " + ";
        if(node.multiplicity != 1) os << node.multiplicity;
        os << node.c;
        first = false;
    }
    os << "]";
    return os;
}
















