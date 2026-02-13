#include "cchain.hpp"

//Cnode constructors
cnode::cnode(){}
cnode::cnode(const cnumber& c){
    value = c;
    multiplicity = 1;
    next = NULL;
    if(c.inLHPlane()) {
        value = c.opposite();
        multiplicity = -1;
    }
}
cnode::cnode(const cnumber& c, int multiplicity_){
    value = c;
    multiplicity = multiplicity_;
    next = NULL;
    if(c.inLHPlane()){
        value = c.opposite();
        multiplicity*=-1;
    }
}

//Cchain constructors and destructors
cchain::cchain(){
    head = NULL;
    tail = NULL;
}

cchain::cchain(const cnumber& cnum){
    auto node = std::make_shared<cnode>(cnum);
    head = node;
    tail = node;
}

cchain::cchain(const cnumber& cnum, int multiplicity_){
    auto node = std::make_shared<cnode>(cnum, multiplicity_);
    head = node;
    tail = node;
}

bool cchain::isZero() const{
    return (head==NULL);
}

//Following operation adds a number to the chain, deleting nodes when necessary, and setting head/tail when necessary
//Keeps ordering of chain by operator< defined in cnumbers.hpp (lexicographical ordering of abs, then arg)
//With precision of equality defined in cnumbers.hpp
void cchain::Add(const cnumber& c, int multiplicity_){
    
    auto node = std::make_shared<cnode>(c, multiplicity_);
    //Two edge cases: NULL Head
    if(head == NULL){
        head = node;
        tail = node;
        return;
    }
    auto curr_node = head;
    auto next_node = curr_node->next;
    //And number is equal to head value, where there are two subcases
    if((head->value) == c){
        auto multi = multiplicity_ + (curr_node->multiplicity);
        if(multi==0){
            head = next_node;
            if(next_node==NULL) tail = NULL;
            return;
        }
        head->multiplicity = multi;
        return;

    }

    /*
     * Now main loop that has three subcases.
     * 1) c is between two nodes: we add a node
     * 2) c is equal to next_node-> value and multiplicities don't cancel out
     * 3) c is equal to next_node->value and multiplicities cancel: remove node
     */
    while((next_node!=NULL) && (c < (next_node->value))){
        curr_node = next_node;
        next_node = next_node->next;
    }
    if(next_node==NULL){
        curr_node->next = node;
        tail=node;
        return;
    }
    if((c == next_node->value)==false){ //dont feel like defining != operator >:(
        curr_node->next = node;
        node->next = next_node;
        return;
    }
    auto mul = multiplicity_ + next_node->multiplicity;
    if(mul!=0){
        next_node->multiplicity = mul;
        return;
    }
    curr_node->next=(next_node->next);
}

cchain Add(const cchain& c1, const cchain& c2){
    auto final_chain = cchain();
    auto node = c1.head;
    while(node!=NULL){
        final_chain.Add(node->value, node->multiplicity);
        node=node->next;
    }
    node = c2.head;
    while(node!=NULL){
        final_chain.Add(node->value, node->multiplicity);
        node = node->next;
    }
    return final_chain;
}
cchain Multiply(const cchain& c1, const cchain& c2){
    auto final_c = cchain();    
    auto node_1 = c1.head;
    while(node_1!=NULL){
        auto node_2 = c2.head;
        while(node_2!=NULL){
            auto mul = (node_1->multiplicity) * (node_2->multiplicity);
            auto cnum = (node_1->value)*(node_2->value);
            final_c.Add(cnum, mul);
            node_2=node_2->next;
        }
        node_1=node_1->next;
    }
    return final_c;
}

// c1 / c2
// I think i have to be careful because nodes can be modified, to be sure deep copy everytime?
//Add, and Multiply deep copy, because they create nodes everytime by use of the .Add(cnumber, multiplicity) method.

cchain Divide(const cchain& c1, const cchain& c2){
    auto reste = cchain();
    reste = Add(reste, c1);
    auto quotient = cchain();
    while(!(reste.isZero())){

        if((reste.tail->value).lessthan((c1.tail->value))) {
            std::cout << "reste is " << reste.tail->value << std::endl;
            std::cout << "C1 is " << c1.tail->value << std::endl;
            std::cout << "C2 is " << c2 << std::endl;
            throw std::runtime_error("Non perfect division");
        }
        
        auto v_1 = reste.head->value; auto m_1 = reste.head->multiplicity;
        auto v_2 = c2.head->value; auto m_2 = c2.head->multiplicity;
        if(m_1%m_2!=0) {
            std::cout << "C1 is " << c1 << std::endl;
            std::cout << "C2 is " << c2 << std::endl;
            throw std::runtime_error("Non perfect division");
        }
        auto mul = m_1/m_2;
        auto v = v_1/v_2;
        quotient.Add(v, mul);
        auto new_term = cchain(v, -mul);
        new_term = Multiply(new_term, c2);
        reste = Add(reste, new_term);
    }
    
    return quotient;
}


std::ostream& operator<<(std::ostream& os, const cchain& c){
    auto curr_node = c.head;
    if(curr_node==NULL) {os << "0"; return os;}
    os <<"(" <<curr_node->multiplicity<< curr_node->value;
    curr_node = curr_node->next;
    while(curr_node!=NULL){
        os <<" + ";
        os <<curr_node->multiplicity<<curr_node->value;
        curr_node = curr_node->next;
    }
    os << ")";
    return os;
}



















