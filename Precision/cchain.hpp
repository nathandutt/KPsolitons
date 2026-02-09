#include "cnumber.hpp"

struct cnode{
    int multiplicity;
    cnumber value;
    cnode* next;

    cnode();
    cnode(const cnumber& cnum);
    cnode(const cnumber& cnum, int multiplicity_);
};
class cchain{
    private:
        cnode* head;
        cnode* tail;
    public:
        cchain();
        cchain(const cnumber& cnum);
        cchain(const cnumber& cnum, int multiplicity_);
        ~cchain();
        bool isZero() const;
        void Add(const cnumber& c, int multiplicity);
        long double Calculate();
        friend cchain Add(const cchain& c1, const cchain& c2);
        friend cchain Multiply(const cchain& c1, const cchain& c2);
        friend cchain Divide(const cchain& c1, const cchain& c2);
};
