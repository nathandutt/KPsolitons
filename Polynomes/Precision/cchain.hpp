#include "cnumber.hpp"

struct cnode{
    int multiplicity;
    cnumber value;
    cnode* next;

    cnode();
    cnode(cnumber cnum);
    cnode(cnumber cnum, int multiplicity_);
};
class cchain{
    private:
        cnumber* head;
        cnumber* tail;
    public:
        cchain();
        cchain(cnumber cnum);
        cchain(cnumber cnum, int multiplicity_);
        ~cchain();
        bool isZero() const;
        void Add(cnode* node);
        long double Calculate();
};
friend cchain Add(const cchain& c1, const cchain& c2);
friend cchain Multiply(const cchain& c1, const cchain& c2);
friend cchain Divide(const cchain& c1, const cchain& c2);
