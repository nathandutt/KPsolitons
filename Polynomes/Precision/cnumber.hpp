#include<complex>
#include <numbers>

typedef long double real;
typedef std::complex<real> comp;

static const long double PI = std::numbers::pi_v<long double>;
static const long double equality_precision=1e-8;
int max(int a, int b);
class cnumber{
    private:
        long double abs;
        long double arg;
    public:
        cnumber();
        cnumber(const comp c);
        cnumber(const real arg, const real abs);
        bool inLHPlane() const;
        cnumber opposite() const;
        bool operator==(const cnumber& c2) const;
        bool operator<(const cnumber& c2) const;
        cnumber operator*(const cnumber& c2) const;
        cnumber operator+(const cnumber& c2) const;
        cnumber operator-(const cnumber& c2) const;
        cnumber operator/(const cnumber& c2) const;

};
