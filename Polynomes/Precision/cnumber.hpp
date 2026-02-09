#include<complex>
#include <numbers>

typedef long double real;
typedef std::complex<double> comp;

static const long double PI = std::numbers::pi_v<long double>;
static const long double equality_precision=1e-8;
class cnumber{
    private:
        long double abs;
        long double mod;
    public:
        cnumber(const comp& c);
        cnumber(const real& arg, const real& abs);
        bool inLHPlane() const;
        cnumber opposite() const;
        cnuber operator==(const cnumber& c2) const;
        cnumber operator<(const cnumber& c2) const;
        cnumber operator*(const cnumber& c2) const;
        cnumber operator+(const cnumber& c2) const;
        cnumber operator-(const cnumber& c2) const;
        cnumber operator/(const cnumber& c2) const;

};
