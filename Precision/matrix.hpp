#include "polynomial.hpp"
#include <vector>

typedef std::vector<std::vector<Polynomial>> Matrix;
std::vector<cnumber> EvolveOffsets(const std::vector<cnumber>& initial_offsets, 
        const std::vector<cnumber>& k_s, const double time);
Polynomial OffDiagonalTerm(const cnumber& k_1, const cnumber& k_2);
Polynomial DiagonalTerm(const cnumber& k, const cnumber& offset);
Matrix Hirota(const std::vector<cnumber>& k_s, const std::vector<cnumber>& offsets);
void Bareiss(Matrix& M);
