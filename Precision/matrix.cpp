#include "matrix.hpp"
#include <iostream>
typedef unsigned int uint;

//Take inital offsets and evolve them: += 12*t*k^2
std::vector<cnumber> EvolveOffsets(const std::vector<cnumber>& initial_offset,
        const std::vector<cnumber>& k_s, const double time){
    auto final_offset = std::vector<cnumber>{};
    
    for(uint i=0u; i<k_s.size(); i++){
        auto new_o = initial_offset[i] + (k_s[i]*k_s[i])*(12.*time);
        final_offset.emplace_back(new_o);
    }
    return final_offset;
}

//Off diagonal term of hirota matrix
Polynomial OffDiagonalTerm(const cnumber& k_i, const cnumber& k_j){
    auto diff = k_j - k_i;
    auto I = cnumber(PI/2., 0.);
    auto constant = cchain(I/diff);

    return Polynomial(constant, 0, 0);
}


//Diagonal terms of Hirota matrix
Polynomial DiagonalTerm(const cnumber& k, const cnumber& offset){
    auto one = cchain(cnumber(0., 0.));
    auto y_coeff = k * (-2.);
    auto constant = Polynomial(offset, 0, 0);
    auto X = Polynomial(one, 1, 1);
    auto Y = Polynomial(y_coeff, 1, 0);
    return (X + Y + constant);
}
Matrix Hirota(const std::vector<cnumber>& k_s, const std::vector<cnumber>& offsets){
   auto M = std::vector<std::vector<Polynomial>>{};
   for(uint i=0u; i<k_s.size(); i++){
       auto line = std::vector<Polynomial>{};
       for(uint j=0u; j<k_s.size(); j++){
           if(i==j) line.emplace_back(DiagonalTerm(k_s[i], offsets[i]));
           else line.emplace_back(OffDiagonalTerm(k_s[i], k_s[j]));
       }
       M.emplace_back(line);
   }
   return M;
}
void Bareiss(Matrix& M){
    std::cout << "Inside Bareiss" << std::endl;
    uint n = M.size();
    for(uint k = 0; k < n; k++){
        for(uint i =k+1; i< n; i++){
            for(uint j=k+1; j<n; j++){
                std::cout << i << " " << j << " " << k << std::endl;
                
                if(i==5 && j==5 && k==4){
                    std::cout << "M[i][j]: " << M[i][j] << std::endl
                              << "M[k][k]: " << M[k][k] << std::endl
                              << "M[i][k]: " << M[i][k] << std::endl
                              << "M[k][j]: " << M[k][j] << std::endl;
                }
                
                auto A = M[i][j]*M[k][k] - M[i][k]*M[k][j];
                auto B = Polynomial(cnumber(0., 0.), 0, 0);
                if(k!=0) B = M[k-1][k-1];
                if(B.isZero()) throw std::runtime_error("Dividing by zero polynomial");;
                M[i][j] = Divide(A, B);
            }
            M[i][k] = Polynomial(0);
        }
    }
    std::cout << "Finished Bareiss";

}
