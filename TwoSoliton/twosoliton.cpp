
#include <string>
#include <vector>
#include <fstream>
#include "pointevolution.hpp"
using namespace std;
typedef vector<vector<Polynomial>> Matrix;

ostream& operator<<(ostream& os, const Matrix& M){
    int N = M.size();
    os << "--------------------" << endl;
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            os << M[i][j] << " -|| ";
        }
        os << endl;
    }
    os << "--------------------" << endl;
    return os;
}



void ReadParameters(fstream& file, vector<ComplexNumber>& k_s, vector<ComplexNumber>& offsets){
    double re_k, im_k, re_o, im_o;
    while(file >> re_k >> im_k >> re_o >> im_o){
        auto n_k = complex<double>(re_k, im_k);
        auto n_o = complex<double>(re_o, im_o);
        k_s.emplace_back(n_k);
        k_s.emplace_back(conj(n_k));
        offsets.emplace_back(n_o);
        offsets.emplace_back(conj(n_o));
    }
}

Polynomial DiagonalTerm(const ComplexNumber& k, const ComplexNumber& offset){
    auto one = ComplexNumber(0., 0.);
    auto X_term = Polynomial(Monomial(1, 0, 0, ComplexChain(one))); 
    auto Y_term = Polynomial(Monomial(0, 1, 0, ComplexChain(ComplexNumber(-2.)*k)));
    auto T_term = Polynomial(Monomial(0, 0, 1, ComplexChain(ComplexNumber(12.)*k*k))); 
    auto Constant  = Polynomial(offset);
    return X_term + Y_term + T_term + Constant;
}

Polynomial OffDiagonal(const ComplexNumber& k_i, const ComplexNumber& k_j){
    auto I = ComplexNumber(complex<double>(0., 1.));
    return Polynomial(Monomial(0, 0, 0, ComplexChain(I/(k_j-k_i))));
}
Matrix Hirota(const vector<ComplexNumber>& k_s, const vector<ComplexNumber>& offsets){
    auto M = Matrix();
    for(unsigned int i =0; i<k_s.size(); i++){
        auto line = std::vector<Polynomial>{};
        for(unsigned int j=0; j<k_s.size(); j++){
            if(i==j) 
                line.emplace_back(DiagonalTerm(k_s[i], offsets[i]));
            else 
                line.emplace_back(OffDiagonal(k_s[i], k_s[j]));
        } 
        M.emplace_back(line);
    }
    return M;
}
ostream& operator<<(ostream& os, const std::set<int>& s){
    for(const auto i : s){
        os << i << " ";
    }
    return os;
}
auto Submatrix(const Matrix& M, const set<int>& i_idx, const set<int>& j_idx){
    auto ret = Matrix();
    for(const auto i : i_idx){
        auto line = std::vector<Polynomial>{};
        for(const auto j: j_idx){
            line.emplace_back(M[i][j]);
        }
        ret.emplace_back(line);
    }
    return ret;
}
auto Minor(const Matrix& M, const set<int>& i_idx, const set<int>& j_idx){
    if(i_idx.size() == 2){
        auto i_it = i_idx.begin(); auto j_it = j_idx.begin();
        auto i1 = *(i_it); auto j1 = *(j_it);
        ++i_it; ++j_it;
        auto i2 = *(i_it); auto j2 = *(j_it);
        return (M[i1][j1]*M[i2][j2]) - (M[i1][j2]*M[i2][j1]);
        
    }
    auto ret = Polynomial();
    auto idx = 1;
    for(const auto i: i_idx){
        auto j = *(j_idx.begin());
        //cout << "Removing 1st column and "<< idx <<"th row to get" << endl;
        auto sign = ComplexNumber(complex<double>(((idx+1)%2==0)?1. : -1.,0.));
        auto polysign = Polynomial(Monomial(0, 0, 0, sign));
        //cout << "Sign is " << polysign << endl;
        auto n_i_idx = i_idx;
        auto n_j_idx = j_idx;
        n_i_idx.erase(n_i_idx.find(i));
        n_j_idx.erase(n_j_idx.find(j));
        //cout << Submatrix(M, n_i_idx, n_j_idx) << endl;
        auto m = polysign*M[i][j]*Minor(M, n_i_idx, n_j_idx);
        //cout << "i, j, m" << i << " " << j << " " << m << endl;
        ret = ret + m;
        idx++;
       
    }
    return ret;
}
Polynomial Determinant(const Matrix& M){
    auto i_idx = set<int>{};
    auto j_idx = set<int>{};
    for(int i = 0; i<(int)(M.size()); i++){
        i_idx.insert(i);
        j_idx.insert(i);
    }
    return Minor(M, i_idx, j_idx);
}

int main(){
    const string filename = "params.txt";
    fstream file;
    file.open(filename, ios::in);
    double timestep = 0.001;
    auto k_s = std::vector<ComplexNumber>{};
    auto offsets = std::vector<ComplexNumber>{};

    ReadParameters(file, k_s, offsets);
    file.close();
    auto M = Hirota(k_s, offsets);
    auto P = Determinant(M);
   // cout << P << endl;
    //cout << P.Simplify() << endl;
    //

    auto A = Derive(Derive(P, 0), 0) - Derive(P, 0)*Derive(P, 0);
    auto B = Derive(Derive(P, 0), 1) - Derive(P, 1)*Derive(P, 0);
    auto Det = Derive(A, 0)*Derive(B, 1) - Derive(B, 0)*Derive(A, 1);
    auto T_x =  Derive(A, 1)*Derive(B, 2) - Derive(B, 1)*Derive(A, 2);
    auto T_y =  Derive(A, 2)*Derive(B, 0) - Derive(B, 2)*Derive(A, 0);
    Det = Det.Simplify();
    T_x = T_x.Simplify();
    T_y = T_y.Simplify();
    auto B_Det = BasicPolynomial(Det);
    auto B_T_x = BasicPolynomial(T_x);
    auto B_T_y = BasicPolynomial(T_y);

//    cout << Det << endl << endl;
//    cout << T_x << endl << endl;
//    cout << T_y << endl << endl;
//
    cout << B_Det << endl;
    const string writeto = "evolution.csv";
    file.open(writeto, ios::out);
    auto pt = Point{-120.5, 0., -10.};
    while(pt.t < 10.){
       pt.Write(file);
       EvolvePoint(pt,timestep, B_Det, B_T_x, B_T_y);
    }
    

}
