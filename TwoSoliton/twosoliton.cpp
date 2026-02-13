
#include <string>
#include <vector>
#include <fstream>
#include "polynomials.hpp"

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



void ReadParameters(fstream& file, vector<complex<long double>>& k_s, vector<complex<long double>>& offsets){
    long double re_k, im_k, re_o, im_o;
    while(file >> re_k >> im_k >> re_o >> im_o){
        auto n_k = complex<long double>(re_k, im_k);
        auto n_o = complex<long double>(re_o, im_o);
        k_s.emplace_back(n_k);
        k_s.emplace_back(conj(n_k));
        offsets.emplace_back(n_o);
        offsets.emplace_back(conj(n_o));
    }
}

Polynomial DiagonalTerm(const complex<long double>& k, const complex<long double>& offset){
    auto X_term = Polynomial(Monomial(1, 0, 0, complex<long double>(1., 0.))); 
    auto Y_term = Polynomial(Monomial(0, 0, 0, (long double)(-2.)*k));
    auto T_term = Polynomial(Monomial(0, 0, 0, (long double)(12.)*k*k)); 
    auto Constant  = Polynomial(offset);
    return X_term + Y_term + T_term + Constant;
}

Polynomial OffDiagonal(const complex<long double>& k_i, const complex<long double>& k_j){
    auto I = complex<long double>(0., 1.);
    return Polynomial(I/(k_j-k_i));
}
Matrix Hirota(const vector<complex<long double>>& k_s, const vector<complex<long double>>& offsets){
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
        cout << "Removing 1st column and "<< idx <<"th row to get" << endl;
        auto sign = complex<long double>(((idx+1)%2==0) ? 1. : -1., 0.);
        cout << "Sign is " << sign << endl;
        auto n_i_idx = i_idx;
        auto n_j_idx = j_idx;
        n_i_idx.erase(n_i_idx.find(i));
        n_j_idx.erase(n_j_idx.find(j));
        cout << Submatrix(M, n_i_idx, n_j_idx) << endl;
        auto m = sign*M[i][j]*Minor(M, n_i_idx, n_j_idx);
        //cout << "i, j, m" << i << " " << j << " " << m << endl;
        ret = ret + m;
        idx++;
       
    }
    return ret;
}
Polynomial Determinant(const Matrix& M){
    cout << "Getting determinant of " << M << endl;
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

    auto k_s = std::vector<complex<long double>>{};
    auto offsets = std::vector<complex<long double>>{};

    ReadParameters(file, k_s, offsets);
    auto M = Hirota(k_s, offsets);
    auto P = Determinant(M);
    cout << M << endl;
    cout << P << endl;
}
