#include "matrix.hpp"
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

void ReadSolitons(fstream& file, vector<cnumber>& k_s, vector<cnumber>& offsets){
    long double k_re, k_im, o_re, o_im;
    while(file >> k_re >> k_im >> o_re >> o_im){
        auto new_k = cnumber(complex<long double>(k_re, k_im));
        auto new_k_bar = cnumber(complex<long double>(k_re, -k_im));
        auto new_o = cnumber(complex<long double>(o_re, o_im));
        auto new_o_bar = cnumber(complex<long double>(o_re, -o_im));
        k_s.emplace_back(new_k); k_s.emplace_back(new_k_bar);
        offsets.emplace_back(new_o); offsets.emplace_back(new_o_bar); 
    }
}

int main(){
    cout << setprecision(3) << endl;
    //Initial definitions
    string filename = "params.txt";
    fstream file;
    auto k_s = vector<cnumber>{};
    auto offsets = vector<cnumber>{};
    double time;
    
    //Read from file
    file.open(filename, ios::in);
    file >> time;
    ReadSolitons(file, k_s, offsets);

    //Evolve Offsets
    offsets = EvolveOffsets(offsets, k_s, time);

    //Create Matrix
    auto M = Hirota(k_s, offsets);
    int N = M.size();

    //Do Bareiss and extract determinant
    Bareiss(M);
    auto P = M[N-1][N-1];

    //If all went well, Print. TODO: Concert cchain to float for final det.

    cout << P << endl;
}

