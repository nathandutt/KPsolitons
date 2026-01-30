#include <stdio.h>
#include <iostream>
#include <vector>
#include <complex>
using namespace std;

const double zero = 1e-10;

typedef pair<int,int> Monomial;

auto max(int a, int b){
    if(a>b) return a;
    return b;
}

struct Polynomial{
    
    int degree;
    vector<vector<complex>> t;
    int leading;
    
    Polynomial(int d){
        degree = d;
        vector<vector<complex>> t_(d*(d+3)/2);
        t=t_;
        leading = t.size();
    }

    //Creates a monomial, i know there are more efficient ways of doing all that i am DOINGNGGNNOIDO starting by not using a multiplication by a monomial to just shift!!!!!!!!!! BUt i am lazy
    Polynomial(int i, int j){
        degree = i+j;
        idx = j*(j+1)/2;
        vector<vector<complex>> t_(d*(d+3)/2);
        t_[idx] = complex<double>(1., 0.);
        t = t_;
    }

    bool isZero() const{
        return (leading >= t.size())
    }
    void SetLeading(){
        int i;
        for(i=0; i<t.size(); i++){
            if(fabs(t[i]) > zero) break;
        }
        leading = i;
    }
    

    //Supposes lexicographical order with x power priority
    complex Get(int i, int j) const{
        if((i+j) > degree) return complex<double>(0., 0.);
        auto idx = (degree - i)*(degree-i+1) /2 + degree - i -j;
        return t[idx];
    }

    void Add(int i, int j, complex val){
        auto idx = (degree - i)*(degree-i+1) /2 + degree - i -j;
        t[idx] += val;
    }
     //Sorry for all the brackets
    Polynomial operator*(const Polynomial& P2) const {
        auto new_d = P2.degree + degree;
        auto new_P = Polynomial(new_d);

        for(int i_1 = 0; i_1 <= degree; i_1++){
            for(int j_1=0; j_1 <= degree; j_1++){
                for(int i_2=0;i_2 <=P2.degree; i_2++){
                    for(int j_2; j_2 <=P2.degree; j_2++){
                        auto new_addition = Get(i_1,j_1)*P2.Get(i_2,j_2);
                        new_P.Add(i_1+i_2, j_1+j_2, new_addition);
                    }
                } 
            }
        }
        new_P.Setx_Leading()
        return new_P;
    }

    Polynomial operator+(const Polynomial& P2) const{
        auto new_d = max(degree, P2.degree);
        auto new_P = Polynomial(new_d);

        for(int i=0; i<=new_d; i++){
            for(int j=0; j<=new_d; j++){
                auto val = Get(i, j) + P2.Get(i,j);
                new_P.Add(i, j, val);
            }
        }
        new_P.SetLeading();
        return new_P;
    }

    
    Polynomial operator-(const Polynomial& P2) const{
        auto new_d = max(degree, P2.degree);
        auto new_P = Polynomial(new_d);

        for(int i=0; i<=new_d; i++){
            for(int j=0; j<=new_d; j++){
                auto val = Get(i, j) - P2.Get(i,j);
                new_P.Add(i, j, val);
            }
        }
        new_P.SetLeading();
        return new_P;
    }

    //Hardest to implement is perfect division
    // And i know all my implementation is cringe in terms of operations, but for now anyways our polynomials won't be too large....
    Polynomial operator/(const Polynomial& P2) const{
        auto P1 = P;
        
    }
};

int main(){

}
