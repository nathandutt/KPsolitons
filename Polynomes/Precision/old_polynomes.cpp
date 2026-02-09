#include <stdio.h>
#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <iostream>
using namespace std;

const double zero = 1e-8;

typedef pair<int,int> Monomial;

typedef long long bigint;
//Supposes a in Z, b > 0
bigint GCD(bigint a, bigint b){

    if(a==0 && b==0) return 1;
    if(a<0) return GCD(-a, b);
    if(a<b) return GCD(b, a);
    if(b==0) return a;
    return GCD(b, a%b);
}

auto MakeCoprime(bigint a, bigint b){
    if(a==0) return pair<bigint, bigint>(a, b);
    if(b==0) return pair<bigint, bigint>(a, b);
    auto c = GCD(a, b);
    return pair<bigint, bigint>(a/c, b/c);
}
struct Rational{
    bigint whole;
    bigint a;
    bigint b;
    Rational(){
        a=0; b=1;
    } 
    Rational(bigint a_, bigint b_){
        auto a_1 = a_; auto b_1 = b_;
        if(b_ < 0) {a_1 = -a_1; b_1 = -b_1;}
        auto [a_2, b_2] = MakeCoprime(a_1, b_1);
        a=a_2; b=b_2;
    }
    
    Rational(bigint a_){
        a =a_; b=1;
    }
    void Simplify(){
        auto [a_, b_] = MakeCoprime(a, b);
        a = a_; b = b_;
    }

    float Floatize() const{
        return static_cast<float>(a)/static_cast<float>(b);
    }
    
    Rational operator+(Rational r) const{
        auto [b_1, b_2]  = MakeCoprime(b, r.b);
        auto n_a = a*b_2 + r.a*b_1;
        auto n_b = b * b_2;
        auto res = Rational(n_a, n_b);
        res.Simplify();
        return res;
    }
    Rational operator-(Rational r) const{
        auto [b_1, b_2]  = MakeCoprime(b, r.b);
        auto n_a = a*b_2 - (r.a)*b_1;
        auto n_b = b * b_2;
        auto res = Rational(n_a, n_b);
        res.Simplify();
        return res;
    }
    Rational operator*(Rational r) const{
        if(a==0 || r.a==0) return Rational(0);
        auto [a_1, b_2] = MakeCoprime(a, r.b);
        auto [a_2, b_1] = MakeCoprime(r.a, b);
        return Rational(a_1*a_2, b_1*b_2);
    }

    Rational operator/(Rational r) const{
        if(r.a==0) {cout << "Division by zero " << endl; return Rational(0, 0);}
        auto [a_1, a_2] = MakeCoprime(a, r.a);
        auto [b_1, b_2] = MakeCoprime(b, r.b);
        return Rational(a_1*b_2, b_1*a_2);
    }

    void operator+=(Rational r){
        auto r1 = Rational(a, b);
        auto res = r1+r;
        a = res.a; b = res.b;
    }


    void operator-=(Rational r){
        auto r1 = Rational(a, b);
        auto res = r1-r;
        a = res.a; b = res.b;
    }

    void Print() const{
        cout << a << "/" << b;
    }



};

struct RComplex{
    Rational re;
    Rational im;
    
    RComplex(){
        re = Rational();
        im = Rational();
    }
    RComplex(Rational a, Rational b){
        re = a; im = b;
    }
    RComplex(bigint a, bigint b){
        re = Rational(a); im = Rational(b);
    }
    Rational sqabs()const {
        return re*re + im*im;
    }
    RComplex conj()const {
        return RComplex(re, Rational(-1)*im);
    }
    RComplex operator*(RComplex r) const {
        auto n_re =re * r.re - im * r.im;
        auto n_im = im*r.re + re*r.im;
        return RComplex(n_re, n_im);
    }
    RComplex operator*(Rational l)const {
        return RComplex(l*re, l*im);
    }
    RComplex operator/(Rational l)const {
        return RComplex(re/l, im/l);
    }

    RComplex inv()const {
        auto m = Rational(1)/sqabs();
        return conj()*m;

    } 

    RComplex operator/(RComplex r)const {
        auto n_r = RComplex(re,im);
        n_r = n_r*r.inv();
        return n_r;
    }

    RComplex operator+(RComplex r)const {
        return RComplex(re+r.re, im + r.im);
    }

    RComplex operator-(RComplex r)const {
        return RComplex(re-r.re, im - r.im);
    }

    void operator+=(RComplex r){
        re+=r.re; im+=r.im;
    }
    void operator-=(RComplex r){
        re-=r.re; im-=r.im;
    }
    



    Rational real() const{
        return re;
    }
    Rational imag() const{
        return im;
    }

    void Print() const{
        cout << "(" << re.a << "/" << re.b << "," << im.a << "/" << im.b << ")";
    }

    bool IsZero() const{
        return ((re.a == 0) && (im.a == 0));
    }

    
};

typedef RComplex cnum;
/*
 *Possible optimisations, add setleading inside add and all loops
 *Vector addition quickly
 *Monomial multipliation instead of simple shift? really?
 */
auto max(int a, int b){
    if(a>b) return a;
    return b;
}

int IdxFromPair(const Monomial& monomial){
    auto d = monomial.first;
    auto n = monomial.second;
    return d*(d+1)/2 +n;
}
Monomial PairFromIdx(int idx){
    auto discriminant = 1+8*idx;
    int d = floor((-1. + sqrt(static_cast<double>(discriminant)))*0.5);
    int n = idx - d*(d+1)/2;
    return Monomial(d, n);
}
    

//at 0 is zero term, leading term is at end of array
struct Polynomial{

    int degree;
    int leading_term;
    vector<cnum> t;
        
    
    Polynomial(int d){
        degree = d;
        leading_term = -1;
        vector<cnum> t_((d+1)*(d+2)/2);
        t = t_;
    }

    //Monomial
    Polynomial(int d, int n, cnum c){
        degree = d;
        leading_term = IdxFromPair(Monomial(d, n));
        vector<cnum> t_((d+1)*(d+2)/2);
        t = t_;
        t[leading_term]  = c;
    }
    
    cnum Get(const Monomial& mon) const{
        auto idx = IdxFromPair(mon);
        return t[idx];
    }

    bool IsZero() const{return (leading_term == -1);}
    
    void SetLeading(){
        leading_term=-1;
        for(int i=0; i<t.size(); i++){
            if(!(t[i].IsZero())) leading_term = i;
        }
    }
    void Reduce(){
        if(leading_term == -1) return;
        auto lt = PairFromIdx(leading_term);
        auto d = lt.first;
        if(d==degree) return;
        t.erase(t.begin()+((d+1)*(d+2)/2), t.end());
        degree=d;
    }
    /*Define operations on Polynomials*/
    Polynomial operator+(const Polynomial& P2) const{
        
        auto d_2 =P2.degree;
        auto d_new = max(degree, d_2);
        auto P_new = Polynomial(d_new);
        for(int i =0; i < P_new.t.size(); i++){
            if(i < (P2.t).size()) P_new.t[i] += P2.t[i];
            if(i<t.size()) P_new.t[i] += t[i];
            if(!((P_new.t[i]).IsZero())) P_new.leading_term = i;
        }
        return P_new;
    } 
    

    Polynomial operator-(const Polynomial& P2) const{
        
        auto d_2 =P2.degree;
        auto d_new = max(degree, d_2);
        auto P_new = Polynomial(d_new);
        for(int i =0; i < P_new.t.size(); i++){
            if(i < (P2.t).size()) P_new.t[i] -= P2.t[i];
            if(i<t.size()) P_new.t[i] += t[i];
            if(!(P_new.t[i]).IsZero()) P_new.leading_term = i;
        }
        P_new.Reduce();
        return P_new;
    } 

    Polynomial operator*(const cnum c) const{
        auto P_new = Polynomial(degree);
        if(c.IsZero()) return P_new;
        for(int i=0; i <P_new.t.size(); i++) P_new.t[i] = c*t[i];
        P_new.leading_term = leading_term;
        return P_new;
    }

    Polynomial operator*(const Polynomial& P2) const{
        if (P2.IsZero() || IsZero()) return Polynomial(0);
        
        auto P_new = Polynomial(degree + P2.degree);
        for(int d_1 = 0; d_1 <= degree; d_1++){
            for(int d_2 = 0; d_2 <= P2.degree; d_2++){
                for(int n_1 = 0; n_1 <=d_1; n_1++){
                    for(int n_2=0; n_2<=d_2; n_2 ++){
                        auto n_new = n_1+n_2;
                        auto d_new = d_1+d_2;
                        auto t1 = Get(Monomial(d_1, n_1));
                        auto t2 = P2.Get(Monomial(d_2, n_2));
                        auto idx_new = IdxFromPair(Monomial(d_new, n_new));
                        P_new.t[idx_new] += t1*t2;
                        if(!((P_new.t[idx_new]).IsZero())) P_new.leading_term = idx_new; // A few redudant sets here.. just lazy
                    }
                }
            }
        }
        return P_new;
    }

    void Print() const{
        cout << "Printing polynomial: ";
        for(int i = 0u; i <= leading_term; i++){
            if(t[i].IsZero()) continue;
            auto [d, n] = PairFromIdx(i);
            t[i].Print();
            cout << "x^" << n << "y^" << (d-n) << " + ";
        }
        cout << endl;
    }

    void PrintReal() const{
        cout << "Re(P): ";
        for(int i = 0u; i <= leading_term; i++){
            if(t[i].IsZero()) continue;
            auto [d, n] = PairFromIdx(i);
            cout << (t[i].re).Floatize() << "x^" << n << "y^" << (d-n) << " + ";
        }
        cout << endl;
    }
    Polynomial Derive(int axis){
        auto new_P = Polynomial(degree-1);

        for(int i = 1; i<t.size(); i++){
            auto [d_, n_] = PairFromIdx(i);
            auto val = t[i];
            if( 
                d_==0 || 
                (axis+n_ == 0) || 
                ((axis!=0) && ((d_-n_)==0))
            ){continue;}
            else if(axis ==0){val = val * Rational(n_); d_ -= 1; n_-=1;}
            else{val = val * Rational((d_-n_)); d_ -= 1;}

            auto new_i = IdxFromPair(Monomial(d_, n_));
            new_P.t[new_i] = val;
        }
        new_P.Reduce();
        new_P.SetLeading();
        return new_P;
    }
};   

Polynomial Divide(const Polynomial& P1, const Polynomial& P2){
    
    auto P_a = P1;
    auto P_b = P2;
    if(P2.IsZero()){ cout<<"Division by zero.." << endl; return Polynomial(0) ;}
    auto res = Polynomial(0);
    while(P_a.leading_term >= P_b.leading_term){
        auto [d_1, n_1] = PairFromIdx(P_a.leading_term);
        auto [d_2, n_2] = PairFromIdx(P_b.leading_term);
        auto divides = (((n_1-n_2) >= 0) && ( (d_1-d_2) >= (n_1 - n_2)));
        //cout << n_1 << "et " << n_2 << "et les d " << d_1 << " " << d_2 << endl;
        //cout << n_1 - n_2 << " " << d_1 - d_2 << endl;
        if(!divides){
            cout << "Error not perfect division..." << endl; 
            cout << n_1 << " " << n_2 << " " << d_1 << " " << d_2 << endl;
            cout << "Tried dividing: " << endl;
            P_a.Print();
            P_b.Print();
            return Polynomial(0);
        }
        auto l_a = P_a.t[P_a.leading_term];
        auto l_b = P_b.t[P_b.leading_term];
        auto c = l_a/l_b;
        auto Mon = Polynomial(d_1-d_2, n_1 - n_2, c);
        P_a = P_a - Mon * P_b;
        res=res+Mon;
    }
    return res;
}

typedef vector<vector<Polynomial>> Matrix;

void Bareiss(Matrix& M){
    int n = M.size();
    for(int k = 0; k < n; k++){
        for(int i =k+1; i< n; i++){
            for(int j=k+1; j<n; j++){
                cout << i << " " << j << " " << k << endl;
                auto A = M[i][j]*M[k][k] - M[i][k]*M[k][j];
                auto B = Polynomial(0, 0, cnum(1, 0));
                if(k!=0) B = M[k-1][k-1];
                if(B.IsZero()) cout << "Divide by zero!!!!!!!! " << endl;
                M[i][j] = Divide(A, B);
            }
            M[i][k] = Polynomial(0);
        }
    }

}

auto f_i(const cnum k, const cnum offset){
    auto P_x = Polynomial(1, 1, cnum(1, 0));
    auto P_y = Polynomial(1, 0, cnum(2, 0)*k);
    auto P_c = Polynomial(0, 0, offset);
    return (P_x-P_y+P_c);
}

auto EvolveOffset(const vector<cnum>& offsets, const vector<cnum>& k_s, const Rational time){
    auto n_offsets = vector<cnum>{};
    for(int i=0; i<k_s.size(); i++){
        auto n_o = offsets[i] + cnum(12, 0)*k_s[i]*k_s[i]*time;
        n_offsets.emplace_back(n_o);
    }
    return n_offsets;
}
auto Hirota(vector<cnum> k_s, vector<cnum> offsets)->Matrix{
    Matrix M;
    int n=k_s.size();
    for(int i=0; i<n; i++){
        auto line = vector<Polynomial>{};
        for(int j=0; j<n; j++){
           if(i==j){
               line.emplace_back(f_i(k_s[i], offsets[i]));
               continue;
           }
           auto c = cnum(0, 1) /(k_s[j] - k_s[i]);
           auto P = Polynomial(0, 0, c);
           line.emplace_back(P);
        }
        M.emplace_back(line);
    }
    return M;
}


int main(){
    
    cout << "Reading parameters" << endl;
    std::string write_to = "coefficients.txt";
    std::string filename = "params.txt";
    fstream file;
    file.open(filename, ios::in);
    bigint k_re, k_im, o_re, o_im;
    double t;
    auto k = vector<cnum>{};
    auto offsets = vector<cnum>{};
    file >> t;
    cout << "Fetching initial data " << endl;
    while(file >> k_re >> k_im >> o_re >> o_im){
        auto n_k = cnum(k_re, k_im);
        auto n_o = cnum(o_re, o_im);
        k.emplace_back(n_k); k.emplace_back((n_k.conj()));
        offsets.emplace_back(n_o); offsets.emplace_back(n_o.conj());
    }
    int n = k.size();
    file.close();
    cout << "Creating list of times" << endl; 
    //Time evolve
    auto times = std::vector<Rational>{};
    for(int i = 0; i < 1; i++){
        auto t_ = Rational(0, 1) + Rational(i, 1);
        times.emplace_back(t_);
    }
    file.open(write_to, ios::out);
     
    for(const auto t: times){
        cout << "t = " << t.Floatize() << endl;
        auto offsets_t = EvolveOffset(offsets, k, t);
        cout << "Creating Matrix" << endl;
        //Create Matrix
        auto M = Hirota(k, offsets_t);

        cout << "Getting determinant" << endl;
        Bareiss(M);
        auto P = M[n-1][n-1];
        P.Print();
        cout << "Extracting final polynomial" << endl;
        //Now calculate final polynomials
        P.Print();
        auto P_x = P.Derive(0);
        auto P_xx = P_x.Derive(0);
        auto P_xy = P_x.Derive(1);
        auto P_y = P.Derive(1); 
        auto Phi_x = P_xx*P - P_x*P_x;
        auto Phi_y = P_xy*P - P_x*P_y;
    
         
        cout << "----------------------------" << endl;
        cout << "Phi_x: ";
        Phi_x.PrintReal();
         
        cout << "----------------------------" << endl;
        cout << "Phi_y: ";
        Phi_y.PrintReal();


        //Write Polynomials to file
        file <<  t.Floatize() << " ";
        for(const auto& coeff: P.t){
            file << (coeff.real()).Floatize() << " ";
        }
        file << endl;
        file << t.Floatize() << " ";
        for(const auto& coeff: P_x.t){
            file << (coeff.real().Floatize()) << " ";
        }
        file << endl;
    }
    file.close();
    




}


