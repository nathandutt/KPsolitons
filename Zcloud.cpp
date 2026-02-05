#include <iostream>
#include <stdio.h>
#include <vector>
#include <complex>
#include <string>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <cmath>
#include <random>
using namespace std;
const int N = 4; //Must be even, to count for complex conjugates

/*
 * Column struct allows us to define
 * addition and multiplication by a
 * scalar to better implement RK4.
*/

struct Column{

    vector<complex<double>> Z;
    
    //Constructors
    Column(){
        Z = std::vector<complex<double>>{};
    }
    Column(vector<complex<double>> Z_){
        Z = Z_;
    }

    /*
     * Below we define additon and
     * multiplication by a scalar
     * on columns
     */

    Column operator+(const Column& C2) const{
        auto new_Z = vector<complex<double>>{};
        for(int i=0; i<N; i++){
            auto val = C2.Z[i]+ Z[i];
            new_Z.emplace_back(val);
            
        }
        return Column(new_Z);
    }


    Column operator*(const double scalar) const {
        auto new_Z = vector<complex<double>>{};
        for(int i=0; i<N; i++){
            auto val = Z[i] * scalar;
            new_Z.emplace_back(val);
            
        }
        return Column(new_Z);
    }
    
	void Print() const{
	    cout << "Printing column: ";
	    for(int i= 0; i < N; i++){
	        cout << Z[i] << ", " ;
	    }
	    cout << endl;
	} //Print function for debugging puproses


};

//Inter-Pole interaction force
auto F(const Column& C)-> Column{

    auto res = vector<complex<double>>{};
	for(int i =0; i < N; i++){
	    complex<double> F_i = 0.;
	    for(int j=0; j<N; j++){
	        if(i==j) continue;
	        F_i += 2./pow((C.Z[i] - C.Z[j]), 3.);
	    }
	    res.emplace_back(F_i);
	}
	return Column(res);
}

double GreatestVariation(const Column& F, const Column& C_dot){

    auto var = 0.;
    for(int i=0; i<N; i++){
        auto local_var = fabs(F.Z[i]/(C_dot.Z[i]+1e-5));
        if(local_var > var) var = local_var;
    }
    return var;
}//Used for variable timestep


/*
 *Define a State struct, for easy integration
 *through RK4
 */


struct PoleState{

    Column C;
    Column C_dot;
    double time;

    PoleState(Column C_, Column C_dot_, double time_){

        C = C_;
        C_dot = C_dot_;
        time = time_;
    }

    /*
     *Update Method implements the RK4 integration.
     *We use a variable timestep, because our force
     *has large singularities. Timestep is adapted
     *so that variation of our state stays beneath
     *a certain threshold percentage.
     */

	void Update(const double eps, const double percentage){

	    // k1
	    Column k1_C     = C_dot;
	    Column k1_Cdot  = F(C);
        
        //Change timestep if variation is too large
        auto var = GreatestVariation(k1_Cdot, C_dot);
        auto var_2 = GreatestVariation(C_dot, C);
        auto tstep = min(percentage/var, eps);
        tstep = min(percentage/var_2, tstep);

	    // k2
	    Column C2       = C     + k1_C    * (0.5 * tstep);
	    Column Cdot2    = C_dot + k1_Cdot * (0.5 * tstep);
	    Column k2_C     = Cdot2;
	    Column k2_Cdot  = F(C2);
	
	    // k3
	    Column C3       = C     + k2_C    * (0.5 * tstep);
	    Column Cdot3    = C_dot + k2_Cdot * (0.5 * tstep);
	    Column k3_C     = Cdot3;
	    Column k3_Cdot  = F(C3);
	
	    // k4
	    Column C4       = C     + k3_C    * tstep;
	    Column Cdot4    = C_dot + k3_Cdot * tstep;
	    Column k4_C     = Cdot4;
	    Column k4_Cdot  = F(C4);
	
	    // combine
	    C     = C     + (k1_C    + k2_C*2.0 + k3_C*2.0 + k4_C)    * (tstep / 6.0);
	    C_dot = C_dot + (k1_Cdot + k2_Cdot*2.0 + k3_Cdot*2.0 + k4_Cdot) * (tstep / 6.0);
	
	    time += tstep;
	}
    
    //For writing to CSV file
    void Write(const std::string& filename) const {
        std::ofstream out(filename, std::ios::app);
        if (!out) {
            throw std::runtime_error("Could not open file: " + filename);
        }

        out << time;
        for (const auto& c : C.Z) {
            out << "," << c.real() << "," << c.imag();
        }
        out << "\n";
        out.close();
    }
};


/*
 *Our aim is to study a Soliton Gas. Therefore we
 *need to be able to generate random ensembles of 
 *solitons. Which is equivalent to generating the k's
 *we only take N/2 because, the other half are the
 *complex conjugates.
 */

//! Recursive function, loops till it has generated a correct value.
//In a circle distribution around z_avg
auto GenerateComplex(uniform_real_distribution<>& d, mt19937 g, double r, complex<double>& z_avg)-> complex<double>{
    complex<double> res;
    auto x = r*d(g);
    auto y = r*d(g);
    cout << "x, y " << x << " " << y << endl;
    auto r_h = sqrt(pow(x,2.)+ pow(y, 2.));
    if(r_h > r) res = GenerateComplex(d,g, r, z_avg);
    else{
        res = complex<double>(x, y) + z_avg;
    }
    return res;

}

//Generate N (global const) poles distributed in a circle of radius r around z_avg
auto GeneratePoles(uniform_real_distribution<>& d,mt19937 g, double y_init, double r, complex<double> z_avg)-> PoleState{
    auto C = vector<complex<double>>{};
    auto C_dot = vector<complex<double>>{};
    for(int i=0; i < N/2; i++){
        cout << "i = " << i << endl;
        auto b_k = GenerateComplex(d, g, r, z_avg);
        auto z_k = 2.*b_k*y_init;
        auto b_k_bar = 2.*conj(b_k);
        auto z_k_bar = 2.*conj(z_k);

        C.emplace_back(z_k_bar);
        C.emplace_back(z_k);
        C_dot.emplace_back(b_k_bar);
        C_dot.emplace_back(b_k);
    }
    return PoleState(C, C_dot, y_init);
}

//Generate an initial pole state given a vector of k values
// !Careful! : For now, you must set N (global const) by hand to be equal to 2*k_values.size()
auto PolesFromK(double y_init, const vector<complex<double>>& k_values, const vector<complex<double>>& offsets){
    auto C = vector<complex<double>>{};
    auto C_dot = vector<complex<double>>{};
    for(int i =0; i<static_cast<int>(k_values.size()); i++){
        auto z_val_1 = 2.*conj(k_values[i])*y_init+offsets[i];
        auto z_val_2 = 2.*k_values[i]*y_init+conj(offsets[i]);
        auto z_dot_1 = 2.*conj(k_values[i]);
        auto z_dot_2 = 2.*(k_values[i]);
        C.emplace_back(z_val_1);
        C.emplace_back(z_val_2);
        C_dot.emplace_back(z_dot_1);
        C_dot.emplace_back(z_dot_2);
    }
    return PoleState(C, C_dot, y_init);
}

//Parameters for Run
const auto y_init = -100.;
const auto y_end = 100.;
const auto eps = 1e-2;
const auto steps = fabs(y_end - y_init) / eps;
const auto write_every = 0.2;
const auto writing_step = 20;
const auto max_percentage = 0.1;
int main(){
    /*
     *cout << "Setting up random seed..." << endl;
	 *std::random_device rd;  // Will be used to obtain a seed for the random number engine
	 *std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
	 *std::uniform_real_distribution<> dis(-1., 1.0);
	 *cout << "Done." << endl;
	 *cout << "Generating Random Solitons..." << endl;
	 *auto Poles = GeneratePoles(dis, gen, y_init, 0.1, 1. + 1.*1i);
    */

	auto filename = "pole_evolution.csv";
    //Values of k
    auto k_vals = vector<complex<double>>{
        complex<double>(3., 4.),
        complex<double>(1., 2.),
    };
    //Offsets
    auto offsets = vector<complex<double>>{
        complex<double>(-10., 0.),
        complex<double>(0., 0.),
    };
    
    //Create PoleState
    cout << "Creating PoleState " << endl;
    auto Poles = PolesFromK(y_init, k_vals, offsets);
    cout << "Done." << endl;

    cout << "Starting Numerical Integration..." << endl;
    cout << "Params" << endl << "-----------------------" << endl;
    cout << "Total estimated steps " << steps << " with Timestep " << eps << endl << "We write every " << write_every << " time interval"<< endl;
    cout << "Max percentage increase is set to " << max_percentage << endl;
    cout << "-----------------------"<< endl;
    double write_time = y_init;
    int curr_step = 0;
    while(Poles.time < y_end){
        if(write_time <= Poles.time) {Poles.Write(filename); write_time +=write_every;}
        /*if(curr_step%(writing_step)==0) {
            Poles.Write(filename); 
            write_time +=write_every;
            //cout << "Step " << curr_step << endl;
        }*/
        Poles.Update(eps, max_percentage);
        curr_step++;
    }
    cout << "Finished Simulation! Total steps: " << curr_step << ", Predicted was " << steps << endl;
    cout << endl;
    return 0;
}
