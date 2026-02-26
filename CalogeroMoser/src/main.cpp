#include "tracker.hpp"
#include "global.hpp"
#include "polestate.hpp"
#include "initialconditions.hpp"
#include "field.hpp"
#include "npy.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <complex>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace std;
/*
 *       Format of params.txt
 *       t_i  t_f  t_step
 *       y_i  y_f  y_step y_write_step
 *       x_i  x_f  x_step 
 *       k1_re  k1_im  o1_re  o1_im
 *       ...
 *       kn_re  kn_im  on_re  on_im
 */

//Read Evolution parameters
struct Config{
    double y_i; double y_f; double y_step; double y_write_step;
    double x_i; double x_f; double x_step;
    double t_i; double t_f; double t_step;

    Config(fstream& file){
	file >> t_i >> t_f >> t_step;
	file >> y_i >> y_f >> y_step >> y_write_step;
	file >> x_i >> x_f >> x_step;
    }
    void Print(){
	cout << "Y parameters are: " << endl
	     << "Initial y: " << y_i << " Final y: " << y_f <<
	     " Pole integration step: " << y_step << " Zero finding step: " << y_write_step << endl
	     << "X parameters are: " << endl
	     << "Initial x: " << x_i << " Final x: " << x_f << " Zero finding step: " << x_step << endl
	     << "T parameters are: " << endl
	     << "Initial t: " << t_i << " Final t: " << t_f << " Timestep:" << t_step << endl;
    }
};

//Read soliton parameters
void ReadParameters(fstream& file, vector<complex<double>>& k_s, vector<complex<double>>& offsets){
    double k_re, o_re, k_im, o_im;
    while(file >> k_re >> k_im >> o_re >> o_im){
       auto k = complex<double>(k_re, k_im); 
       auto o = complex<double>(o_re, o_im); 
       k_s.emplace_back(k);
       offsets.emplace_back(o);
       offsets.emplace_back(conj(o));
       k_s.emplace_back(conj(k));
    }
}

//Evolve soliton offsets for a fixed time
vector<complex<double>> EvolveOffsets(double t,const vector<complex<double>>& k_s, const vector<complex<double>>& offsets){
    auto new_offsets = std::vector<complex<double>>{};
    if(k_s.size()!=offsets.size()) throw std::runtime_error("Nonequal k_s and offsets");
    for(unsigned int i = 0; i < k_s.size(); i++){
        auto new_o = offsets[i] - 12. * k_s[i]*k_s[i]*t;
        new_offsets.emplace_back(new_o);
    }
    return new_offsets;
}


//Given a fixed time, and soliton parameters, get all critical points
std::vector<Point> GetZeros(const double time, const Config& config, const vector<complex<double>>& k_s, const vector<complex<double>>& offsets){
    //First, evolve soliton parameters
    auto zero_list = std::vector<Point>{}; 
    auto t_offsets = EvolveOffsets(time, k_s, offsets);

    //Initialize poles
    auto y_curr = config.y_i;
    auto points = PoleSet(config.y_i, k_s, t_offsets);
    auto y_write = config.y_i;

    //Initialize Zero Trackers;
    auto trackers = TrackerList{};
    auto written = 0;
    while(y_curr < config.y_f){

	//If at zero searching resolution, look for zeros and add them to trackers
	int write_every = floor(config.y_write_step/config.y_step);
	if(written%write_every==0){
	    y_write+=config.y_write_step;

	    //Get zeros for fixed y
	    auto zeros = FindZero(points, config.x_i, config.x_f, config.x_step);
	    //Cycle through zeros and add them to trakcer
	    for(const auto& x_coord : zeros){
		auto new_pt = Point(x_coord, y_curr);
		auto new_val = Phi_y(points, x_coord);
		TrackPoint(trackers, new_pt, new_val, 0.3*config.y_write_step);
//		zero_list.emplace_back(new_pt);
	    }
	}

	//Evolve points with integration step
	points.Evolve(config.y_step);
	y_curr = points.y;
	written++;
    }
    return ExtractZeros(trackers);
}

// !!! We suppose params.txt is well formatted, see above
int main(){
    const string input_file = "params.txt";
    const string output_dir = "Output/";
    fstream file;
    auto k_s = vector<complex<double>>{};
    auto offsets = vector<complex<double>>{};
     
    file.open(input_file, ios::in);
    //Get config
    auto config = Config(file);

    cout << "Take a moment to check parameters please. " << endl;
    config.Print();

    //Get soliton parameters
    ReadParameters(file, k_s, offsets);
    file.close();

    int pole_number = k_s.size();
    assert(pole_number == N);
    
    //for now only simulate one time
    auto current_time = config.t_i;
    auto t_offsets = EvolveOffsets;


    //Initialize poles
    auto [p, v] = InitialConditions(config.y_i, k_s, t_offsets);
    PoleState poles(config.y_i, p, v);
    
    //Vector to save
    int total_pt_estimate = floor((config.y_f - config.y_i)/config.y_write_step);
    std::vector<std::unique_ptr<SavedState>> v;
    v.reserve(1.5*total_pt_estimate);

    //Initialize loop parameters

    double current_y = config.y_i;
    double next_write_y = config.y_i;

    while(current_y < config.y_f){
	if(current_y >= next_write_y){
	    poles.Insert(v);
	    next_write_y = config.y_write_step;
	}
	poles.Evolve(config.y_step);
	current_y += config.y_step;
    }
    







//    auto current_time = config.t_i;
//    while(current_time < config.t_f){
//
//	//Define write path for .npy output
//	ostringstream ss;
//	ss << fixed << setprecision(5) << current_time;
//	string write_path = output_dir + ss.str() + ".npy";
//
//	//Get critical points
//	cout << "Time = " << current_time << endl;
//	auto critical_points = GetZeros(current_time, config, k_s, offsets);
//	cout << "Got " << critical_points.size() << " critical points" <<endl;
//	current_time += config.t_step;
//	cout << "Writing" << endl;    
//	//Put data in correct format for .npy writing
//	const std::vector<unsigned long> shape{critical_points.size(), 2};
//	auto data = std::vector<double>{};
//	//write points to data vector
//	for(const auto& pt : critical_points){
//	    data.emplace_back(pt.first); data.emplace_back(pt.second);
//	}
//	//initialize ptr
//	const npy::npy_data_ptr<double> data_ptr{data.data(), shape, false};
//	//write to .npy
//	write_npy(write_path, data_ptr);
//    }
}
