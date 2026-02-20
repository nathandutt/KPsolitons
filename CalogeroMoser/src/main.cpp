#include <fstream>
#include <iostream>
#include <vector>
#include <complex>
#include "integration.hpp"
#include "npy.hpp"
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <iomanip>
using namespace std;
/*
 *       Format of params.txt
 *       t_i  t_f  t_step
 *       y_i  y_f  y_step y_write_step
 *       k1_re  k1_im  o1_re  o1_im
 *       ...
 *       kn_re  kn_im  on_re  on_im
 */
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

vector<complex<double>> EvolveOffsets(double t,vector<complex<double>>& k_s, vector<complex<double>>& offsets){
    auto new_offsets = std::vector<complex<double>>{};
    if(k_s.size()!=offsets.size()) throw std::runtime_error("Nonequal k_s and offsets");
    for(unsigned int i = 0; i < k_s.size(); i++){
        auto new_o = offsets[i] - 12. * k_s[i]*k_s[i]*t;
        new_offsets.emplace_back(new_o);
    }
    return new_offsets;
}

void PrintPoints(const Points& p){
    cout << "Y = " << p.y << endl;
    for(const auto& pole : p.poles){
        cout << "Point: " << pole.pole << ", Velocity: " << pole.velocity << endl;
    }
}
// !!! We suppose params.txt is well formatted, see above
int main(){
    const string input_file = "params.txt";
    const string output_dir = "Output/";
    fstream file;
    auto k_s = vector<complex<double>>{};
    auto offsets = vector<complex<double>>{};
     
    file.open(input_file, ios::in);
    
    //Get integration parameters
    double t_i, t_f, t_step;
    file >> t_i >> t_f >> t_step;
    double y_i, y_f, y_step, y_write_step;
    file >> y_i >> y_f >>  y_step >> y_write_step;
    
    //Get soliton parameters
    ReadParameters(file, k_s, offsets);
    file.close();

    auto T = floor((t_f - t_i)/t_step);
    auto Y = floor((y_f - y_i)/y_step);
    int pole_number = k_s.size();

    cout << "So we have " << pole_number/2 << " solitons" << endl
         << "And we are looking at " << T << " different time values" << endl
         << "And also integrating on y from " << y_i << " to " << y_f << endl;
    for(int i = 0; i < T; i++){
        //Here we iterate over different t values, for each t value we will write a new .npy data file
        //formatted as a numpy array of shape (Y, pole_number, 2)
        
        double t = t_i + t_step * i;
        auto t_data = std::vector<double>{};
        ostringstream ss;
        ss << fixed << setprecision(3) << t;
        string write_path = output_dir + ss.str() + ".npy";

        //Evolve offsets and create initial points
        auto evolved_offsets = EvolveOffsets(t, k_s, offsets); 
        auto points = Points(y_i, k_s, evolved_offsets);
        PrintPoints(points);
        double y_write = y_i;
        double written_y = 0;
        auto curr_y = points.y;
        while(curr_y < y_f){
            if(curr_y >= y_write) {
                AddPoints(t_data, points);
                //PrintPoints(points);
                //cout << "Wrote t=" << t << " and y= " << curr_y << endl;
                y_write += y_write_step;
                written_y++;
            }
            points.Evolve(y_step);
            curr_y = points.y;
        }
        const std::vector<unsigned long> shape{written_y, pole_number, 2};
        const npy::npy_data_ptr<double> data_ptr{t_data.data(), shape, false};
        write_npy(write_path, data_ptr);
    }
    
}
