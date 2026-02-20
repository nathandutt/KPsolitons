#include "integration.hpp"
#include <stdexcept>
#include <iostream>

Points::Points(double y_init, const std::vector<std::complex<double>>& k_s, const std::vector<std::complex<double>>& offsets){
    if(k_s.size() != offsets.size()) throw std::runtime_error("Not as many ks as offsets");
    y = y_init;
    poles = std::vector<Pole>{};
    for(unsigned int i = 0; i < k_s.size(); i++){
        auto new_pole = Pole{
           2.*k_s[i]*y_init + offsets[i],
           2.*k_s[i]
        };
        poles.emplace_back(new_pole);
    }
}

//TODO Make this more efficient by spatial selection if i want many points
//for now it is O(N^2)
//Okay for because i'm looking at only 4 points
//Plus all these operations on Points could be sped up if many points

Points Combination(const Points& p1, const Points& p2, const double coeff){ //returns p1 + coeff*p2
    auto res = Points();
    if(p1.poles.size() != p2.poles.size()) 
        throw std::runtime_error("Summng two point list with non equal number of points");
    
    for(unsigned int i = 0; i < p1.poles.size(); i++){
        const auto& pole_1 = p1.poles[i];
        const auto& pole_2 = p2.poles[i];
        res.poles.emplace_back(
                Pole{
                    pole_1.pole + coeff*pole_2.pole,
                    pole_1.velocity + coeff*pole_2.velocity
                });
    }
    res.y = p1.y + coeff*p2.y;
    return res;
}

Points F(const Points& points, double& max_f){
    auto res = Points();
    for(unsigned int i = 0; i < points.poles.size(); i++){
        res.poles.emplace_back(Pole(points.poles[i].velocity, 0.));
    }
    
    for(unsigned int i = 0; i<points.poles.size(); i++){
        for(unsigned int j = 0; j<i; j++){

            //Derivative
            if(std::fabs(points.poles[i].pole - points.poles[j].pole) < 1e-2) throw std::runtime_error("Numerical problem");
            auto f_ij = 1./std::pow((points.poles[i].pole - points.poles[j].pole), 3);
            res.poles[i].velocity += f_ij;
            res.poles[j].velocity -= f_ij;
            if(std::abs(f_ij) > max_f) max_f = std::abs(f_ij);
        }
    }
    return res;
}
double PowerOfTen(double k){
    int a = std::floor(std::log(k)/std::log(10.));
    a++;
    return std::pow(10., -2*a);
}
//RK4 evolution
void Points::Evolve(const double timestep){
    const auto& p1 = *(this);
    double max_f = 0.;
    auto k1 = F(p1, max_f);
    auto adaptative_timestep = timestep;
    auto l = std::pow(max_f, 1.);
    if(max_f > 1.) adaptative_timestep = timestep/l;
    auto p2 = Combination(p1, k1, adaptative_timestep/2.);
    auto k2 = F(p2, max_f);
    auto p3 = Combination(p1, k2, adaptative_timestep/2.);
    auto k3 = F(p3, max_f);
    auto p4 = Combination(p1, k3, adaptative_timestep);
    auto k4 = F(p4, max_f);
    auto res = Points();
    for(unsigned int i =0; i<p1.poles.size(); i++){
        const auto& pole1 = k1.poles[i];
        const auto& pole2 = k2.poles[i];
        const auto& pole3 = k3.poles[i];
        const auto& pole4 = k4.poles[i];
        
        poles[i].pole += adaptative_timestep/6. * (pole1.pole + 2.*pole2.pole + 2.*pole3.pole + pole4.pole);
        poles[i].velocity += adaptative_timestep/6. * (pole1.velocity + 2.*pole2.velocity + 2.*pole3.velocity + pole4.velocity);
    }
    y+=adaptative_timestep;
}
