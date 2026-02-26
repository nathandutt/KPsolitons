#include <complex>
#include <vector>

struct Pole{
    std::complex<double> pole;
    std::complex<double> velocity;

    Pole(const std::complex<double>& pole, const std::complex<double>& velocity):
        pole(pole), velocity(velocity) {};
};
struct PoleSet{
    double y;
    std::vector<Pole> poles;

    PoleSet(const double y_init, const std::vector<std::complex<double>>& k_s, const std::vector<std::complex<double>>& offsets);
    PoleSet(){
        y = 1.;
        poles = std::vector<Pole>{};
    }
    double Evolve(const double timestep); 
};

//For now we don't need to write time idx, we just want to plot trajectories
inline void AddPoleSet(std::vector<double>& data, const PoleSet& points){
    for(auto pole : points.poles){
        data.emplace_back(pole.pole.real());
        data.emplace_back(pole.pole.imag());
    }
}
    
