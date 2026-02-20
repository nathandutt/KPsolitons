#include <complex>
#include <vector>

struct Pole{
    std::complex<double> pole;
    std::complex<double> velocity;

    Pole(const std::complex<double>& pole, const std::complex<double>& velocity):
        pole(pole), velocity(velocity) {};
};
struct Points{
    double y;
    std::vector<Pole> poles;

    Points(const double y_init, const std::vector<std::complex<double>>& k_s, const std::vector<std::complex<double>>& offsets);
    Points(){
        y = 1.;
        poles = std::vector<Pole>{};
    }
    void Evolve(const double timestep); 
};

//For now we don't need to write time idx, we just want to plot trajectories
inline void AddPoints(std::vector<double>& data, const Points& points){
    for(auto pole : points.poles){
        data.emplace_back(pole.pole.real());
        data.emplace_back(pole.pole.imag());
    }
}
    
