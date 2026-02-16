#include "pointevolution.hpp"

typedef std::pair<double, double> vec;

void Point::Write(std::fstream& file){
    file << t << ',' << x << ',' << y << std::endl;
}

Point operator+(const Point& p1, const Point& p2){
    return Point{p1.x+p2.x, p1.y + p2.y, p1.t + p2.t};
}
Point operator*(const Point& p1, const double l){
    return Point{p1.x*l, p1.y*l, p1.t*l};
}

vec operator+(const vec& v1, const vec& v2){
    return {v1.first+v2.first, v1.second+v2.second};
}
vec operator*(const vec& v, const double l){
    return {v.first*l, v.second*l};
}
vec F(const Point& p, const BasicPolynomial& Det, const BasicPolynomial& T_x, const BasicPolynomial& T_y){
    auto d = Det.Evaluate(p.x, p.y, p.t);
    auto t_x = T_x.Evaluate(p.x, p.y, p.t);
    auto t_y = T_y.Evaluate(p.x, p.y, p.t);
    return {t_x/(d+1e-8), t_y/(d+1e-8)};
}

void EvolvePoint(Point& p, const double timestep, const BasicPolynomial& Det, const BasicPolynomial& T_x, const BasicPolynomial& T_y){
    auto d = Det.Evaluate(p.x, p.y, p.t);
    double adaptative_timestep = timestep * std::min(1.0, std::abs(d)/1e-3);
    auto k_1 = F(p, Det, T_x, T_y);
    auto p2 = p + Point{k_1.first, k_1.second, 1.}*(adaptative_timestep/2.);
    auto k_2 = F(p2, Det, T_x, T_y);
    auto p3 = p + Point{k_2.first, k_2.second, 1.}*(adaptative_timestep/2.);
    auto k_3 = F(p3, Det, T_x, T_y);
    auto p4 = p + Point{k_3.first, k_3.second, 1.}*adaptative_timestep;
    auto k_4 = F(p4, Det, T_x, T_y);
    auto step = (k_1 + k_2*2. + k_3*2.+k_4)*(adaptative_timestep/6.);
    auto pt = Point{step.first, step.second, adaptative_timestep};
    p = p + pt;
}
