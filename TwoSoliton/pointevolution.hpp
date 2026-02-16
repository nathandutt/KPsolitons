#include "basicpolynomial.hpp"
#include <fstream>
struct Point{
    double x;
    double y;
    double t;

    void Write(std::fstream& file);
};
Point operator+(const Point& p1, const Point& p2);
Point operator*(const Point& p1, const double l);

void EvolvePoint(Point& p, const double timestep, const BasicPolynomial& Det, const BasicPolynomial& T_x, const BasicPolynomial& T_y);
