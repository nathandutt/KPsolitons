#pragma once
#include <Eigen/Dense>
using Column = Eigen::Matrix<double, 2, N>;  

static const int N = 2;

struct SavedState{
    double y;
    Column poles;

    SavedState(const double y, const Column& p) :
	y(y), poles(p) {}
};

struct Point{
    double x;
    double y;
};
