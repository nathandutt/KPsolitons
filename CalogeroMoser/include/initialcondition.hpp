#pragma once
#include <complex>
#include <vector>

auto InitialConditions(const std::vector<std::complex<double>>& k_s, const std::vector<std::complex<double>>& offsets, const double y_i)
    ->std::pair<std::vector<std::complex<double>>, std::vector<std::complex<double>>>;


