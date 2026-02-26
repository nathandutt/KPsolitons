#include "polestate.hpp"

Column F(const Column& m, double& max_c){
    //Takes a column c of poles, and returns a column of forces
    // where F_i_x = \sum_{j \neq i} r_ij^-6 * (x_ij^3 - 3x_ijy_ij^2)
    // and F_i_y = \sum_{j\neqi} r_ij^-6 * (x_ij^3 - 3x_ij^2y_ij)

    //We want to take advantage of vectorization possibilities
    
    Column f;

    Eigen::Array<double, 1, N> r_squared =
        m.row(0).array().square() + m.row(1).array().square();

    if(max_c < 0.)
	max_c = r_squared.maxCoeff();

    Eigen::Array<double, 1, N> r_cubed_inv = r_squared.pow(-1.5);

    Eigen::Array<double, 1, N> num0 =
        m.row(0).array().cube() - 3 * m.row(0).array() * m.row(1).array().square();

    Eigen::Array<double, 1, N> num1 =
        m.row(1).array().cube() - 3 * m.row(1).array() * m.row(0).array().square();

    Eigen::Array<double, 2, N> numerators;
    numerators.row(0) = r_cubed_inv * num0;
    numerators.row(1) = r_cubed_inv * num1;

    for (int i = 0; i < N; ++i)
    {
        Eigen::Array<bool, 1, N> mask = Eigen::Array<bool, 1, N>::Constant(N, true);
        mask(i) = false;

        f(0, i) = (numerators.row(0) * mask.cast<double>()).sum();
        f(1, i) = (numerators.row(1) * mask.cast<double>()).sum();
    }

    return f;
}

constexpr double adaptative_pow = 0.8;
void PoleState::Evolve(const double timestep){
    //RK4 evolution with adaptative timestep
    double max_c = -1.;

    //first k
    Column k1_p = velocity;
    Column k1_v = F(poles, max_c);
    //Define adaptative timestep
    double astep = (max_c > 1.) ? timestep : timestep*std::pow(max_c, adaptative_pow*0.5);

    Column k2_p = velocity + (0.5*astep)*k1_v;
    Column k2_v = F(poles + (0.5*astep)*k1_p, max_c);

    Column k3_p = velocity + (0.5*astep)*k2_v;
    Column k3_v = F(poles+(0.5*astep)*k2_p, max_c);

    Column k4_p = velocity + astep*k3_v;
    Column k4_v = F(poles + astep*k3_p, max_c);

    velocity += astep/6. * (k1_v + 2.*k2_v + 2.*k3_v + k4_v);
    poles += astep/6. * (k1_p + 2.*k2_p + 2.*k3_p + k4_p);
    y += astep;
}
