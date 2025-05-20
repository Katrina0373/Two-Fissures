#ifndef SOLVE_EQS_H
#define SOLVE_EQS_H

#include<vector>
#include<cmath>
#include<complex>
#include<armadillo>;

using arma::cx_double;

std::vector<double> find_ratio_roots(const double k, double start, double end,
    const double h, const double eps, const double max_step);

arma::cx_vec find_complex_roots(double start, const double end, const double eps, const double k, int max_steps = 200, double h = 0.1);

double RungeKutta1(const double alpha, const double x0, const double U1_0, const double Sigma1_0, const double x_end, const double k, double h = 1e-5);
cx_double RungeKutta3(const double x0, const double U2_0, const double sigma2_0, const double U2_alpha_0, const double sigma2_alpha_0, const double x_end, const double h, const arma::cx_double alpha, const double k);

std::pair<cx_double, cx_double> RungeKuttaI(const cx_double alpha, const double x0, const double U2_0, const double sigma2_0, const double x_end, const double k, double h = 1e-5);


#endif // !
