#include "solve_eqs.h"

using namespace std;
using namespace arma;

double find_kernel_dihotomii(const double aa_r, const double bb_r, const double aa_i, const double bb_i,
    const double eps, const int max_step, const double k)
{
    int s = 0;
    double f_a, f_b, f_xm, xm_r, aaa_r, bbb_r, xm_i, aaa_i, bbb_i;
    f_a = RungeKutta1(aa_r, 0, 0, 1, 1, k);
    f_b = RungeKutta1(bb_r, 0, 0, 1, 1, k);
    aaa_r = aa_r;
    bbb_r = bb_r;
    aaa_i = aa_i;
    bbb_i = bb_i;
    while (eps < (double)(bbb_r - aaa_r) || eps < (double)(bbb_i - aaa_i))
    {
        xm_r = (double)(aaa_r + bbb_r) / 0.2e1;
        xm_i = (double)(aaa_i + bbb_i) / 0.2e1;
        f_xm = RungeKutta1(xm_r, 0, 0, 1, 1, k);
        if ((double)(f_a * f_xm) <= 0.0e0)
        {
            bbb_r = xm_r;
            bbb_i = xm_i;
            f_b = f_xm;
        }
        else
        {
            aaa_r = xm_r;
            aaa_i = xm_i;
            f_a = f_xm;
        }
    }
    return (double)(aaa_r + bbb_r) / 0.2e1;
}

vector<double> find_ratio_roots(const double k, const  double start, const double end,
    const double h, const double eps, const double max_step)
{
    vector<double> roots;
    double x = start;
    double prev_fun = RungeKutta1(x, 0, 0, 1, 1, k);
    while (x < end) {
        x += h;
        double cur_sigma = RungeKutta1(x, 0, 0, 1, 1, k);
        if (cur_sigma * prev_fun <= 0) {
            roots.push_back(find_kernel_dihotomii(x - h, x, 0, 0, eps, max_step, k));
        }
        prev_fun = cur_sigma;
    }

    return roots;
}

// Функция mu(x)
double mu(const double x) {
    return 1.0 + x * x / 10.0;
}

// Производная Sigma1
double dSigma(const double alpha, const  double k, const  double U1, const double mu) {
    return (alpha * alpha * mu - k * k) * U1;
}

// Производная U1
double dU(const double Sigma1, const double mu) {
    return Sigma1 / mu;
}


// Метод Рунге-Кутта 4-го порядка для решения системы ОДУ
double RungeKutta1(const double alpha, const double x0, const double U1_0, const double Sigma1_0, const double x_end, const double k, const double h) {
    double U1 = U1_0;
    double Sigma1 = Sigma1_0;
    double x = x0;
    double k1, k2, k3, k4, l1, l2, l3, l4;

    while (x < x_end) {
        
        k1 = h * dU(Sigma1, mu(x));
        l1 = h * dSigma(alpha, k, U1, mu(x));

        k2 = h * dU(Sigma1 + 0.5 * l1, mu(x + 0.5 * h));
        l2 = h * dSigma(alpha, k, U1 + 0.5 * k1, mu(x + 0.5 * h));

        k3 = h * dU(Sigma1 + 0.5 * l2, mu(x + 0.5 * h));
        l3 = h * dSigma(alpha, k, U1 + 0.5 * k2, mu(x + 0.5 * h));

        k4 = h * dU(Sigma1 + l3, mu(x + h));
        l4 = h * dSigma(alpha, k, U1 + k3, mu(x + h));

        U1 += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        Sigma1 += (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0;

        if (abs(x - 0.5) < 0.001)
            int num = 10;

        x += h;
    }

    return Sigma1;
}



//====================================  МНИМОЕ  ================================

// Производная Σ2 по x
cx_double dSigma2_i(const cx_double alpha, const double k, const cx_double U2, const double mu) {
    return (alpha * alpha * mu - k * k) * U2;
}

// Производная U2 по x
cx_double dU2_i(const cx_double sigma2, const double mu) {
    return sigma2 / mu;
}

pair<cx_double, cx_double> RungeKuttaI(const cx_double alpha, const double x0, const double U2_0, const double sigma2_0, const double x_end, const double k, const double h) {

    cx_double U2 = U2_0;
    cx_double sigma2 = sigma2_0;
    cx_double k1, k2, k3, k4, l1, l2, l3, l4;

    for (double x = x0; x < x_end; x+=h) {
        double mu_value = mu(x);

        // Для U2 и Σ2
        l1 = h * dSigma2_i(alpha, k, U2, mu_value);
        k1 = h * dU2_i(sigma2, mu_value);

        k2 = h * dU2_i(sigma2 + 0.5 * l1, mu_value);
        l2 = h * dSigma2_i(alpha, k, U2 + 0.5 * k1, mu_value);

        k3 = h * dU2_i(sigma2 + 0.5 * l2, mu_value);
        l3 = h * dSigma2_i(alpha, k, U2 + 0.5 * k2, mu_value);

        k4 = h * dU2_i(sigma2 + l3, mu_value);
        l4 = h * dSigma2_i(alpha, k, U2 + k3, mu_value);

        U2 += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        sigma2 += (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0;
    }

    return { sigma2, U2 };
}

// Производная Σ2,α по x
cx_double dSigma2_alpha(const cx_double alpha, const double mu, const double k, const cx_double U2, const cx_double U2_alpha) {
    return 2.0 * alpha * mu * U2 + (alpha * alpha * mu - k * k) * U2_alpha;
}

cx_double RungeKutta3(const double x0, const double U2_0, const double sigma2_0, const double U2_alpha_0, const double sigma2_alpha_0, const double x_end, const double h, const cx_double alpha, const double k) {
    cx_double k1, k2, k3, k4, l1, l2, l3, l4;
    cx_double m1, m2, m3, m4, n1, n2, n3, n4;

    cx_double U2 = U2_0;
    cx_double sigma2 = sigma2_0;
    cx_double U2_alpha = U2_alpha_0;
    cx_double sigma2_alpha = sigma2_alpha_0;

    for (double x = x0; x < x_end; x+=h) {
        double mu_value = mu(x);

        // Для U2 и Σ2
        k1 = h * dU2_i(sigma2, mu_value);
        l1 = h * dSigma2_i(alpha, k, U2, mu_value);

        k2 = h * dU2_i(sigma2 + 0.5 * l1, mu_value);
        l2 = h * dSigma2_i(alpha, k, U2 + 0.5 * k1, mu_value);

        k3 = h * dU2_i(sigma2 + 0.5 * l2, mu_value);
        l3 = h * dSigma2_i(alpha, k, U2 + 0.5 * k2, mu_value);

        k4 = h * dU2_i(sigma2 + l3, mu_value);
        l4 = h * dSigma2_i(alpha, k, U2 + k3, mu_value);

        sigma2 += (l1 + 2.0 * l2 + 2.0 * l3 + l4) / 6.0;

        // Для U2,α и Σ2,α
        m1 = h * dU2_i(sigma2_alpha, mu_value);
        n1 = h * dSigma2_alpha(alpha, mu_value, k, U2, U2_alpha);

        m2 = h * dU2_i(sigma2_alpha + 0.5 * n1, mu_value);
        n2 = h * dSigma2_alpha(alpha, mu_value, k, U2 + 0.5 * k1, U2_alpha + 0.5 * m1);

        m3 = h * dU2_i(sigma2_alpha + 0.5 * n2, mu_value);
        n3 = h * dSigma2_alpha(alpha, mu_value, k, U2 + 0.5 * k2, U2_alpha + 0.5 * m2);

        m4 = h * dU2_i(sigma2_alpha + n3, mu_value);
        n4 = h * dSigma2_alpha(alpha, mu_value, k, U2 + k3, U2_alpha + m3);

        U2_alpha += (m1 + 2.0 * m2 + 2.0 * m3 + m4) / 6.0;
        sigma2_alpha += (n1 + 2.0 * n2 + 2.0 * n3 + n4) / 6.0;
        U2 += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }

    return sigma2_alpha;
}


cx_double newton_method(const double x, const cx_double alpha0, const double eps, const double k, int max_iter = 100, double h = 0.001) {
    cx_double alpha = alpha0;

    for (int i = 0; i < max_iter; i++) {
        cx_double f_val = RungeKuttaI(alpha, 0.0, 0.0, 1.0, x, k, h).first;
        cx_double df_val = (RungeKuttaI(alpha + cx_double(0, 1e-8), 0.0, 0.0, 1.0, x, k, h).first - f_val) / cx_double(0, 1e-8);
        cx_double delta_alpha = f_val / df_val;

        alpha -= delta_alpha;

        if (abs(delta_alpha) < eps) {
            return alpha;
        }
    }
    throw runtime_error("Метод Ньютона не сошелся");
}

cx_vec find_complex_roots(const double start, const double end, const double eps, const double k, const int max_steps, const double h)
{
    cx_vec roots(50); size_t ind = 0;
    cx_double x = cx_double(0, start);
    cx_double prev_fun = RungeKuttaI(x, 0, 0, 1, 1, k).first;
    while (x.imag() < end) {
        x += cx_double(0,h);
        cx_double cur_fun = RungeKuttaI(x, 0, 0, 1, 1, k).first;
        if (cur_fun.real() * prev_fun.real() <= 0) {
            roots(ind) = newton_method(1, x, eps, k);
            ind++;
        }
        prev_fun = cur_fun;
    }
    return roots.head(ind);

}