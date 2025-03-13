#ifndef FIND_FISSURES_H
#define FIND_FISSURES_H
#include <complex>
#include <vector>
#include <math.h>
#include <armadillo>
#include <omp.h>
#include <map>

using std::complex; using arma::mat; using arma::vec; using arma::cx_vec; using arma::cx_double;

class Fissures {

public:
	double l1, d1, l2, d2, k30, k;
	int N, p;
	vec t;
	vec t0;
	cx_vec roots_sigma2, xi;
	//äëÿ ìàòðèöû Areg
	static std::vector<double> k3_integral;
	double alphaM;

	Fissures();
	Fissures(double l1, double l2, double d1, double d2, int N = 20, double k = 5, int p = 1);
	void set_coordinates(double l1, double d1, double l2, double d2);

	//==========================ÏÐßÌÀß ÇÀÄÀ×À============================

	void fill_t();
	void fill_Ac(mat& A11, mat& A21, mat& A12);
	void fill_Areg(mat& A11, mat& A21, mat& A12, mat& A22);
	void fill_F(cx_vec& F1, cx_vec& F2);
	void fill_k3_integral();
	double k3reg(const double x1, const double x2);
	double k3reg1(const double x1,const double x2, vec k3_integral,const double widht);
	double k3(const double alpha);

	cx_vec solve_xi();
	double solve_k30();

	//==========================ÎÁÐÀÒÍÀß ÇÀÄÀ×À============================
	
	cx_vec sigma2, sigma1, u1, u2;
	cx_vec sigma2_alpha;
	size_t N1;

	cx_double number_field(const double x);
	void fill_u_sigma();
	cx_double alpha_x3_rj(const cx_double r);
};


#endif // !FIND_FISSURES_H