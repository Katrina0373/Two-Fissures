#ifndef FIND_FISSURES_H
#define FIND_FISSURES_H
#include <complex>
#include <vector>
#include <math.h>
#include <armadillo>
#include <map>

using std::complex; using arma::mat; using arma::vec; using arma::cx_vec; using arma::cx_double;

class Fissures {

public:

	Fissures();
	Fissures(double l1, double d1, double l2, double d2, int N = 20, double k = 5, int p = 1);
	void set_parameters(double l1, double d1, double l2, double d2);

	//���������� ������� ����������
	cx_vec solve_xi();
	//���������� ���� �������� � �������� �����
	cx_double number_field(const double x);
	//���������� ��������� ��������
	void eval_static_vecs();
	//������ � ��������� ���� ������ ���� ��������
	void write_field_to_csv(std::string file_name);
	//������ � ��������� ���� ������ ������� ����������
	void write_bundle_fun_to_csv(std::string file_name);

private:
	double l1, d1, l2, d2;
	int p;
	size_t N;
	vec t;
	vec t0;
	static double k30, k;
	static cx_vec roots_sigma2;
	cx_vec xi;
	static cx_vec sigma2_alpha;
	//��� ������� Areg
	static std::vector<double> k3_integral;
	double alphaM = 350;


	//==========================������ ������============================

	void check_parameters();
	void fill_t();
	void fill_Ac(mat& A11, mat& A21, mat& A12);
	void fill_Areg(mat& A11, mat& A21, mat& A12, mat& A22);
	void fill_F(cx_vec& F1, cx_vec& F2);
	void fill_k3_integral();
	double k3reg(const double x1, const double x2);
	double k3(const double alpha);

	double solve_k30();

	//==========================�������� ������============================
	
	static cx_vec sigma2, sigma1, u1, u2;
	static size_t N1;

	void fill_u_sigma();
	void fill_sigma2_alpha();
	cx_double L_x3_rj(const double r);

};


#endif // !FIND_FISSURES_H