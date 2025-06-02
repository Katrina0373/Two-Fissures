#include"find_fissure.h"
#include"solve_eqs.h"
#include<iostream>
#include<chrono>

using std::complex; using std::vector; using std::cout; using std::endl;
using namespace arma;
using namespace std::complex_literals;

std::vector<double> Fissures::k3_integral;
arma::cx_vec Fissures::sigma2_alpha;
arma::cx_vec Fissures::roots_sigma2;
arma::cx_vec Fissures::sigma1, Fissures::sigma2, Fissures::u1, Fissures::u2;
double Fissures::k30;
size_t Fissures::N1;
double (*mu)(double) = { [](double x) {return
		1 + std::pow(x, 2) / 10; } };



void Fissures::fill_Ac(mat& A11, mat& A21, mat& A12) {
	
	A11.clear();
	A21.clear();
	A12.clear();
	A11.resize(N, N);
	A21.resize(N, N);
	A12.resize(N, N);

	for (size_t i = 0; i < N; i++)
	{ 
			for (int j = 0; j < N; j++)
			{
				A11(i, j) = 2 * k30 * (1 / (t[i + 1] - t0[j]) - 1 / (t[i] - t0[j]));

				A21(i, j) = 2 * k30 * (
					1 / (d2 + l2 * t[i + 1] - d1 - l1 * t0[j]) -
					1 / (d2 + l2 * t[i] - d1 - l1 * t0[j]));

				A12(i, j) = 2 * k30 * (
					1 / (d1 + l1 * t[i + 1] - d2 - l2 * t0[j]) -
					1 / (d1 + l1 * t[i] - d2 - l2 * t0[j]));

			}
	}
}

void Fissures::fill_Areg(mat& A11, mat& A21, mat& A12, mat& A22) {
	
	
	A11.clear();
	A21.clear();
	A12.clear();
	A22.clear();
	A11.resize(N, N);
	A21.resize(N, N);
	A12.resize(N, N);
	A22.resize(N, N);

	for (size_t i = 0; i < N; i++)
	{
			for (int j = 0; j < N; j++)
			{
				A11(i, j) = k3reg(
					l1 * (t[i + 1] - t0[j]),
					l1 * (t[i] - t0[j]));

				A21(i, j) = k3reg(
					d2 + l2 * t[i + 1] - d1 - l1 * t0[j],
					d2 + l2 * t[i] - d1 - l1 * t0[j]);

				A12(i, j) = k3reg(
					d1 + l1 * t[i + 1] - d2 - l2 * t0[j],
					d1 + l1 * t[i] - d2 - l2 * t0[j]);

				A22(i, j) = k3reg(
					l2 * (t[i + 1] - t0[j]),
					l2 * (t[i] - t0[j]));
			}
	}
}



double Fissures::k3reg(const double x1, const double x2) {
	
	double h = alphaM / k3_integral.size();
	double integral = 0;
	auto Dij = [x1, x2, this](double alpha, size_t i) {
		return k3_integral[i] *	(sin(alpha * x1) - sin(alpha * x2));
	};

	double step = h / 2;
	for (size_t i = 0; i < k3_integral.size(); i++)
	{
		integral += Dij(step, i);
		step += h;
	}

	//cout << 2 * h * integral << endl;
	return 2 * h * integral;
}

double Fissures::k3reg1(const double y1, const double y2) {

	double width = alphaM / k3_integral.size();

	auto Dij = [y1, y2, this](double alpha, int i) {
		
		return k3_integral[i] * (sin(alpha * y1) - sin(alpha * y2));
	};

	double simpson_integral = 0;
	
	int start = 0;
	double x1 = width / 2 + start * width;
	double x2 = width / 2 + (start + 1) * width;
	double x3 = (x2 - x2) / 2;

	for (int step = start; step < k3_integral.size(); step++) {
		simpson_integral += width / 6.0 * (Dij(x1, step) + 4.0 * Dij(x3, step) + Dij(x2, step + 1));
		x1 += width;
		x2 += width;
		x3 += width;
	}

	//cout << simpson_integral << endl;
	return simpson_integral;
}

//x = 1
double Fissures::k3(const double alpha) {
	double sigma1 = RungeKutta1(alpha, 0, 1, 0, 1, k, 1e-5);
	double sigma2 = RungeKutta1(alpha, 0, 0, 1, 1, k, 1e-5);
	return -sigma1 / sigma2;
}

void Fissures::fill_F(cx_vec& F1, cx_vec& F2) {

	F1.clear();
	F2.clear();
	F1.resize(N);
	F2.resize(N);

	auto F = [this](double x) {
		cx_double sum = 0;
		int sign = x > 0 ? 1 : -1;
		
		for (size_t i = 0; i < N1; i++)
		{
			if (sign == 1)
				sum += exp(roots_sigma2(i) * x * 1.i) / sigma2_alpha(N1 + i);
			else 
				sum += exp(roots_sigma2(i) * -x * 1.i) / sigma2_alpha(i);
			//cout << sum << endl;
		}
		
		return sign * p * 2 * datum::pi * sum * 1.i;
	};


	for (size_t i = 0; i < N; i++)
	{
		F1(i) = F(d1 + l1 * t0[i]);
		F2(i) = F(d2 + l2 * t0[i]);
	}

 }

void Fissures::fill_k3_integral()
{
	double h = 0.1;
	size_t steps = alphaM / h;
	k3_integral.resize(steps);
	
	double alpha = h / 2;

	for (int i = 0; i < steps; i++) //2646?
	{
		k3_integral[i] = (k3(alpha) - k30 * alpha) / alpha;
		alpha += h;
	}

}

bool is_symmetric(const arma::cx_mat& A, double tol = 1e-10) {
	return arma::approx_equal(A, A.t(), "absdiff", tol);
}

cx_vec Fissures::solve_xi()
{
	mat Ac11, Ac12, Ac21;
	cx_vec F1, F2;
	mat Areg11, Areg12, Areg21, Areg22;

	//auto start_time = std::chrono::high_resolution_clock::now();
	fill_Ac(Ac11, Ac21, Ac12);
	/*auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end_time - start_time;
	std::cout << "Время Ac: " << duration.count() << std::endl;
	start_time = std::chrono::high_resolution_clock::now();*/
	
	
	fill_Areg(Areg11, Areg12, Areg21, Areg22);
	
	//end_time = std::chrono::high_resolution_clock::now();
	//duration = end_time - start_time;
	//std::cout << "Время Areg: " << duration.count() << std::endl;

	//start_time = std::chrono::high_resolution_clock::now();	
	/*end_time = std::chrono::high_resolution_clock::now();
	duration = end_time - start_time;
	std::cout << "Время Areg: " << duration.count() << std::endl;

	start_time = std::chrono::high_resolution_clock::now();
	*/
	fill_F(F1, F2);

	/*end_time = std::chrono::high_resolution_clock::now();
	duration = end_time - start_time;
	std::cout << "Время F: " << duration.count() << std::endl << std::endl;*/

	cx_mat A(2*N, 2*N, fill::zeros);
	cx_vec B(2*N);

	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			A(i, j) = Ac11(j, i) / l1 + l1 * Areg11(j, i);
			A(i, N + j) = l2 * Ac21(j, i) + l2 * Areg21(j, i);
			A(i + N, j) = l1 * Ac12(j, i) + l1 * Areg12(j, i);
			A(i + N, j + N) = Ac11(j, i) / l2 + l2 * Areg22(j, i);
		}
		B(i) = F1(i);
		B(i + N) = F2(i);
	}
	//cout << B << endl << endl;
	
	xi = solve(A, B);
	
	//cout << xi << endl;
	return xi;
}

double Fissures::solve_k30() {
	return k3(600) / 600;
}

cx_double Fissures::number_field(const double x)
{
	cx_double sum1 = 0, sum2 = 0;
	for (size_t i = 0; i < N; i++)
	{
		double r = d1 + l1 * t0[i] - x;
		sum1 += xi(i) * L_x3_rj(r);
		r = d2 + l2 * t0[i] - x;
		sum2 += xi(N + i) * L_x3_rj(r);
	}
	return (l1 * sum1 + l2 * sum2) / (datum::pi * N);
}

//учитывая что количество корней счётно и достаточно мало, а обращение к задачам Коши будет частым, 
// то лучше создать под них собственные массивы
void Fissures::fill_u_sigma()
{
	sigma1.resize(2 * N1);
	sigma2.resize(2 * N1);
	u1.resize(2 * N1);
	u2.resize(2 * N1);
	for (size_t i = 0; i < N1; i++)
	{
		//заполняем значения sigma1 u1 
		auto res_cauchy = RungeKuttaI(roots_sigma2(i), 0, 1, 0, 1.0, k);
		sigma1(i) = res_cauchy.first;
		u1(i) = res_cauchy.second;
		res_cauchy = RungeKuttaI(-roots_sigma2(i), 0, 1, 0, 1.0, k);
		sigma1(i + N1) = res_cauchy.first;
		u1(i + N1) = res_cauchy.second;

		//заполняем значения sigma2 u2
		res_cauchy = RungeKuttaI(roots_sigma2(i), 0, 0, 1, 1.0, k);
		sigma2(i) = res_cauchy.first;
		u2(i) = res_cauchy.second;
		res_cauchy = RungeKuttaI(-roots_sigma2(i), 0, 0, 1, 1.0, k);
		sigma2(i + N1) = res_cauchy.first;
		u2(i + N1) = res_cauchy.second;
	}

}
//для всех моделей с одинаковым k значения sigma2_alpha одни и те же
void Fissures::fill_sigma2_alpha() {
	vec roots_sigma2re(find_ratio_roots(k, 0, 100, 0.1, 1.0e-8, 1000));
	roots_sigma2.resize(roots_sigma2re.size());
	for (size_t i = 0; i < roots_sigma2re.size(); i++)
	{
		roots_sigma2(i) = cx_double(roots_sigma2re[i], 0);
	}
	roots_sigma2 = join_cols(roots_sigma2, find_complex_roots(0, 100, 1.0e-8, k, 300, 0.1));
	N1 = roots_sigma2.size();

	sigma2_alpha.resize(2 * N1);
	for (int i = 0; i < N1; i++)
	{
		sigma2_alpha(i) = RungeKutta3(0, 0, 1, 0, 0, 1, 1.0e-5, roots_sigma2[i], k);
		sigma2_alpha(N1 + i) = RungeKutta3(0, 0, 1, 0, 0, 1, 1.0e-5, -roots_sigma2[i], k);
	}
}

cx_double Fissures::L_x3_rj(const double r)
{
	cx_double sum = 0;
	if (r > 0) {
		for (size_t i = 0; i < N1; i++)
		{
			sum += (sigma2(i) * u1(i) - sigma1(i) * u2(i)) 
				/ sigma2_alpha(i) * exp(1.i * roots_sigma2(i) * r);
		}
	}
	else {
		for (size_t i = 0; i < N1; i++)
		{
			sum += (sigma2(N1 + i) * u1(N1 + i) - sigma1(N1 + i) * u2(N1 + i)) 
				/ sigma2_alpha(N1 + i) * exp(1.i * roots_sigma2(i) * -r);
		}
		sum *= -1;
	}
	return 2 * datum::pi * sum * 1.i;
}

void Fissures::eval_static_vecs()
{
	k30 = solve_k30();
	fill_k3_integral();
	fill_sigma2_alpha();
	fill_u_sigma();
}

void Fissures::check_parameters()
{
	double eps = 1e-7;
	if (l1 < eps || l2 < eps) {
		throw("The half-length of the bundle is too small");
	}
}

void Fissures::fill_t()
{
	t = linspace(-1, 1, N+1);
	t0 = linspace((t[0] + t[1]) / 2, (t[N] + t[N - 1]) / 2, N);
}

Fissures::Fissures()
{
	l1 = 0.1; l2 = 0.1;
	d1 = 0.5; d2 = 1;
	N = 20; k = 5; p = 1;
	fill_t();
}

Fissures::Fissures(double l1, double d1, double l2, double d2, int N, double k, int p)
{
	if (N == 0)
		throw("The number of nodes must be greater than 0");
	this->l1 = l1;
	this->l2 = l2;
	this->d1 = d1;
	this->d2 = d2;
	this->k = k;
	this->N = N;
	this->p = p;
	check_parameters();
	fill_t();
}

void Fissures::set_parameters(double l1, double d1, double l2, double d2)
{
	this->l1 = l1;
	this->l2 = l2;
	this->d1 = d1;
	this->d2 = d2;
	check_parameters();
}

double Fissures::get_k()
{
	return k;
}
double Fissures::get_l1()
{
	return l1;
}
double Fissures::get_d1()
{
	return d1;
}
double Fissures::get_l2()
{
	return l2;
}
double Fissures::get_d2()
{
	return d2;
}
double Fissures::get_N()
{
	return N;
}
double Fissures::get_p()
{
	return p;
}
