#include"Algorithms.h"
#include<iostream>
#include<random>
using std::vector;
//=======================Вспомогательные функции===========================

//меняет местами значение концов отрезка, если надо
void a_b_correction(double& a, double& b) {
	if (a > b) {
		std::swap(a, b);
	}
}

//выбор начальных точек (метод многогранника)
vector<vector<double>> choose_x0(vector<double>& x0, double& l) {
	int n = x0.size();
	double s, r;
	s = l * (sqrt(n + 1) - 1 + n) / (sqrt(2) * n);
	r = l * (sqrt(n + 1) - 1) / (sqrt(2) * n);

	vector<vector<double>> x(n + 1);
	for (size_t i = 0; i < n + 1; i++)
		x[i].resize(n);

	x[0] = x0;
	for (size_t i = 1; i < n + 1; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			if ((j + 1) == i)
				x[i][j] = x0[j] + r;
			else
				x[i][j] = x0[j] + s;
		}
	}

	return x;
}

//упорядочивает точки начиная с номера s (метод многогранника)
void arrange(int s, vector<vector<double>>& x, vector<double>& fi) {

	int n = x.size();

	for (size_t i = s; i < n; i++)
	{
		for (size_t j = i; j > s - 1; j--)
		{
			if (fi[j - 1] > fi[j]) {
				std::swap(x[j - 1], x[j]);
				std::swap(fi[j - 1], fi[j]);
			}
		}
	}


}

vector<double> gradient(vector<double> x, double(*f)(vector<double>), double delta) {
	int n = x.size();
	vector<double> y = x, g(n);
	double f_y, f_x;
	f_x = f(x);
	for (size_t i = 0; i < n; i++)
	{
		y[i] += delta;
		f_y = f(y);
		g[i] = (f_y - f_x) / delta;
		y[i] = x[i];
	}
	return g;
}

//вычисление Евклидовой нормы вектора
double Euclidean_norm(vector<double> v) {
	double norm = 0;
	for (size_t i = 0; i < v.size(); i++)
		norm += pow(abs(v[i]), 2);
	return sqrt(norm);
}

//скалярное произведение векторов
double scalar(vector<double> x1, vector<double> x2) {
	return std::inner_product(std::begin(x1), std::end(x1), std::begin(x2), 0.0);
}

//номер строки наибольшего числа в стобце [col]  
int col_max(const vector<vector<double>>& matrix, int col, int n) {
	double max = std::abs(matrix[col][col]);
	int maxPos = col;
	for (int i = col + 1; i < n; ++i) {
		double element = abs(matrix[i][col]);
		if (element > max) {
			max = element;
			maxPos = i;
		}
	}
	return maxPos;
}

int triangulation(vector<vector<double>>& matrix, int n) {
	unsigned int swapCount = 0;
	if (0 == n)
		return swapCount;
	const int num_cols = matrix[0].size();
	for (int i = 0; i < n - 1; ++i) {
		unsigned int imax = col_max(matrix, i, n);
		if (i != imax) {
			swap(matrix[i], matrix[imax]);
			++swapCount;
		}
		for (int j = i + 1; j < n; ++j) {
			double mul = -matrix[j][i] / matrix[i][i];
			for (int k = i; k < num_cols; ++k) {
				matrix[j][k] += matrix[i][k] * mul;
			}
		}
	}
	return swapCount;
}

//Вычисление определителя по Гауссу
double gauss_determinant(std::vector<std::vector<double> >& matrix, int n) {
	unsigned int swapCount = triangulation(matrix, n);
	double determinanit = 1;
	if (swapCount % 2 == 1)
		determinanit = -1;
	for (int i = 0; i < n; ++i) {
		determinanit *= matrix[i][i];
	}
	return determinanit;
}



void composition(std::vector<double>& i, std::vector<double>& j, double k) { //  i = j *k +i
	if (i.size() != j.size()) {
		throw "Ошибка! Несоответсивие размеров двух векторов!";
	}
	for (int d = 0; d < i.size(); d++) {
		i[d] = i[d] + j[d] * k;
	}
}
void divStr(std::vector<double>& d, double k) {
	if (abs(k) < 0.00000001) {
		throw "Ошибка! На 0 делить нельзя!";
	}

	for (int i = 0; i < d.size(); i++) {
		d[i] = d[i] / k;
	}
}
void  solve_LES(std::vector<std::vector<double> >& a) {
	for (int k = 0; k < a.size(); k++) {
		if (abs(a[k][k]) < 0.00000001) {
			for (int i = k + 1; i < a.size(); i++) {
				if (a[i][k] > 0.00000001) {
					std::swap(a[i], a[k]);
					break;
				}
			}
		}

		if (a[k][k] < 0.999999999 || a[k][k] > 1.000000001) {
			divStr(a[k], a[k][k]);
		}

		for (int i = k + 1; i < a.size(); i++) {
			composition(a[i], a[k], -a[i][k]);
		}
	}

	for (int k = a.size() - 1; k >= 0; k--) {
		for (int i = k - 1; i >= 0; i--) {
			composition(a[i], a[k], -a[i][k]);
		}
	}
}

//подходит только для гессиана
bool is_positive_H(vector<vector<double>> matr) {

	if (matr[0][0] <= 0)
		return false;

	size_t n = matr.size();

	double det;
	for (size_t i = 1; i < n; i++)
	{	
		det = gauss_determinant(matr, i);
		if (det <= 0)
			return false;
	}
	return true;
}



//f - указателей на функцию, для которой вычисляется гессиан
//args - точка, в которой вычисляется гессиан
//step - шаг метода вычисления
vector<vector<double>> ComputeHessian(double (*f)(vector<double>), vector<double> args, double step)
{
	
	int varCount = args.size();
	vector<vector<double>> hessian(varCount);
	for (size_t i = 0; i < varCount; i++)
		hessian[i].resize(varCount);

	//Вспомогательные 
	vector<double> fPlus(varCount);
	vector<double> fMinus(varCount);
	vector<double> args1 = args;


	double curFValue = f(args); //значение функции в текущей точке 
	for (int i = 0; i < varCount; i++)
	{
		/* В случае производной по одной переменной берём в качестве аппроксимации вторую центральную производную */
		args1[i] += step;
		fPlus[i] = f(args1);
		args1[i] = args[i];
		args1[i] -= step;
		fMinus[i] = f(args1);
		args1[i] = args[i];

		hessian[i][i] = (fPlus[i] - 2 * curFValue + fMinus[i]) / (step * step);
	}

	for (int i = 0; i < varCount; i++)
		//Так как гессиан симметричен, вычисляем только его половину
		for (int j = i + 1; j < varCount; j++)
		{
			args1[i] += step;
			args1[j] += step;
			hessian[i][j] = (f(args1) - fPlus[i] - fPlus[j] + curFValue) / (step * step);
			args1[i] = args[i];
			args1[j] = args[j];
		}

	for (int i = 0; i < varCount; i++)
		for (int j = 0; j < i; j++)
			hessian[i][j] = hessian[j][i];


	return hessian;
}

//======================Методы минимизации функции=========================

std::pair<double, double> local_segment_Svenn(double x0, double delta, double (*f)(double))
{
	double a = 0, b = 0, x, y;
	if (f(x0) <= f(x0 + delta) && f(x0) <= f(x0 - delta)) {
		a = x0 - delta;
		b = x0 + delta;
		return std::make_pair(a, b);
	}
	else if (f(x0) >= f(x0 + delta) && f(x0) >= f(x0 - delta))
		return std::make_pair(a, b);

	if (f(x0) <= f(x0 - delta) && f(x0) >= f(x0 + delta)) {
		a = x0; x = x0; y = x0 + delta;
	}
	else {
		b = x0; x = x0; y = x0 - delta; delta *= -1;
	}

	do {
		delta *= 2; x = y; y = x + delta;
		if (f(x) <= f(y)) {
			if (delta > 0)
				b = y;
			else
				a = y;
		}
		else {
			if (delta > 0)
				a = x;
			else
				b = x;
		}
	} while (f(x) > f(y));

	return { a, b };
}

//NOTE: не работает на функции pow(8.7 * pow(sin(x),2) + 1.1, 1.0/3) + 4.0 * exp(x - 3.4)  [-5.8, -2.1]
double passive_serach(double a, double b, double(*f)(double), double eps)
{
	a_b_correction(a, b);
	int n = 1;
	while ((b - a) / (n + 1) >= eps / 2)
		n += 2;

	double h = (b - a) / (n + 1); //шаг
	int k = 0;					//количество проб
	double f1, f2 = DBL_MAX;
	double c = a - 2 * h, d = b, x = a;
	do {
		k++;
		x += h;
		f1 = f2;
		f2 = f(x);
		c += h;
		if (f1 <= f2) d = x;
	} while (k != n && f1 > f2);

	return (c + d) / 2;
}

double dichotomy(double a, double b, double(*f)(double), double eps, double delta)
{
	a_b_correction(a, b);
	if (eps < delta) {
		delta = eps / 10;
	}
	double c = a, d = b, x;
	while (abs(d - c) >= eps) {
		x = (d + c - delta) / 2;
		if (f(x) <= f(x + delta))
			d = x + delta;
		else c = x;
	}
	return (c + d) / 2;
}

double dichotomy_2(double a, double b, std::function<double(double)> f, double eps, double delta)
{
	a_b_correction(a, b);
	if (eps < delta) {
		delta = eps / 10;
	}
	double c = a, d = b, x;
	while (abs(d - c) >= eps) {
		x = (d + c - delta) / 2;
		if (f(x) <= f(x + delta))
			d = x + delta;
		else c = x;
	}
	return (c + d) / 2;
}


double divide_by_half(double a, double b, double(*f)(double), double eps)
{
	a_b_correction(a, b);
	double c = a, d = b, x1 = (c + d) / 2, x2, x3;

	do
	{
		x2 = (c + x1) / 2;
		x3 = (d + x1) / 2;

		if (f(x1) >= f(x2)) {
			d = x1; x1 = x2;
		}
		else if (f(x1) >= f(x3)) {
			c = x1; x1 = x3;
		}
		else {
			c = x2; d = x3;
		}

	} while (d - c >= eps);

	return x1;
}

double Fibonachi_method(double a, double b, double(*f)(double), double eps, double delta)
{
	a_b_correction(a, b);
	int n = 2;
	double L;
	std::vector<int> F(3);
	F[0] = 1, F[1] = 1, F[2] = 2;
	L = (b - a) / F[2] + (F[0] / F[2]) * delta;
	while (L >= eps) {
		int t = F[2];
		F[2] = F[2] + F[1];
		F[0] = F[1];
		F[1] = t;
		n++;
		L = (b - a) / F[2] + (F[0] / F[2]) * delta;
	}

	double x1 = b - (F[1] * (b - a) / F[2] + pow(-1, n) / F[2] * delta);
	double x2 = a + b - x1;

	for (size_t i = 2; i < n; i++)
	{
		if (f(x1) <= f(x2)) {
			b = x2; x2 = x1; x1 = a + b - x2;
		}
		else {
			a = x1; x1 = x2; x2 = a + b - x1;
		}
	}
	return (a + b) / 2;
}

double gold_sector(double a, double b, double(*f)(double), double eps)
{
	a_b_correction(a, b);
	double tau = 1.618034;
	double x1, x2, c = a, d = b;
	x1 = b - 1 / tau * (b - a);
	x2 = a + b - x1;

	while (d - c >= eps) {
		if (f(x1) <= f(x2)) {
			d = x2;
			x2 = x1;
			x1 = c + d - x2;
		}
		else {
			c = x1;
			x1 = x2;
			x2 = c + d - x1;
		}
	}
	return (x2 + x1) / 2;
}

//не трогать, а то сломается
double quadr_interpolation(double a, double b, double(*f)(double), double eps, double delta)
{
	a_b_correction(a, b);
	bool ok = true;
	double x1, x2, x3, _x, x_min, f1, f2, f3, f_min, _f;
	do
	{
		if (ok) {
			x2 = (b + a) / 2;
			x1 = x2 - delta;
			x3 = x2 + delta;
		}
		f1 = f(x1); f2 = f(x2); f3 = f(x3);

		//нахождение минимума из трёх
		if (f1 > f2) {
			if (f2 > f3) {
				f_min = f3;
				x_min = x3;
			}
			else {
				f_min = f2;
				x_min = x2;
			}
		}
		else {
			if (f1 > f3) {
				f_min = f3;
				x_min = x3;
			}
			else {
				f_min = f1;
				x_min = x1;
			}
		}

		if ((f1 - f3) * (x2 - x3) == (f2 - f3) * (x1 - x3)) {  // три точки на одной прямой 
			if (f_min < f3)
				b = x2;
			else
				a = x2;

			_x = x2; _f = f2; ok = true;
		}
		else {
			_x = (f1 * (x3 * x3 - x2 * x2) + f2 * (x1 * x1 - x3 * x3) + f3 * (x2 * x2 - x1 * x1))
				/ (2 * (f1 * (x3 - x2) + f2 * (x1 - x3) + f3 * (x2 - x1)));
			_f = f(_x);
			if (abs((f_min - _f) / _f) >= eps || abs((x_min - _x) / _x) >= eps) {
				if (_x >= x1 && _x <= x3) {
					if (_f <= f3)
						b = x3;
					if (f1 >= _f) {
						a = x1;
						if (f2 >= _f)
							a = x2;
					}
					if (_x < x2) {
						x3 = x2; x2 = _x;
					}
					else {
						x1 = x2; x2 = _x;
					}
					ok = false;
				}
				else {
					if (f1 <= f2)
						b = x2;
					else {
						a = x1;
						if (f2 >= f3)
							a = x2;
					}
					ok = true;
				}
			}
		}
	} while ((abs(f_min - _f)/_f) > eps && abs((x_min - _x) / _x) >= eps);

	return x_min;
}

//======================Методы минимизации n-го порядка=========================
vector<double> polyhedrom_method(double(*f)(vector<double>), double eps,
	vector<double> x0, double alfa, double beta, double gamma, double l, int max_iter) {

	int k = 0, s = 1, n = x0.size();
	double atc = 0.0;
	vector<double> c(n), u(n), w(n), v(n), fi(n+1);

	auto x = choose_x0(x0, l);
	for (size_t i = 0; i <= n; i++)
	{
		fi[i] = f(x[i]);
	}
	arrange(s, x, fi);
	int iter = 0;
	do {
		
		for (size_t i = 0; i < n; i++)
		{
			c[i] = std::accumulate(x[i].begin(), x[i].end() - 1, 0.0) / n;
		}


		for (size_t i = 0; i < n; i++)
		{
			u[i] = c[i] + alfa * (c[i] - x[n][i]);
		}

		s = n;

		double fu = f(u);

		if (fu < fi[0]) {
			for (size_t i = 0; i < n; i++)
			{
				v[i] = c[i] + beta * (u[i] - c[i]);
			}
			double fv = f(v);
			if (fv < fu) {
				x[n] = v;
				fi[n] = fv;
			}
			else {
				x[n] = u;
				fi[n] = fu;
			}
		}
		else {
			if (fi[0] <= fu && fu <= fi[n - 1]) {
				fi[n] = fu;
				x[n] = u;
			}
			else {
				if (fu <= fi[n]) {
					for (size_t i = 0; i < n; i++)
					{
						w[i] = c[i] + gamma * (u[i] - c[i]);
					}
				}
				else {
					for (size_t i = 0; i < n; i++)
					{
						w[i] = c[i] + gamma * (x[n][i] - c[i]);
					}
				}
				double fw = f(w);
				if (fw < fmin(fi[n], fu)) {
					fi[n] = fw;
					x[n] = w;
				}
				else {
					for (size_t i = 1; i < n + 1; i++)
					{
						for (size_t j = 0; j < n; j++)
						{
							x[i][j] = (x[0][j] + x[i][j]) / 2;
						}
					}
					s = 1;
				}
			}
		}
		k++;
		atc = 0;
		for (size_t i = 1; i <= n; i++)
		{
			atc += pow(f(x[i]) - f(x[0]), 2);
		}
		atc = sqrt(atc / n);

		arrange(s, x, fi);
		iter++;
	} while (atc >= eps && iter < max_iter);

	vector<double> res(n);
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n + 1; j++)
		{
			res[i] += x[j][i];
		}
		res[i] = res[i] / (n + 1);
	}

	return res;
}

std::vector<double> polyhedrom_method(double(*f)(std::vector<double>), double eps, std::vector<double> x0, double l)
{
	double alfa = 1, beta = 2, gamma = 0.5;
	return polyhedrom_method(f, eps, x0, alfa, beta, gamma, l);
}

//осторожно выбирать delta, желательно меньше 0.01
std::vector<double> gradient_method_of_steepest_descent(double(*f)(std::vector<double>),
	double eps, std::vector<double> x0, double eps1, double delta, double delta2) {

	int n = x0.size();
	vector<double> x2 = x0, x1, s1, arg1(n), arg2(n);
	double delta3 = 0.00001, a, b, lambda, lambda1;

	while (Euclidean_norm(gradient(x2, f, delta3)) >= eps1) {
		x1 = x2;
		s1 = gradient(x1, f, delta3);

		a = 0; lambda = delta; b = 2 * delta;
		for (size_t i = 0; i < n; i++)
		{
			s1[i] = -s1[i];
			arg1[i] = x1[i] + lambda * s1[i];
			arg2[i] = x1[i] + (lambda + delta) * s1[i];
		}

		while (f(arg1) > f(arg2)) {
			lambda += delta;
			delta *= 2;
			b += delta;
			for (size_t i = 0; i < n; i++)
			{
				arg1[i] = x1[i] + lambda * s1[i];
				arg2[i] = x1[i] + (lambda + delta) * s1[i];
			}
		}

		std::function<double(double)> fi = [f, s1, x1](double l) {
			vector<double> arg(x1.size());
			for (size_t i = 0; i < x1.size(); i++)
				arg[i] = x1[i] - l * s1[i];
			return f(arg); };

		lambda1 = dichotomy_2(a, b, fi, eps, delta2);

		for (size_t i = 0; i < n; i++)
		{
			x2[i] = x1[i] + lambda1 * s1[i];
		}

	}
	return x2;
}

//не чувствителен к параметрам; рекомендуется эпсилон = 0.1, лямбда = 1, бета = 0.95
std::vector<double> gradient_method_with_step_splitting(double(*f)(std::vector<double>), double eps, std::vector<double> x0, double eps1, double lambda, double beta)
{
	size_t n = x0.size();
	vector<double> x2 = x0, x1, s1, arg1(n);
	double delta = 0.00000001, lambda1, norm = 0;

	while (Euclidean_norm(gradient(x2, f, delta)) >= eps1) {
		x1 = x2;
		s1 = gradient(x1, f, delta);
		lambda1 = lambda;
		for (size_t i = 0; i < n; i++)
		{
			s1[i] = -s1[i];
			arg1[i] = x1[i] + lambda1 * s1[i];
		}
		double y = f(x1), z = f(arg1);
		norm = Euclidean_norm(s1);
		while ((z - y) > -eps * lambda1 * norm * norm) {
			lambda1 *= beta;
			for (size_t i = 0; i < n; i++)
			{
				arg1[i] = x1[i] + lambda1 * s1[i];
			}
			z = f(arg1);
		}

		for (size_t i = 0; i < n; i++)
			x2[i] = x1[i] + lambda1 * s1[i];
	}

	return x2;
}

std::vector<double> conjugate_gradient(double(*f)(std::vector<double>), double eps, std::vector<double> x0) {

	size_t n = x0.size();
	vector<double> x1(n), x2 = x0, x3 = x0, s;
	int k = 0;
	double beta, delta = 0.000001, lambda;

	while (Euclidean_norm(gradient(x2, f, delta)) >= eps) {
		x1 = x2;

		if (k % n == 0) {
			beta = 0;
			s = gradient(x1, f, delta);
			for (size_t i = 0; i < n; i++)
				s[i] = -s[i];
		}
		else {
			auto s1 = gradient(x1, f, delta);
			auto s3 = gradient(x3, f, delta);
			beta = scalar(s1, s1) / scalar(s3, s3);
			auto s_add = gradient(x1, f, delta);
			for (size_t i = 0; i < n; i++)
				s[i] = -s_add[i] + beta * s[i];
		}

		std::function<double(double)> fi = [f, s, x1](double l) {
			vector<double> arg(x1.size());
			for (size_t i = 0; i < x1.size(); i++)
				arg[i] = x1[i] - l * s[i];
			return f(arg); };

		lambda = dichotomy_2(0, 1, fi, eps, eps / 10);

		for (size_t i = 0; i < n; i++)
			x2[i] = x1[i] + lambda * s[i];
		x3 = x1;
		k++;
	}
	return x2;
}

std::vector<double> newton_method(double (*f)(std::vector<double>), double eps, std::vector<double> x0) {
	size_t n = x0.size();
	vector<double> x2 = x0, x1(n), s(n), s1(n);
	double delta = 0.0000001, lambda;
	while (Euclidean_norm(gradient(x2, f, delta)) >= eps)
	{
		x1 = x2;
		auto H = ComputeHessian(f, x1, delta);
		s1 = gradient(x1, f, delta);
		for (size_t i = 0; i < n; i++)
			s1[i] = -s1[i];
		if (is_positive_H(H)) {
		
			for (size_t i = 0; i < n; i++)
				H[i].push_back(s1[i]);

			solve_LES(H);

			for (size_t i = 0; i < n; i++)
				s[i] = H[i][n];
			
			lambda = 1;
		}
		else {
			s = s1;
			std::function<double(double)> fi = [f, s, x1](double l) {
				vector<double> arg(x1.size());
				for (size_t i = 0; i < x1.size(); i++)
					arg[i] = x1[i] - l * s[i];
				return f(arg); };

			lambda = dichotomy_2(0, 1, fi, eps, eps / 10);
		}
		for (size_t i = 0; i < n; i++)
		{
			x2[i] = x1[i] + lambda * s[i];
		}
	}
	return x2;
}

std::vector<double> newton_raphson_method(double(*f)(std::vector<double>), double eps, std::vector<double> x0)
{
	size_t n = x0.size();
	vector<double> x2 = x0, x1(n), s(n), s1(n);
	double delta = 0.00000001, lambda;
	while (Euclidean_norm(gradient(x2, f, delta)) >= eps)
	{
		x1 = x2;
		auto H = ComputeHessian(f, x1, delta);
		s1 = gradient(x1, f, delta);
		for (size_t i = 0; i < n; i++)
			s1[i] = -s1[i];
		if (is_positive_H(H)) {

			for (size_t i = 0; i < n; i++)
				H[i].push_back(s1[i]);
			solve_LES(H);

			for (size_t i = 0; i < n; i++)
				s[i] = H[i][n];

			std::function<double(double)> fi = [f, s, x1](double l) {
				vector<double> arg(x1.size());
				for (size_t i = 0; i < x1.size(); i++)
					arg[i] = x1[i] + l * s[i];
				return f(arg); };

			lambda = dichotomy_2(0, 1, fi, eps, eps / 10);

			
		}
		else {
			s = s1;
			std::function<double(double)> fi = [f, s, x1](double l) {
				vector<double> arg(x1.size());
				for (size_t i = 0; i < x1.size(); i++)
					arg[i] = x1[i] - l * s[i];
				return f(arg); };

			lambda = dichotomy_2(0, 1, fi, eps, eps / 10);
		}
		for (size_t i = 0; i < n; i++)
		{
			x2[i] = x1[i] + lambda * s[i];
		}
	}
	return x2;
}

//==================================================Генетический алгоритм====================================================


double random_state0(const double min, const double max) {
	std::random_device rd; // Источник случайных чисел
	std::mt19937 gen(rd()); // Генератор случайных чисел с инициализацией
	std::uniform_real_distribution<> dis(min, max); // Распределение на [min, max)
	return dis(gen);
}

//Мутация
//vector<double> mutation(const vector<double> x, const vector<double> x0, const vector<double> x1, int cur_generation, int max_generation) {
//	
//	vector<double> y(x0.size());
//	double max_delta, delta, mutation;
//	for (size_t i = 0; i < x0.size(); i++)
//	{
//		max_delta = (x1[i] - x0[i]) * 0.07;
//		delta = max_delta * (1.0 - (double)cur_generation / max_generation);
//		mutation = ((static_cast<double>(rand()) / RAND_MAX) * 2 - 1) * delta;
//		y[i] = x[i] + mutation;
//	}
//	return y;
//}

vector<double> mutation(const vector<double> point, const vector<double> x0, const vector<double> x1, const double fi)  // мутация: генерация случайной величины
{
	const int NUM = 100000000;
	vector<double> x(point);
	for (size_t i = 0; i < x0.size(); i++)
	{
		if (fi < 0.001 && i % 2 == 0)				//если функция уже меньше какого-то значения то возможно l у неё более менее точные но с d всегда проблема
			continue;
		x[i] = fabs((double)((rand() * NUM) % (int)((x1[i] - x0[i]) * NUM) + 1) / NUM) + x0[i];
	}
	return x;
}

vector<double> inversion(vector<double> x, const double eps, const double fi) {
	static int sign = 0;
	sign++;
	sign %= 2;
	if (sign == 0) {
		for (size_t i = 0; i < x.size(); i++) {
			if (fi < 0.01 && i % 2 == 0)
				continue;
			x[i] -= eps;
		}
	}
	else {
		for (size_t i = 0; i < x.size(); i++) {
			if (fi < 0.01 && i % 2 == 0)
				continue;
			x[i] += eps;
		}
	}
	return x;
}

vector<double> avg_gens(const vector<double> x1, const vector<double> x2) {
	vector<double> sum(x1.size());
	for (size_t i = 0; i < x1.size(); i++)
		sum[i] = (x1[i] + x2[i]) / 2;
	
	return sum;
}

vector<double> mutation_around(vector<double> x, const vector<double> x0, const vector<double> x1){
	for (size_t i = 0; i < x.size(); i++)
	{
		double eps = (x1[i] - x0[i]) / 5;
		x[i] = random_state0(x[i] - eps, x[i] + eps);
	}
	return x;
}

vector<double> avg_gen(const vector<double> x1, const vector<double > x2) {
	return vector<double>();
}

//Скрещивание и формирование новой популяции
vector<vector<double>> crossover(const vector<vector<double>> p, const double eps, const vector<double> x0, const vector<double> x1, const int iter, const int max_iter, const vector<double> fi)	{
	int k = p.size()-1, n = p.size();
	vector<vector<double>> new_p(n);

	for (int i = 0; i < 4; i++)
		for (int j = i + 1; j < 4; j++)
		{
			new_p[k] = avg_gens(p[i], p[j]);
			k--;
	 	}

	new_p[k] = p[0];
	k--; 
	new_p[k] = p[1];
	k--;
	for (size_t i = 0; i < 5; i++)
	{
		new_p[k] = mutation_around(p[0], x0, x1);			//сделаем мутации в окрестностях точки
		k--;
	}

	for (size_t i = 0; i < 5; i++)
	{
		new_p[k] = inversion(p[i], eps, fi[i]); k--;
	}
	int rest = k / 4;
	
	for (size_t i = 0; i < rest; i++)
	{
		new_p[k] = mutation_around(p[i], x0, x1);			//сделаем мутации в окрестностях точки
		k--;
	}

	for (size_t i = 0; i <= k; i++) {
		new_p[i] = mutation(p[i], x0, x1, fi[i]);//, x1, iter, max_iter);
	}
	return new_p;
}

//сортировка от наименьших значений к наибольшим
void sort(vector<vector<double>>& p, vector<double>& fi) {
	int size = p.size();
	for (int i = 0; i < size; i++)
		for (int j = i + 1; j < size; j++)
			if (fi[j] < fi[i]) {
				std::swap(p[i], p[j]);
				std::swap(fi[i], fi[j]);
			}
}

std::vector<double> genetic_alg(std::function<double(vector<double>)> f, const double eps, const std::vector<double> x0, const std::vector<double> x1)
{
	int iter = 0, n = x0.size(), size = 50, max_iter = 500;
	vector<vector<double>> population(size);
	vector<double> fi(size);					//значения функций
	for (size_t i = 0; i < size; i++)
		population[i].resize(n);

	srand(time(NULL));
	for (int i = 0; i < size; i++)     // Формирование начальной популяции
	{
		for (size_t j = 0; j < n; j++)
		{
			population[i][j] = random_state0(x0[j], x1[j]);
		}
		fi[i] = f(population[i]);
	}
	sort(population, fi);

	while(f(population[0]) > eps && iter < max_iter) {
		iter++;
		std::cout << iter << " ";
		population = crossover(population, eps, x0, x1, iter, max_iter, fi);
		for (int i = 0; i < size; i++) {
			if (population[i][0] < 0.00000001 || population[i][2] < 0.0000001)  
				fi[i] = 1000000;
			else
				fi[i] = f(population[i]);
		}
		sort(population, fi);
	}
	std::cout << "\nКоличество итераций: " << iter << std::endl;
	return population[0];	
}

//==================================================НОВЫЙ МЕТОД НЕЛДЕРА-МИДА==============================================================

struct Point {
	vector<double> coords;
	double value;
};

bool compare(const Point& a, const Point& b) {
	return a.value < b.value;
}

Point centroid(const vector<Point>& simplex) {
	size_t dim = simplex[0].coords.size();
	vector<double> center(dim, 0.0);
	for (size_t i = 0; i < simplex.size() - 1; ++i) {
		for (size_t j = 0; j < dim; ++j) {
			center[j] += simplex[i].coords[j];
		}
	}
	for (size_t j = 0; j < dim; ++j) {
		center[j] /= (simplex.size() - 1);
	}
	return { center, 0.0 };
}

Point reflect(const Point& centroid, const Point& worst, const double alpha) {
	size_t dim = centroid.coords.size();
	vector<double> reflected(dim);
	for (size_t i = 0; i < dim; ++i) {
		reflected[i] = centroid.coords[i] + alpha * (centroid.coords[i] - worst.coords[i]);
	}
	return { reflected, 0.0 };
}

Point expand(const Point& centroid, const Point& reflected, const double gamma) {
	size_t dim = centroid.coords.size();
	vector<double> expanded(dim);
	for (size_t i = 0; i < dim; ++i) {
		expanded[i] = centroid.coords[i] + gamma * (reflected.coords[i] - centroid.coords[i]);
	}
	return { expanded, 0.0 };
}

Point contract(const Point& centroid, const Point& worst, const double rho) {
	size_t dim = centroid.coords.size();
	vector<double> contracted(dim);
	for (size_t i = 0; i < dim; ++i) {
		contracted[i] = centroid.coords[i] + rho * (worst.coords[i] - centroid.coords[i]);
	}
	return { contracted, 0.0 };
}

void shrink(vector<Point>& simplex, const double sigma) {
	for (size_t i = 1; i < simplex.size(); ++i) {
		for (size_t j = 0; j < simplex[i].coords.size(); ++j) {
			simplex[i].coords[j] = simplex[0].coords[j] + sigma * (simplex[i].coords[j] - simplex[0].coords[j]);
		}
	}
}

vector<Point> generate_simplex(const vector<double>& x0, double l) {
	int dim = x0.size();
	double n = static_cast<double>(dim);
	double s = l*(sqrt(n + 1.0) - 1.0 + n);
	s /= sqrt(2.0) * n;
	double r = l*(sqrt(n + 1.0) - 1.0);
	r /= sqrt(2.0) * n;

	vector<Point> simplex(dim + 1);
	for (int i = 0; i <= dim; ++i) {
		vector<double> coords = x0;
		for (int j = 0; j < dim; ++j) {
			if (j == i - 1) {
				coords[j] += r;
			}
			else if (i > 0) {
				coords[j] += s;
			}
		}
		simplex[i] = { coords, 0.0 };
	}
	return simplex;
}

vector<double> nelder_mead(function<double(vector<double>)> f, const vector<double> x0, const double l, const double eps,
	const double alpha, const double gamma, const double rho, const double sigma, const int max_iter) {

	auto simplex = generate_simplex(x0, l);

	for (auto& p : simplex) p.value = f(p.coords);
	int iter = 0;
	double atc = 1;
	while (iter < max_iter && atc > eps) {
		sort(simplex.begin(), simplex.end(), compare);
		atc = 0;
		for (size_t i = 1; i < simplex.size(); i++)
		{
			atc += pow(simplex[i].value - simplex[0].value, 2);
		}
		atc = sqrt(atc / simplex.size());
		Point c = centroid(simplex);
		Point r = reflect(c, simplex.back(), alpha);
		r.value = f(r.coords);

		if (r.value < simplex.front().value) {
			Point e = expand(c, r, gamma);
			e.value = f(e.coords);
			simplex.back() = (e.value < r.value) ? e : r;
		}
		else if (r.value < simplex[simplex.size() - 2].value) {
			simplex.back() = r;
		}
		else {
			Point con = contract(c, simplex.back(), rho);
			con.value = f(con.coords);
			if (con.value < simplex.back().value) {
				simplex.back() = con;
			}
			else {
				shrink(simplex, sigma);
				for (auto& p : simplex) p.value = f(p.coords);
			}
		}
		iter++;
	}
	sort(simplex.begin(), simplex.end(), compare);

	return simplex.front().coords;
}
