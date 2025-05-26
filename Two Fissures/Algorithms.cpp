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


//==================================================Генетический алгоритм====================================================


double random_state0(const double min, const double max) {
	std::random_device rd; // Источник случайных чисел
	std::mt19937 gen(rd()); // Генератор случайных чисел с инициализацией
	std::uniform_real_distribution<> dis(min, max); // Распределение на [min, max)
	return dis(gen);
}


vector<double> mutation(const vector<double> point, const vector<double> x0, const vector<double> x1, const double fi)  // мутация: генерация случайной величины
{
	const int NUM = 100000000;
	vector<double> x(point);
	for (size_t i = 0; i < x0.size(); i++)
	{
		if (fi < 0.0001 && i % 2 == 0)				//если функция уже меньше какого-то значения то возможно l у неё более менее точные но с d всегда проблема
			continue;
		x[i] = fabs((double)((rand() * NUM) % (int)((x1[i] - x0[i]) * NUM) + 1) / NUM) + x0[i];
	}
	return x;
}

vector<double> inversion(vector<double> x, const double fi) {
	static int sign = 0;
	sign++;
	sign %= 2;
	double eps1 = 0.1;
	double eps2 = 0.01;
	if (sign == 0) {
		for (size_t i = 0; i < x.size(); i++) {
			if (i % 2 == 0) {
				if (fi > 0.0001)
					x[i] -= eps2;
				continue;
			}
			x[i] -= eps1;
		}
	}
	else {
		for (size_t i = 0; i < x.size(); i++) {
			if (i % 2 == 0) {
				if (fi > 0.0001)
					x[i] += eps2;
				continue;
			}
			x[i] += eps1;
		}
	}
	return x;
}

vector<double> mutation_around(vector<double> x, const vector<double> x0, const vector<double> x1, const double fi){
	for (size_t i = 0; i < x.size(); i++)
	{
		if (i % 2 == 0 && fi < 0.0001) {
			continue;
		}
		double eps = (x1[i] - x0[i]) / 5;
		x[i] = random_state0(x[i] - eps, x[i] + eps);
	}
	return x;
}

vector<double> mutation_around2(vector<double> x, const vector<double> x0, const vector<double> x1, const double fi) {
	const double k = 1.0;
	for (size_t i = 0; i < x.size(); i++)
	{
		if (i % 2 == 0 && fi < 0.0001) {
			continue;
		} 
		double base_eps = (x1[i] - x0[i]) / 5.0;
		//double decay = exp(-k / fi);
		double rand_factor = random_state0(0.5, 2);  // Множитель случайности
		double eps = base_eps * rand_factor;

		x[i] = random_state0(x[i] - eps, x[i] + eps);
	}
	return x;
}


vector<double> intermediate_recombination(const vector<double>& p1, const vector<double>& p2) {
	vector<double> x(p1.size());
	for (size_t i = 0; i < p1.size(); i++)
	{
		double alpha = random_state0(-0.25, 1.25);
		x[i] = p1[i] + alpha * (p2[i] - p1[i]);
	}
	return x;
}

vector<double> blend_crossover(const vector<double> p1, const vector<double > p2, const vector<double > x0, const vector<double> x1) {

	double d, min_val, max_val;
	double alpha = 0.3;
	vector<double> x_new(p1.size());
	for (size_t i = 0; i < p1.size(); i++)
	{
		d = abs(p1[i] - p2[i]);
		min_val = std::max(std::min(p1[i], p2[i]) - alpha * d, x0[i]);
		max_val = std::min(std::max(p1[i], p2[i]) + alpha * d, x1[i]);
		x_new[i] = random_state0(min_val, max_val);
	}
	return x_new;
}

size_t roulette_wheel_selection(const vector<vector<double>>& p, const vector<double>& fi, double total_sum) {
	
	double pick = random_state0(0, total_sum);

	double current = 0;
	for (size_t i = 0; i < p.size(); i++)
	{
		current += (1 / fi[i]) / total_sum;
		if (current >= pick)
			return i;
	}
	return 0;
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


//Скрещивание и формирование новой популяции
void crossover(std::function<double(vector<double>)> f, vector<vector<double>>& p, const double eps, const vector<double>& x0, const vector<double>& x1,
	const int iter, const int max_iter, vector<double>& fi)	{
	int n = p.size();
	vector<vector<double>> new_p(n/2);
	vector<double> new_fi(n / 2);

	double total_sum = 0;
	for (size_t i = 0; i < n; i++)
	{
		total_sum += 1 / fi[i];
	}
	//промежуточная рекомбинация с выбором родителей 
	for (size_t i = 0; i < n/2; i++)
	{
		auto p1 = roulette_wheel_selection(p, fi, total_sum);
		auto p2 = roulette_wheel_selection(p, fi, total_sum);
		if (p1 == p2) {
			i--;
			continue;
		}
		new_p[i] = intermediate_recombination(p[p1],p[p2]);
		new_fi[i] = f(new_p[i]);
	}
	for (size_t i = 0; i < n/2; i++)
	{
		if (random_state0(0, 1) <= 0.5) {
			new_p[i] = mutation_around2(new_p[i], x0, x1, new_fi[i]);
			new_fi[i] = f(new_p[i]);
		}
	}

	sort(p, fi);

	new_p.insert(new_p.end(), p.begin(), p.begin() + n/2);
	new_fi.insert(new_fi.end(), fi.begin(), fi.begin()+n/2);
	p = new_p;
	fi = new_fi;
}



//перемешка
void shuffle(vector<vector<double>>& p, vector<double>& fi) {
	// Создаем индексы
	std::vector<size_t> indices(p.size());
	for (size_t i = 0; i < indices.size(); ++i)
		indices[i] = i;

	// Перемешиваем индексы
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(indices.begin(), indices.end(), g);

	// Создаем новые вектора по перемешанным индексам
	std::vector<vector<double>>    data_shuffled(p.size());
	std::vector<double> labels_shuffled(fi.size());

	for (size_t i = 0; i < indices.size(); ++i) {
		data_shuffled[i] = p[indices[i]];
		labels_shuffled[i] = fi[indices[i]];
	}
	p = data_shuffled;
	fi = labels_shuffled;
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

	while(f(population[0]) > eps && iter < max_iter) {
		iter++;
		std::cout << iter << " ";
		shuffle(population, fi);
		crossover(f, population, eps, x0, x1, iter, max_iter, fi);
	}
	sort(population, fi);
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
