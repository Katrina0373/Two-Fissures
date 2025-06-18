#include"Algorithms.h"
#include<iostream>
#include<random>
using std::vector;

//==================================================Генетический алгоритм====================================================


double random_state0(const double min, const double max) {
	std::random_device rd; // Источник случайных чисел
	std::mt19937 gen(rd()); // Генератор случайных чисел с инициализацией
	std::uniform_real_distribution<> dis(min, max); // Распределение на [min, max)
	return dis(gen);
}

Genetic_Alg::Genetic_Alg(function<double(vector<double>)> f, vector<double> lower_bound, vector<double> upper_bound, double eps) {
	this->f = f;
	n = lower_bound.size();
	x0 = lower_bound;
	x1 = upper_bound;
	this->eps = eps;
}

vector<double> Genetic_Alg::mutation(const vector<double> point, const double fi)  // мутация: генерация случайной величины
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

void Genetic_Alg::mutation_around(vector<double>& x) {
	for (size_t i = 0; i < x.size(); i++)
	{
		double eps = (x1[i] - x0[i]) / 5;
		x[i] = random_state0(x[i] - eps, x[i] + eps);
	}
}

vector<double> Genetic_Alg::mutation_around2(vector<double> x, const double fi) {
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


vector<double> Genetic_Alg::intermediate_recombination(const vector<double>& p1, const vector<double>& p2) {
	vector<double> x(p1.size());
	for (size_t i = 0; i < p1.size(); i++)
	{
		double alpha = random_state0(-0.25, 1.25);
		x[i] = p1[i] + alpha * (p2[i] - p1[i]);
	}
	return x;
}

vector<double> Genetic_Alg::blend_crossover(const vector<double> p1, const vector<double > p2) {

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

size_t Genetic_Alg::roulette_wheel_selection(const vector<vector<double>>& p, const vector<double>& fi, double total_sum) {

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
void Genetic_Alg::sort(vector<vector<double>>& p, vector<double>& fi) {
	for (size_t i = 0; i < pop_size; i++)
		for (size_t j = i + 1; j < pop_size; j++)
			if (fi[j] < fi[i]) {
				std::swap(p[i], p[j]);
				std::swap(fi[i], fi[j]);
			}
}

double Genetic_Alg::square_dist(const vector<double> p1, const vector<double> p2) {
	double dist = 0;
	for (size_t i = 0; i < n; i++)
	{
		dist += pow(abs((p1[i] - p2[i]) / (x0[i] - x1[i])), d);
	}
	return dist / n;
}

//Скрещивание и формирование новой популяции
void Genetic_Alg::crossover(vector<vector<double>>& p,	const int iter, vector<double>& fi) {
	size_t new_pop_size = 0.75 * pop_size;
	vector<vector<double>> new_p(new_pop_size);
	vector<double> new_fi(new_pop_size);

	double total_sum = 0;
	for (size_t i = 0; i < pop_size; i++)
	{
		total_sum += 1 / fi[i];
	}
	for (size_t i = 0; i < new_pop_size; i++)
	{

		//промежуточная рекомбинация с выбором родителей 
		auto p1_ind = roulette_wheel_selection(p, fi, total_sum);
		auto p2_ind = roulette_wheel_selection(p, fi, total_sum);
		if (p1_ind == p2_ind) {
			p2_ind = rand() % pop_size;
		}
		auto child = intermediate_recombination(p[p1_ind], p[p2_ind]);
		//вычисляем вероятность мутации
		double dist = square_dist(p[p1_ind], p[p2_ind]);
		if (random_state0(0, 1) <= (1 - dist) * Mm) {
			mutation_around(child);
		}
		double child_fi = f(child);

		new_p[i] = child;
		new_fi[i] = child_fi;
	}

	sort(p, fi);
	new_p.insert(new_p.begin(), p.begin(), p.begin() + pop_size - new_pop_size);
	new_fi.insert(new_fi.begin(), fi.begin(), fi.begin() + pop_size - new_pop_size);

	p = new_p;
	fi = new_fi;
}



//перемешка
void Genetic_Alg::shuffle(vector<vector<double>>& p, vector<double>& fi) {
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

std::vector<double> Genetic_Alg::genetic_alg()
{
	int iter = 0;
	vector<vector<double>> population(pop_size);
	vector<double> fi(pop_size);					//значения функций
	for (size_t i = 0; i < pop_size; i++)
		population[i].resize(n);

	srand(time(NULL));
	for (int i = 0; i < pop_size; i++)     // Формирование начальной популяции
	{
		for (size_t j = 0; j < n; j++)
		{
			population[i][j] = random_state0(x0[j], x1[j]);
		}
		fi[i] = f(population[i]);
	}

	while (*(std::min_element(fi.begin(), fi.end())) > eps && iter < max_iter) {
		iter++;
		std::cout << iter << " ";
		shuffle(population, fi);
		crossover(population, iter, fi);
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
