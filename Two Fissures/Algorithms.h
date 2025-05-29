#ifndef ALGTHS_H
#define ALGTHS_H

#include<math.h>
#include<tuple>
#include<float.h>
#include<vector>
#include<iostream>
#include<functional>
#include<numeric>
#include<ctime>
#include<functional>

using std::function;
using std::vector;

/// <summary>
/// ћетод Ќелдера-ћида
/// </summary>
/// <param name="f">- функци€ минимизации</param>
/// <param name="x0">- начальна€ точка</param>
/// <param name="l">- длина симплекса</param>
/// <param name="alpha">- отражение</param>
/// <param name="gamma">- растежение</param>
/// <param name="rho"></param>
/// <param name="sigma">- сжатие</param>
/// <param name="max_iter">- максимальное количество итераций</param>
/// <returns>точка минимума функции</returns>
vector<double> nelder_mead(const function<double(vector<double>)> f, const vector<double> x0,const double l,const double eps,
	const double alpha = 1.0,const double gamma = 2.0, const double rho = 0.5, const double sigma = 0.5, const int max_iter = 200);


class Genetic_Alg {
private:
	const size_t pop_size = 50; //размер попул€ции
	size_t n; //количество генов у одной особи
	function<double(vector<double>)> f; //текуща€ функци€ поиска экстремума
	vector<double> x0, x1; //нижн€€ и верхн€€ граница изменени€ значений
	const size_t max_iter = 500; //максимальное количество итераций
	double eps;
	//дл€ самоадаптации
	const double d = 0.2, Mm = 0.9;
	
	// мутаци€ по всей границе
	vector<double> mutation(const vector<double> point, const double fi);  

	vector<double> inversion(vector<double> x, const double fi);

	void mutation_around(vector<double>& x);

	//самоадаптирующийс€
	vector<double> mutation_around2(vector<double> x, const double fi);

	//виды скрещиваний
	vector<double> intermediate_recombination(const vector<double>& p1, const vector<double>& p2);
	vector<double> blend_crossover(const vector<double> p1, const vector<double> p2);

	//выбор родителем методом рулетки
	size_t roulette_wheel_selection(const vector<vector<double>>& p, const vector<double>& fi, double total_sum);

	//сортировка от наименьших значений к наибольшим
	void sort(vector<vector<double>>& p, vector<double>& fi);

	//рассто€ние по площади
	double square_dist(const vector<double> p1, const vector<double> p2);

	//—крещивание и формирование новой попул€ции
	void crossover(vector<vector<double>>& p, const int iter, vector<double>& fi);

	//перемешка
	void shuffle(vector<vector<double>>& p, vector<double>& fi);
public:
	Genetic_Alg(function<double(vector<double>)> f, vector<double> lower_bound, vector<double> upper_bound, double eps);
	//запуск
	std::vector<double> genetic_alg();

};

#endif // !ALGTHS_H

