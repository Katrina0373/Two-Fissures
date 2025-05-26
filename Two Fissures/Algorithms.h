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
/// Генетический алгоритм для минимизации
/// </summary>
/// <param name="f">- функция</param>
/// <param name="eps">- точность</param>
/// <param name="x0, x1">- нижняя и верхняя граница возможных значний</param>
/// <returns></returns>
std::vector<double> genetic_alg(const std::function<double(std::vector<double>)> f, const double eps, const std::vector<double> x0, const std::vector<double> x1);

/// <summary>
/// Метод Нелдера-Мида
/// </summary>
/// <param name="f">- функция минимизации</param>
/// <param name="x0">- начальная точка</param>
/// <param name="l">- длина симплекса</param>
/// <param name="alpha">- отражение</param>
/// <param name="gamma">- растежение</param>
/// <param name="rho"></param>
/// <param name="sigma">- сжатие</param>
/// <param name="max_iter">- максимальное количество итераций</param>
/// <returns>точка минимума функции</returns>
vector<double> nelder_mead(const function<double(vector<double>)> f, const vector<double> x0,const double l,const double eps,
	const double alpha = 1.0,const double gamma = 2.0, const double rho = 0.5, const double sigma = 0.5, const int max_iter = 200);

#endif // !ALGTHS_H

