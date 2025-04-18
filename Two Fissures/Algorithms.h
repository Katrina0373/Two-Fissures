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

//Mетод Свенна для построения начального отрезка локализации
std::pair<double, double> local_segment_Svenn(double x0, double delta, double (*f)(double));

//=============================Методы минимизации нулевого порядка===============================================

/// <summary>
/// Метод равномерного пассивного поиска
/// </summary>
/// <param name="a, b">- начало и конец отрезка</param>
/// <param name="eps">- точность</param>
/// <returns>отрезок локализации заданной точности</returns>
double passive_serach(double a, double b, double (*f)(double), double eps);

/// <summary>
/// Метод дихотомии
/// </summary>
/// <param name="a, b">- начало и конец отрезка</param>
/// <param name="eps">- точность</param>
/// <param name="delta">- шаг проб, должно быть меньше eps</param>
/// <returns>отрезок локализации заданной точности</returns>
double dichotomy(double a, double b, double (*f)(double), double eps, double delta);

/// <summary>
/// метод половинного деления
/// </summary>
/// <param name="a, b">- начало и конец отрезка</param>
/// <param name="eps">- точность</param>
/// <returns>оптимальное решение x*</returns>
double divide_by_half(double a, double b, double (*f)(double), double eps);

/// <summary>
/// Метод активного поиска, основаный на последовательности Фибоначчи
/// </summary>
/// <param name="a, b">- начало и конец отрезка</param>
/// <param name="f">- функция</param>
/// <param name="eps">- погрешность</param>
/// <param name="delta">- малая величина</param>
/// <returns>отрезок локализации</returns>
double Fibonachi_method(double a, double b, double (*f)(double), double eps, double delta);

/// <summary>
/// Метод золотого сечения
/// </summary>
/// <param name="a, b">- начало и конец отрезка</param>
/// <param name="f">- функция</param>
/// <param name="eps">- погрешность</param>
/// <returns>оптимальное решение x</returns>
double gold_sector(double a, double b, double (*f)(double), double eps);

/// <summary>
/// Метод квадратичной интеполяции (метод Пауэлла)
/// </summary>
/// <param name="a, b">- начало и конец отрезка</param>
/// <param name="f">- функция</param>
/// <param name="eps1, eps2">- малые велечины погрешности</param>
/// <param name="delta">- шаг проб, если функция не работает то стоит поменять это значение</param>
/// <returns>оптимальное решение x</returns>
double quadr_interpolation(double a, double b, double (*f)(double), double eps, double delta);

//=============================Методы минимизации n-ого порядка===============================================

//Метод прямого поиска

/// <summary>
/// Метод многогранника или метод Недлера-Мида. Основывается на построении симплекса
/// </summary>
/// <param name="f">- функция минимизации</param>
/// <param name="eps">- точность</param>
/// <param name="alfa">- коэффициент отражения > 0</param>
/// <param name="beta">- коэффициент растяжения > 1</param>
/// <param name="gamma">- коэффициент сжатия (0, 1)</param>
/// <param name="x0">- начальная точка</param>
/// <param name="l">- длина ребра симплекса</param>
/// <returns></returns>
std::vector<double> polyhedrom_method(double(*f)(std::vector<double>), double eps,
	std::vector<double> x0, double alfa, double beta, double gamma, double l, int max_iter = 500);

/// <summary>
/// Градиентный метод наискорейшего спуска
/// </summary>
/// <param name="eps">- точность шага</param>
/// <param name="eps1">- точность значений точки</param>
/// <param name="delta">- желательно меньше 0.01</param>
/// <param name="delta2">- шаг для выбора шага</param>
/// <returns></returns>
std::vector<double> gradient_method_of_steepest_descent(double(*f)(std::vector<double>),
	double eps, std::vector<double> x0, double eps1, double delta, double delta2);

/// <summary>
/// Градиентный метод с дроблением шага. Не особо чувствителен к параметрам
/// </summary>
/// <param name="eps">- точность шага (рекомендуется 0.1)</param>
/// <param name="eps1">- точность значений функции</param>
/// <param name="lambda">- шаг (1.0)</param>
/// <param name="beta">- дробление шага (0.95)</param>
/// <returns></returns>
std::vector<double> gradient_method_with_step_splitting(double(*f)(std::vector<double>),
	double eps, std::vector<double> x0, double eps1, double lambda, double beta);

/// <summary>
/// Метод сопряженных градиентов
/// </summary>
/// <returns></returns>
std::vector<double> conjugate_gradient(double(*f)(std::vector<double>), double eps, std::vector<double> x0);

/// <summary>
/// Метод Ньютона для минимизации.
/// Использует производную второго порядка
/// </summary>
/// <param name="f">- функция минимизации</param>
/// <param name="eps">- точность</param>
/// <param name="x0">- начальная точка</param>
/// <returns></returns>
std::vector<double> newton_method(double (*f)(std::vector<double>), double eps, std::vector<double> x0);

/// <summary>
/// Метод Ньютона-Рафсона
/// </summary>
/// <param name="f">- функция минимизации</param>
/// <param name="eps">- точность</param>
/// <param name="x0">- начальная точка</param>
/// <returns></returns>
std::vector<double> newton_raphson_method(double (*f)(std::vector<double>), double eps, std::vector<double> x0);

/// <summary>
/// Генетический алгоритм для минимизации
/// </summary>
/// <param name="f">- функция</param>
/// <param name="eps">- точность</param>
/// <param name="x0, x1">- верхняя и нижняя граница возможных значний</param>
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

