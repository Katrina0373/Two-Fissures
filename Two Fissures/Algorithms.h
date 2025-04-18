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

//M���� ������ ��� ���������� ���������� ������� �����������
std::pair<double, double> local_segment_Svenn(double x0, double delta, double (*f)(double));

//=============================������ ����������� �������� �������===============================================

/// <summary>
/// ����� ������������ ���������� ������
/// </summary>
/// <param name="a, b">- ������ � ����� �������</param>
/// <param name="eps">- ��������</param>
/// <returns>������� ����������� �������� ��������</returns>
double passive_serach(double a, double b, double (*f)(double), double eps);

/// <summary>
/// ����� ���������
/// </summary>
/// <param name="a, b">- ������ � ����� �������</param>
/// <param name="eps">- ��������</param>
/// <param name="delta">- ��� ����, ������ ���� ������ eps</param>
/// <returns>������� ����������� �������� ��������</returns>
double dichotomy(double a, double b, double (*f)(double), double eps, double delta);

/// <summary>
/// ����� ����������� �������
/// </summary>
/// <param name="a, b">- ������ � ����� �������</param>
/// <param name="eps">- ��������</param>
/// <returns>����������� ������� x*</returns>
double divide_by_half(double a, double b, double (*f)(double), double eps);

/// <summary>
/// ����� ��������� ������, ��������� �� ������������������ ���������
/// </summary>
/// <param name="a, b">- ������ � ����� �������</param>
/// <param name="f">- �������</param>
/// <param name="eps">- �����������</param>
/// <param name="delta">- ����� ��������</param>
/// <returns>������� �����������</returns>
double Fibonachi_method(double a, double b, double (*f)(double), double eps, double delta);

/// <summary>
/// ����� �������� �������
/// </summary>
/// <param name="a, b">- ������ � ����� �������</param>
/// <param name="f">- �������</param>
/// <param name="eps">- �����������</param>
/// <returns>����������� ������� x</returns>
double gold_sector(double a, double b, double (*f)(double), double eps);

/// <summary>
/// ����� ������������ ����������� (����� �������)
/// </summary>
/// <param name="a, b">- ������ � ����� �������</param>
/// <param name="f">- �������</param>
/// <param name="eps1, eps2">- ����� �������� �����������</param>
/// <param name="delta">- ��� ����, ���� ������� �� �������� �� ����� �������� ��� ��������</param>
/// <returns>����������� ������� x</returns>
double quadr_interpolation(double a, double b, double (*f)(double), double eps, double delta);

//=============================������ ����������� n-��� �������===============================================

//����� ������� ������

/// <summary>
/// ����� ������������� ��� ����� �������-����. ������������ �� ���������� ���������
/// </summary>
/// <param name="f">- ������� �����������</param>
/// <param name="eps">- ��������</param>
/// <param name="alfa">- ����������� ��������� > 0</param>
/// <param name="beta">- ����������� ���������� > 1</param>
/// <param name="gamma">- ����������� ������ (0, 1)</param>
/// <param name="x0">- ��������� �����</param>
/// <param name="l">- ����� ����� ���������</param>
/// <returns></returns>
std::vector<double> polyhedrom_method(double(*f)(std::vector<double>), double eps,
	std::vector<double> x0, double alfa, double beta, double gamma, double l, int max_iter = 500);

/// <summary>
/// ����������� ����� ������������� ������
/// </summary>
/// <param name="eps">- �������� ����</param>
/// <param name="eps1">- �������� �������� �����</param>
/// <param name="delta">- ���������� ������ 0.01</param>
/// <param name="delta2">- ��� ��� ������ ����</param>
/// <returns></returns>
std::vector<double> gradient_method_of_steepest_descent(double(*f)(std::vector<double>),
	double eps, std::vector<double> x0, double eps1, double delta, double delta2);

/// <summary>
/// ����������� ����� � ���������� ����. �� ����� ������������ � ����������
/// </summary>
/// <param name="eps">- �������� ���� (������������� 0.1)</param>
/// <param name="eps1">- �������� �������� �������</param>
/// <param name="lambda">- ��� (1.0)</param>
/// <param name="beta">- ��������� ���� (0.95)</param>
/// <returns></returns>
std::vector<double> gradient_method_with_step_splitting(double(*f)(std::vector<double>),
	double eps, std::vector<double> x0, double eps1, double lambda, double beta);

/// <summary>
/// ����� ����������� ����������
/// </summary>
/// <returns></returns>
std::vector<double> conjugate_gradient(double(*f)(std::vector<double>), double eps, std::vector<double> x0);

/// <summary>
/// ����� ������� ��� �����������.
/// ���������� ����������� ������� �������
/// </summary>
/// <param name="f">- ������� �����������</param>
/// <param name="eps">- ��������</param>
/// <param name="x0">- ��������� �����</param>
/// <returns></returns>
std::vector<double> newton_method(double (*f)(std::vector<double>), double eps, std::vector<double> x0);

/// <summary>
/// ����� �������-�������
/// </summary>
/// <param name="f">- ������� �����������</param>
/// <param name="eps">- ��������</param>
/// <param name="x0">- ��������� �����</param>
/// <returns></returns>
std::vector<double> newton_raphson_method(double (*f)(std::vector<double>), double eps, std::vector<double> x0);

/// <summary>
/// ������������ �������� ��� �����������
/// </summary>
/// <param name="f">- �������</param>
/// <param name="eps">- ��������</param>
/// <param name="x0, x1">- ������� � ������ ������� ��������� �������</param>
/// <returns></returns>
std::vector<double> genetic_alg(const std::function<double(std::vector<double>)> f, const double eps, const std::vector<double> x0, const std::vector<double> x1);

/// <summary>
/// ����� �������-����
/// </summary>
/// <param name="f">- ������� �����������</param>
/// <param name="x0">- ��������� �����</param>
/// <param name="l">- ����� ���������</param>
/// <param name="alpha">- ���������</param>
/// <param name="gamma">- ����������</param>
/// <param name="rho"></param>
/// <param name="sigma">- ������</param>
/// <param name="max_iter">- ������������ ���������� ��������</param>
/// <returns>����� �������� �������</returns>
vector<double> nelder_mead(const function<double(vector<double>)> f, const vector<double> x0,const double l,const double eps,
	const double alpha = 1.0,const double gamma = 2.0, const double rho = 0.5, const double sigma = 0.5, const int max_iter = 200);

#endif // !ALGTHS_H

