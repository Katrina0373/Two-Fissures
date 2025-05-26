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
/// ������������ �������� ��� �����������
/// </summary>
/// <param name="f">- �������</param>
/// <param name="eps">- ��������</param>
/// <param name="x0, x1">- ������ � ������� ������� ��������� �������</param>
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

