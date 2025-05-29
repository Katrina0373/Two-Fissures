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


class Genetic_Alg {
private:
	const size_t pop_size = 50; //������ ���������
	size_t n; //���������� ����� � ����� �����
	function<double(vector<double>)> f; //������� ������� ������ ����������
	vector<double> x0, x1; //������ � ������� ������� ��������� ��������
	const size_t max_iter = 500; //������������ ���������� ��������
	double eps;
	//��� �������������
	const double d = 0.2, Mm = 0.9;
	
	// ������� �� ���� �������
	vector<double> mutation(const vector<double> point, const double fi);  

	vector<double> inversion(vector<double> x, const double fi);

	void mutation_around(vector<double>& x);

	//������������������
	vector<double> mutation_around2(vector<double> x, const double fi);

	//���� �����������
	vector<double> intermediate_recombination(const vector<double>& p1, const vector<double>& p2);
	vector<double> blend_crossover(const vector<double> p1, const vector<double> p2);

	//����� ��������� ������� �������
	size_t roulette_wheel_selection(const vector<vector<double>>& p, const vector<double>& fi, double total_sum);

	//���������� �� ���������� �������� � ����������
	void sort(vector<vector<double>>& p, vector<double>& fi);

	//���������� �� �������
	double square_dist(const vector<double> p1, const vector<double> p2);

	//����������� � ������������ ����� ���������
	void crossover(vector<vector<double>>& p, const int iter, vector<double>& fi);

	//���������
	void shuffle(vector<vector<double>>& p, vector<double>& fi);
public:
	Genetic_Alg(function<double(vector<double>)> f, vector<double> lower_bound, vector<double> upper_bound, double eps);
	//������
	std::vector<double> genetic_alg();

};

#endif // !ALGTHS_H

