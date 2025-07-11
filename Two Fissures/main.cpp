#include<iostream>
#include"find_fissure.h"
#include"Algorithms.h"
#include<chrono>
#include<format>
#include"solve_eqs.h"

using namespace std;

vector<cx_double> true_u;
vector<double> view_points;

double resid_func(vector<double> parameters) {
	if (parameters[0] < 1e-7 || parameters[2] < 1e-7
		|| parameters[1] < 0 || parameters[3] < 0)
		return 1000000;
	Fissures f = Fissures();
	f.set_parameters(parameters[0], parameters[1], parameters[2], parameters[3]);
	f.solve_xi();
	double res = 0;
	for (size_t i = 0; i < view_points.size(); i++)
	{
		cx_double cur_u = f.number_field(view_points[i]);
		res += pow(cur_u.real() - true_u[i].real(), 2) * 100 + pow(cur_u.imag() - true_u[i].imag(), 2) * 100;
	}
	return res;
}
//
//void four_fields() {
//
//	Fissures f = Fissures();
//	double l1, l2, d1, d2;
//	l1 = 0.1;
//	d1 = 0.5;
//	l2 = 0.1;
//	d2 = 1.0;
//	f.set_parameters(l1, d1, l2, d2);
//
//	f.eval_static_vecs();
//	f.solve_xi();
//
//	auto x = arma::linspace(0.0, 10.0, 2000);
//	std::ofstream file("D:/VS Projects/Two Fissures/results/three_fields.csv");
//	if (!file.is_open()) {
//		std::cerr << "�� ������� ������� ���� ��� ������!" << std::endl;
//		return;
//	}
//	file << f.l1 << endl;
//	for (size_t i = 0; i < x.size(); i++)
//	{
//		file << f.number_field(x(i)).real() << ';';
//	}
//	file << endl;
//
//	f.l1 = 0.05;
//	f.solve_xi();
//	file << f.l1 << endl;
//	for (size_t i = 0; i < x.size(); i++)
//		file << f.number_field(x(i)).real() << ';';
//	file << endl;
//
//	f.l1 = 0.1;
//	f.solve_xi();
//	file << f.l1 << endl;
//	for (size_t i = 0; i < x.size(); i++)
//		file << f.number_field(x(i)).real() << ';';
//	file << endl;
//
//	file.close();
//
//}

void find_coordinates() {

	double l1, l2, d1, d2;
	cout << "������� ������� ����������:\nl1 = "; cin >> l1;
	cout << "d1 = "; cin >> d1;
	cout << "l2 = "; cin >> l2;
	cout << "d2 = "; cin >> d2;

	cout << "������ �������� ��� ���������(����� ������ � ��� �� �������): ";
	vector<double> x0;
	for (size_t i = 0; i < 4; i++)
	{
		double a;
		cin >> a;
		x0.push_back(a);
	}
	cout << "������� �������� ��� ���������: ";
	vector<double> x1;
	for (size_t i = 0; i < 4; i++)
	{
		double a;
		cin >> a;
		x1.push_back(a);
	}
	double eps;
	cout << "eps = "; cin >> eps;

	Fissures f;
	f.set_parameters(l1, d1, l2, d2);
	f.eval_static_vecs();
	f.solve_xi();
	cout << "������� �������� x1, x2 ��� ���� ��������\nx1 = ";
	double xx1, xx2;
	cin >> xx1;
	cout << "x2 = "; cin >> xx2;

	view_points = { xx1, xx2 };
	true_u.resize(view_points.size());
	true_u = { f.number_field(xx1), f.number_field(xx2) };

	auto start_time = std::chrono::high_resolution_clock::now();
	auto gen = Genetic_Alg(resid_func, x0, x1, eps);
	auto res = gen.genetic_alg();
	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end_time - start_time;
	std::cout << "\n����� ������ ��������� � ���: " << duration.count() << std::endl << endl;
	
	printf("��������� ���������: %.5f; %.5f; %.5f; %.5f\n", res[0], res[1], res[2], res[3]);
	cout << "�������� �������: " << resid_func(res) << endl;

	system("pause");
}

void Minimize_with_Nelder_Mid(const double l1, const double d1, const double l2, const double d2, const double k,
	const vector<double> point, const double eps, const double l) {
	
	Fissures f = Fissures(l1, d1, l2, d2, 20, k);
	f.eval_static_vecs();
	f.solve_xi();
	true_u.resize(view_points.size());
	for (size_t i = 0; i < view_points.size(); i++)
		true_u[i] = f.number_field(view_points[i]);
	cout << resid_func(point) << endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	//cout << resid_func({ l1, d1, l2, d2 }) << endl;
	auto res = nelder_mead(resid_func, point, l, eps);
	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end_time - start_time;
	std::cout << "����� ������ ���������: " << duration.count() << std::endl;
	printf("��������� ���������: %.5f; %.5f; %.5f; %.5f\n",
		res[0], res[1], res[2], res[3]);
	cout << endl;
	cout << "�������� �������: " << resid_func(res) << endl;
}


void task4(const double l1, const double d1, const double l2, const double d2, const double k,
	const vector<double> floor, const vector<double> roof, const double eps1,
	const double eps2, const double l,
	string field_file_name, string report_file_name
) {

	Fissures f = Fissures(l1, d1, l2, d2, 20, k);
	f.eval_static_vecs();
	f.solve_xi();
	f.write_bundle_fun_to_csv(std::format("../results/������/������� ����������/k{}l1{}d1{}l2{}d2{}.csv", k, l1, d1, l2, d2));
	f.write_field_to_csv(field_file_name);
	
	//cout << "������� ����� ����������: ";
	/*view_points.clear();
	while (true) {
		double a;
		cin >> a;
		if (a == 0)
			break;
		view_points.push_back(a);
	}*/

	true_u.resize(view_points.size());
	for (size_t i = 0; i < view_points.size(); i++)
		true_u[i] = f.number_field(view_points[i]);

	//������� � ������� �������������
	cout << "������ ������ ������������� ���������" << endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	auto gen = Genetic_Alg(resid_func, floor, roof, eps1);
	auto res1 = gen.genetic_alg();
	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end_time - start_time;
	std::cout << "����� ������ ���������: " << duration.count() << std::endl;
	printf("��������� ��������� � ������� ������������� ���������: %.5f; %.5f; %.5f; %.5f\n",
		res1[0], res1[1], res1[2], res1[3]);
	cout << endl;
	cout << "�������� �������: " << resid_func(res1) << endl;

	//�������� � ������� �������-����
	cout << "������ ������ ��������� �������-����" << endl;
	start_time = std::chrono::high_resolution_clock::now();
	auto res2 = nelder_mead(resid_func, res1, l, eps2);
	end_time = std::chrono::high_resolution_clock::now();
	duration = end_time - start_time;
	std::cout << "����� ������ ���������: " << duration.count() << std::endl;
	printf("��������� ���������: %.5f; %.5f; %.5f; %.5f\n",
		res2[0], res2[1], res2[2], res2[3]);
	cout << endl;
	cout << "�������� �������: " << resid_func(res2) << endl;

	std::ofstream file(report_file_name, std::ios::app);
	if (!file.is_open()) {
		std::cerr << "�� ������� ������� ���� ��� ������!" << std::endl;
		return;
	}
	file << endl<< "�������� ���������: " << endl;
	file << "l1 = " << l1 << endl;
	file << "d1 = " << d1 << endl;
	file << "l2 = " << l2 << endl;
	file << "d2 = " << d2 << endl;
	file << "��������� � ������� �������������: " << endl;
	file << "l1 = " << res1[0] << endl;
	file << "d1 = " << res1[1] << endl;
	file << "l2 = " << res1[2] << endl;
	file << "d2 = " << res1[3] << endl;
	file << "��������� � ������� �������-����: " << endl;
	file << "l1 = " << res2[0] << endl;
	file << "d1 = " << res2[1] << endl;
	file << "l2 = " << res2[2] << endl;
	file << "d2 = " << res2[3] << endl;
	file << "���. ����������: " << endl;
	file << "����� ����������:";
	for (size_t i = 0; i < view_points.size(); i++)
		file << " x" << i + 1 << " = " << view_points[i];
	file << endl;
	file << "������ ����� ����������: " << floor[0] << ", " << floor[1] << ", " << floor[2] << ", " << floor[3] << endl;
	file << "������� ����� ����������: " << roof[0] << ", " << roof[1] << ", " << roof[2] << ", " << roof[3] << endl;
	file << "eps ������������� = " << eps1 << endl;
	file << "eps �.-�. = " << eps2 << endl;
	file << "����� �������� = " << l << endl;
	file << "������� k = " << k << endl;
	file << "�������� ������� �������������� �� ������������ = " << resid_func(res1) << endl;
	file << "�������� ������� �������������� �� �.-�. = " << resid_func(res2) << endl << endl;

	file.close();
}



int main() {
	setlocale(LC_ALL, "Russian");
	
	double l1, l2, d1, d2;
	l1 = 0.05;
	d1 = 3.5;
	l2 = 0.1;
	d2 = 1.0;
	vector<double> roof = { 0.2, 10.0, 0.2, 10.0 };
	vector<double> floor = { 0.000001, 0.0, 0.000001, 0.0 };
	double k = 5, l = 0.1;
	double eps1 = 1e-6, eps2 = 1e-10;
	view_points = { 2, 2.4, 2.7, 3, 3.4, 3.7, 4 };

	Fissures f = Fissures(l1, d1, l2, d2, 20, k);

	task4(l1, d1, l2, d2, k, floor, roof, eps1, eps2, l,
		std::format("../results/������/���� ��������/k{}l1{}d1{}l2{}d2{}.csv", k, l1, d1, l2, d2),
		"../results/report4.0.txt");
	 
	

	//Minimize_with_Nelder_Mid(f, { 0.0977921, 0.525886, 0.0223358, 1.99891 }, eps2, 0.1);
	return 0;
}