#include<iostream>
#include"find_fissure.h"
#include"Algorithms.h"
#include<chrono>
#include<format>

using namespace std;

vector<cx_double> true_u;
vector<double> view_points;

void write_to_file(cx_vec x, double l1, double d1, double l2, double d2, double k, string file_name) {
	std::ofstream file(file_name);
	if (!file.is_open()) {
		std::cerr << "Не удалось открыть файл для записи!" << std::endl;
		return;
	}

	file << l1 << endl;
	file << d1 << endl;
	file << l2 << endl;
	file << d2 << endl;
	file << k << endl;

	for (size_t i = 0; i < x.size(); ++i) {
		file << x[i];
		if (i != x.size() - 1) {
			file << ';';
		}
	}
	file.close();

}

double resid_func(vector<double> coordinates) {
	Fissures f = Fissures();
	f.set_coordinates(coordinates[0], coordinates[1], coordinates[2], coordinates[3]);
	f.solve_xi();
	double res = 0;
	for (size_t i = 0; i < view_points.size(); i++)
	{
		auto cur_u = f.number_field(view_points[i]);
		res += pow(cur_u.real() - true_u[i].real(), 2) * 100 + pow(cur_u.imag() - true_u[i].imag(), 2) * 100;
	}

	return res;
}

void field_grahp(Fissures f, string file_name) {
	//график для поля
	f.eval_static_vecs();
	f.solve_xi();
	auto x = arma::linspace(0.0, 10.0, 1000);
	cx_vec u_vals(x.size());
	for (size_t i = 0; i < x.size(); i++)
	{
		u_vals(i) = f.number_field(x(i));
	}
	write_to_file(u_vals, f.l1, f.d1, f.l2, f.d2, f.k, file_name);
}

void four_fields() {

	Fissures f = Fissures();
	double l1, l2, d1, d2;
	l1 = 0.1;
	d1 = 0.5;
	l2 = 0.1;
	d2 = 1.0;
	f.set_coordinates(l1, d1, l2, d2);

	f.eval_static_vecs();
	f.solve_xi();

	auto x = arma::linspace(0.0, 10.0, 2000);
	std::ofstream file("D:/VS Projects/Two Fissures/results/three_fields.csv");
	if (!file.is_open()) {
		std::cerr << "Не удалось открыть файл для записи!" << std::endl;
		return;
	}
	file << f.l1 << endl;
	for (size_t i = 0; i < x.size(); i++)
	{
		file << f.number_field(x(i)).real() << ';';
	}
	file << endl;

	f.l1 = 0.05;
	f.solve_xi();
	file << f.l1 << endl;
	for (size_t i = 0; i < x.size(); i++)
		file << f.number_field(x(i)).real() << ';';
	file << endl;

	f.l1 = 0.1;
	f.solve_xi();
	file << f.l1 << endl;
	for (size_t i = 0; i < x.size(); i++)
		file << f.number_field(x(i)).real() << ';';
	file << endl;

	/*f.l1 = 0.5;
	f.solve_xi();
	f.fill_u_sigma();
	file << f.l1 << endl;
	for (size_t i = 0; i < x.size(); i++)
		file << f.number_field(x(i)).real() << ';';*/
	file.close();

}

void find_coordinates() {

	double l1, l2, d1, d2;
	cout << "Задайте искомые координаты:\nl1 = "; cin >> l1;
	cout << "d1 = "; cin >> d1;
	cout << "l2 = "; cin >> l2;
	cout << "d2 = "; cin >> d2;

	cout << "Нижний диапазон для координат(через пробел в том же порядке): ";
	vector<double> x0;
	for (size_t i = 0; i < 4; i++)
	{
		double a;
		cin >> a;
		x0.push_back(a);
	}
	cout << "Верхний диапазон для координат: ";
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
	f.set_coordinates(l1, d1, l2, d2);
	f.eval_static_vecs();
	f.solve_xi();
	cout << "Задайте значения x1, x2 для поля смещения\nx1 = ";
	double xx1, xx2;
	cin >> xx1;
	cout << "x2 = "; cin >> xx2;

	view_points = { xx1, xx2 };
	true_u.resize(view_points.size());
	true_u = { f.number_field(xx1), f.number_field(xx2) };

	auto start_time = std::chrono::high_resolution_clock::now();
	auto res = genetic_alg(resid_func, eps, x0, x1);
	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end_time - start_time;
	std::cout << "\nВремя работы алгоритма в сек: " << duration.count() << std::endl << endl;
	
	printf("Найденные параметры: %.5f; %.5f; %.5f; %.5f\n", res[0], res[1], res[2], res[3]);
	cout << "Значение функции: " << resid_func(res) << endl;

	system("pause");
}

void Minimize_with_Nelder_Mid(Fissures f, const vector<double> point, const double eps, const double l) {
	
	f.eval_static_vecs();
	f.solve_xi();
	true_u = { f.number_field(view_points[0]), f.number_field(view_points[1])};
	cout << resid_func(point) << endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	//cout << resid_func({ f.l1, f.d1, f.l2, f.d2 }) << endl;
	auto res = nelder_mead(resid_func, point, l, eps);
	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end_time - start_time;
	std::cout << "Время работы алгоритма: " << duration.count() << std::endl;
	printf("Найденные параметры: %.5f; %.5f; %.5f; %.5f\n",
		res[0], res[1], res[2], res[3]);
	cout << endl;
	cout << "Значение функции: " << resid_func(res) << endl;
}

double fun(vector<double> x) {
	return x[0] * x[0] + x[0] * x[1] + x[1] * x[1] - 6 * x[0] - 9 * x[1];
}

void task4(const double l1, const double d1, const double l2, const double d2, const double k,
	const vector<double> floor, const vector<double> roof, const double eps1,
	const double eps2, const double l,
	string field_file_name, string report_file_name
) {

	Fissures f = Fissures(l1, l2, d1, d2, 20, k);
	f.eval_static_vecs();
	f.solve_xi();
	auto x = arma::linspace(0.0, 10.0, 1000);
	//Вычисляем поле
	cx_vec u_vals(x.size());
	for (size_t i = 0; i < x.size(); i++)
	{
		u_vals(i) = f.number_field(x(i));
	}
	write_to_file(u_vals, f.l1, f.d1, f.l2, f.d2, f.k, field_file_name);

	cout << "Введите точки наблюдения: ";
	view_points.clear();
	while (true) {
		double a;
		cin >> a;
		if (a == 0)
			break;
		view_points.push_back(a);
	}

	true_u.resize(view_points.size());
	for (size_t i = 0; i < view_points.size(); i++)
		true_u[i] = f.number_field(view_points[i]);

	//cout << resid_func({ 0.0998931, 0.499328, 0.100106, 0.999252 }) << endl;
	//return;
	//находим с помощью генетического
	cout << "Начало работы генетического алгоритма" << endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	auto res1 = genetic_alg(resid_func, eps1, floor, roof);
	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end_time - start_time;
	std::cout << "Время работы алгоритма: " << duration.count() << std::endl;
	printf("Найденные параметры с помощью генетического алгоритма: %.5f; %.5f; %.5f; %.5f\n",
		res1[0], res1[1], res1[2], res1[3]);
	cout << endl;
	cout << "Значение функции: " << resid_func(res1) << endl;

	//уточняем с помощью Нелдера-Мида
	cout << "Начало работы алгоритма Нелдера-Мида" << endl;
	start_time = std::chrono::high_resolution_clock::now();
	auto res2 = nelder_mead(resid_func, res1, l, eps2);
	end_time = std::chrono::high_resolution_clock::now();
	duration = end_time - start_time;
	std::cout << "Время работы алгоритма: " << duration.count() << std::endl;
	printf("Найденные параметры: %.5f; %.5f; %.5f; %.5f\n",
		res2[0], res2[1], res2[2], res2[3]);
	cout << endl;
	cout << "Значение функции: " << resid_func(res2) << endl;

	std::ofstream file(report_file_name, std::ios::app);
	if (!file.is_open()) {
		std::cerr << "Не удалось открыть файл для записи!" << std::endl;
		return;
	}
	file << endl<< "Истинные параметры: " << endl;
	file << "l1 = " << l1 << endl;
	file << "d1 = " << d1 << endl;
	file << "l2 = " << l2 << endl;
	file << "d2 = " << d2 << endl;
	file << "Найденное с помощью генетического: " << endl;
	file << "l1 = " << res1[0] << endl;
	file << "d1 = " << res1[1] << endl;
	file << "l2 = " << res1[2] << endl;
	file << "d2 = " << res1[3] << endl;
	file << "Найденное с помощью Нелдера-Мида: " << endl;
	file << "l1 = " << res2[0] << endl;
	file << "d1 = " << res2[1] << endl;
	file << "l2 = " << res2[2] << endl;
	file << "d2 = " << res2[3] << endl;
	file << "Доп. информация: " << endl;
	file << "Точки наблюдения:";
	for (size_t i = 0; i < view_points.size(); i++)
		file << " x" << i + 1 << " = " << view_points[i];
	file << endl;
	file << "Нижний порог параметров: " << floor[0] << ", " << floor[1] << ", " << floor[2] << ", " << floor[3] << endl;
	file << "Верхний порог параметров: " << roof[0] << ", " << roof[1] << ", " << roof[2] << ", " << roof[3] << endl;
	file << "eps генетического = " << eps1 << endl;
	file << "eps Н.-М. = " << eps2 << endl;
	file << "Длина симлекса = " << l << endl;
	file << "Частота k = " << k << endl;
	file << "Значение функции резистентности на генетическом = " << resid_func(res1) << endl;
	file << "Значение функции резистентности на Н.-М. = " << resid_func(res2) << endl << endl;

	file.close();
}

int main() {
	setlocale(LC_ALL, "Russian");
	
	Fissures f = Fissures();
	double l1, l2, d1, d2;
	l1 = 0.1;
	d1 = 0.5;
	l2 = 0.1;
	d2 = 1.0;
	f.set_coordinates(l1, d1, l2, d2);
	vector<double> roof = { 0.1, 3.0, 0.1, 3.0 };
	vector<double> floor = { 0.0, 0.0, 0.0, 0.0 };
	double k = 5, l = 0.1;
	double eps1 = 1e-7, eps2 = 1e-10;
	view_points = { 2.0, 2.7 };
	

	task4(l1, d1, l2, d2, k, floor, roof, eps1, eps2, l,
		std::format("D:\\VS Projects\\Two Fissures\\results\\task4\\k{}l1{}d1{}l2{}d2{}.csv", k, l1, d1, l2, d2),
		"D:\\VS Projects\\Two Fissures\\results\\report4.0.txt");

	//view_points = { 2.0, 2.7 };
	
	
	/*f.fill_k3_integral();
	f.fill_sigma2_alpha();
	f.solve_xi();
	f.fill_u_sigma();
	true_u = { f.number_field(view_points[0]), f.number_field(view_points[1]) };
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 10; j++)
		{
			switch (i)
			{
			case 0:
				l1 = 0.1 + j * 0.01;
				break;
			case 1:
				d1 = 0.5 + j * 0.01;
				break;
			case 2:
				l2 = 0.1 + j * 0.01;
				break;
			case 3:
				d2 = 1 + j * 0.01;
				break;
			default:
				break;
			}
			printf("l1 = %.5f, d1 = %.5f, l2 = %.5f, d2 = %.5f, resid_fun = %.5f\n", l1, d1, l2, d2, resid_func({l1, d1, l2, d2}));
			l1 = 0.1; d1 = 0.5; l2 = 0.1; d2 = 1;
		}
	}*/
	//Minimize_with_Nelder_Mid(f, { 0.0864128, 0.332124, 0.112353, 0.933864 }, 1e-8, 0.1);

	return 0;
}