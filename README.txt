Данный вычислительный модуль содержит класс Fissures, который решает прямую геометрическую задачу. Этот класс вычисляет функцию расслоения и волновое поле смещения.
Конструктор Fissures(double l1, double l2, double d1, double d2, int N = 20, double k = 5, int p = 1); позволяет задать параметры моделируемой трещины. Кроме этого конструкта можно поменять параметры трещины не меняя прочих значений используя метод Fissures::set_parameters(l1, d1, l2, d2).

После создания экземпляра класса и перед началом самой минимизации используйте метод Fissures::evalf_static_vecs() а затем метод Fissures::solve_xi(). Первый метод заполнить общие данные, которые не меняются во время оптимизации, а второй метод вычислит функцию расслоения. Не пытайтесь найти поле смещения, не вызвав эти методы, а также не меняйте их порядок вызова.

Порядок вызова методов, если обращение к классу происходит впервые:
1. Конструктор
2. Fissures::evalf_static_vecs()
3. Fissures::solve_xi()
4. Поле смещения или прочие действия

В дальнейшем при работе программы рекомендуется не вызывать метод Fissures::evalf_static_vecs(), но при создании нового экземпляра класса необходимо будет снова вызвать метод Fissures::solve_xi().

За вычисление поля смещения в конкретной точке отвечает функция Fissures::number_field(const double x). Она принимает вещественную точку, в которой надо вычислить поле и возвращает переменную типа cx_double. 

Также класс Fissure содержит метод write_field_to_csv(std::string file_name), который записывает данные поля смещения. В файл попадают следующие данные:
l1
d1
l2
d2
k
1000 значений поля смещения от 0 до 10 перечисленные через ';'
Рекомендуется назначать файл формата .csv для корректной обработки скриптом, создающем график волнового поля по указанным данным.

в заголовке "Algorithms.h" содержаться класс Genetic_Alg и функция nelder_mid. Конструктор Genetic_Alg  принимает минимизируемую функцию f, вектор lower_bound элементов типа double, которые являются нижней допустимой границей переменных функции f, 
Genetic_Alg(function<double(vector<double>)> f, vector<double> lower_bound, vector<double> upper_bound, double eps);
f - минимизируемая функция
lower_bound - нижняя граница допустимых значений переменных
upper_bound - верхняя граница допустимых значений переменных
eps - точность

vector<double> nelder_mead(const function<double(vector<double>)> f, const vector<double> x0, const double l, const double eps,
	const double alpha = 1.0, const double gamma = 2.0, const double rho = 0.5, const double sigma = 0.5, const int max_iter = 200);
f - минимизируемая функция
x0 - начальная точка, вокруг которой строится симплекс
l - длина симплекса
eps - точность

Файл main.cpp содержит функцию резистентности resid_func(vector<double> parameters). Именно она минимизируется в коде. parameters передаются в следующем порядке: {l1, d1, l2, d2}

###Быстрое начало
В файле main.cpp находится функция task4, которая решает обратную задачу, записывает данные поля смещения и результат проведенного эксперимента в указанных файлах.

void task4(const double l1, const double d1, const double l2, const double d2, const double k,
	const vector<double> floor, const vector<double> roof, const double eps1,
	const double eps2, const double l,
	string field_file_name, string report_file_name)
l1, d1, l2, d2 - параметры
k - частота
floor - нижняя граница допустимых значений параметров
roof - верхняя граница допустимых значений параметров
eps1 - точность генетического алгоритма
eps2 - точность алгоритма нелдера-мида
l - длина стороны симплекса 
field_file_name - имя файла с данными поля смещения
report_file_name - имя файла, в который запишутся результаты эксперимента

Отдельно в функции main задаются точки наблюдения в переменную view_points

Функция task4 решает прямую задачу для трещин с заданными параметрами, записывает данные поля смещения, решает обратную задачу в указанных точках наблюдения сначала с помощью генетического алгоритма а потом с помощью алгоритма нелдера-мида и записывает результат эксперимента в указанный файл.

Примечание. На результат минимизации также сильно влияет и длина симплекса. Если результат кажется неудовлетворительным, то стоит попробовать вызвать метод Minimize_with_Nelder_Mid задав в качестве начальной точки точку, полученную на генетическом алгоритме, но изменив l(длина симплекса).
