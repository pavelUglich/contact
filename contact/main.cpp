#include "layer.hpp"
#include <numeric>
#include"plots.h"
#include "Integral.h"
#include "constant_element.h"
#include "linear_element.h"

const double pi = acos(-1.0);


const double L = 4.065e-6;
const auto first_law = [](double x) {return 210.0 + 430.0 * x / L; };

const double x2_star = 0.9 * L;
const double A2 = (-40000.0 + 295900.0 / 8.1) / L / L / L;
const double B2 = -29590 / (8.1 * L * L);

const auto second_law = [=](double x) {
	const auto xi = x - x2_star;
	return A2 * xi * xi * xi + B2 * xi * xi + 640;
	};

const double x3_star = 0.7 * L;
const double A3 = (-40000.0 + 2279700.0 / 44.1) / 27.0 / L / L / L;
const double B3 = -25330.0 / 44.1 / L / L;
const auto third_law = [=](double x) {
	const auto xi = x - x3_star;
	return A3 * xi * xi * xi + B3 * xi * xi + 640;
	};

const double nu = 0.3;
const double rho = 7.8e3;

const auto mu_1 = [=](double x) {
	const double xi = x * L;
	return first_law(xi) / 2.0 / (1.0 + nu);
	};

const auto mu_2 = [=](double x) {
	const double xi = x * L;
	return second_law(xi) / 2.0 / (1.0 + nu);
	};

const auto mu_3 = [=](double x) {
	const double xi = x * L;
	return third_law(xi) / 2.0 / (1.0 + nu);
	};

const auto mu_0 = mu_1(0.0);

const auto mu_1_dimentionless = [=](double x) {
	return mu_1(x) / mu_0;
	};

const auto mu_2_dimentionless = [=](double x) {
	return mu_2(x) / mu_0;
	};

const auto mu_3_dimentionless = [=](double x) {
	return mu_3(x) / mu_0;
	};


/// <summary>
/// построение дисперсионного множества
/// </summary>
/// <typeparam name="Fun">дисперсионное уравнение</typeparam>
/// <param name="the_highest_frequency">наибольшая частота</param>
/// <param name="rho">плотность</param>
/// <param name="mu">модуль сдвига</param>
/// <returns>ассоциативный контейнер</returns>
template<class Fun>
map<double, vector<complex<double>>> dispersional_set(
	double the_highest_frequency, Fun rho, Fun c_44, Fun c_55, bottom_c bc = CLAMPED)
{
	map<double, vector<complex<double>>> result;
	double kappa = 0.1;
	const double step = 0.05;
	layer l = { c_44, c_55, rho, kappa, 5, 1.0e-5, bc };
	auto roots = l.get_roots();
	while (kappa < the_highest_frequency)
	{
		cout << kappa << endl;
		result[kappa] = roots;
		kappa += step;
		l = { c_44, c_55, rho, kappa, 7, 1.0e-6, bc };
		roots = l.get_roots();
	}
	return result;
}

/// <summary>
/// построение дисперсионного множества на графике
/// </summary>
/// <param name="dispersion_set">дисперсионное множество</param>
/// <param name="file_name">имя файла</param>
void plot_the_dispersional_curves(map<double,
	vector<complex<double>>>& dispersion_set, const string& file_name)
{
	ofstream os(file_name);
	os << "\\begin{tikzpicture}[scale=1.5]\n";
	os << "\\begin{axis}[grid]\n";
	const auto number_of_curves = dispersion_set.rbegin()->second.size();
	for (auto& [x, y] : dispersion_set) {
		for (auto& i : y) {
			if (i.real() < DBL_EPSILON)
				i = { -i.imag(),0 };
		}
		std::ranges::sort(y.begin(), y.end(),
			[](auto xx, auto yy) { return xx.real() < yy.real(); });
	}
	for (size_t i = 0; i < number_of_curves; i++) {
		os << "\\addplot[smooth, black] plot coordinates{\n";
		for (auto [x, y] : dispersion_set) {
			if (i < y.size())
			{
				os << "(" << y[i].real() << ", " << x << ") ";
			}
		}
		os << "};\n";
	}
	os << "\\end{axis}\n";
	os << "\\end{tikzpicture}\n";
	os.close();
}

/// <summary>
/// создание текстового файла для построения диаграммы
/// </summary>
/// <param name="x_1">вектор координат по x</param>
/// <param name="x_2">вектор координат по y</param>
/// <param name="solution">решение контактной задачи</param>
/// <param name="l">слой</param>
/// <param name="file_name">имя файла</param>
void surface_plot(const std::vector<double>& x_1, const std::vector<double>& x_2,
	const std::vector<std::complex<double>>& solution, const layer& l,
	const std::string& file_name)
{
	const auto matrix = l.displacement(x_1, x_2, solution);
	std::ofstream ofs(file_name);
	std::ranges::for_each(x_1.begin(), x_1.end(),
		[&ofs](auto x) {ofs << x << " "; });
	ofs << std::endl;
	std::ranges::for_each(x_2.begin(), x_2.end(),
		[&ofs](auto x) {ofs << x << " "; });
	ofs << std::endl;
	for (size_t i = 0; i < x_2.size(); i++)
	{
		for (size_t ii = 0; ii < x_1.size(); ii++)
		{
			ofs << abs(matrix[ii][i + 1]) << " ";
		}
		ofs << std::endl;
	}
	ofs.close();
}

/// <summary>
/// вывод транспонированной матрицы в файл
/// </summary>
/// <param name="u">матрица</param>
/// <param name="ofs">поток вывода</param>
void matrix_t_out(const std::vector<std::vector<double>>& u,
	std::ofstream& ofs)
{
	for (size_t i = 0; i < u.front().size(); i++)
	{
		std::ranges::for_each(u.begin(), u.end(),
			[&ofs, i](auto x) {ofs << x[i] << " "; });
		ofs << std::endl;
	}
}

/// <summary>
/// построение волнового поля вектора Умова-Пойнтинга
/// </summary>
/// <param name="x_1">вектор координат по x</param>
/// <param name="x_2">вектор координат по y</param>
/// <param name="solution">решение контактной задачи</param>
/// <param name="l">слой</param>
/// <param name="file_name">имя текстового файла</param>
void umoff_pointing1(const std::vector<double>& x_1,
	const std::vector<double>& x_2,
	const std::vector<std::complex<double>>& solution, const layer& l,
	const std::string& file_name)
{
	std::ofstream ofs(file_name);
	std::ranges::for_each(x_1.begin(), x_1.end(),
		[&ofs](auto x) {ofs << x << " "; });
	ofs << std::endl;
	std::ranges::for_each(x_2.begin(), x_2.end(),
		[&ofs](auto x) {ofs << x << " "; });
	ofs << std::endl;
	std::vector<std::vector<double>> u(x_1.size(),
		std::vector<double>(x_2.size()));
	std::vector<std::vector<double>> v(x_1.size(),
		std::vector<double>(x_2.size()));
	for (size_t i = 0; i < x_2.size(); i++)
	{
		for (size_t ii = 0; ii < x_1.size(); ii++)
		{
			auto up = l.umoff_pointing(x_1[ii], x_2[i], solution);
			auto norm = up[0] * up[0] + up[1] * up[1];
			norm = sqrt(norm);
			up = up / norm;
			u[ii][i] = up[0];
			v[ii][i] = up[1];
		}
	}
	matrix_t_out(u, ofs);
	matrix_t_out(v, ofs);
	ofs.close();
}

/// <summary>
/// создаёт набор узловых точек
/// </summary>
/// <param name="quantity">количество</param>
/// <param name="a">граница отрезка [-a,a]</param>
/// <returns>вектор узловых точек</returns>
std::vector<double> make_nodes(size_t quantity, double a) {
	/*
	const auto fun = [=](double x) {return a * a - x * x; };
	std::map<size_t, double> lengths;
	const double h = 2.0 * a / static_cast<double>(quantity);
	double sum = 0;
	for (size_t i = 0; i < quantity; i++)
	{
		lengths[i] = fun(-1 + (i + 0.5) * h);
		sum += lengths[i];
	}
	sum /= 2 * a;
	for (auto& [x, y] : lengths) {
		y /= sum;
	}
	std::vector<double> nodes(quantity + 1);
	nodes[0] = -a;
	for (size_t i = 0; i < quantity; i++)
	{
		nodes[i + 1] = nodes[i] + lengths[i];
	}//*/
	const auto h = 2 * a / quantity;
	std::vector<double> nodes(quantity + 1);
	nodes[0] = -a;
	for (size_t i = 0; i < quantity; i++)
	{
		nodes[i + 1] = nodes[i] + h;
	}//*/
	return nodes;
}

/// <summary>
/// создаёт набор элементов
/// </summary>
/// <typeparam name="Elem">тип элемента</typeparam>
/// <param name="nodes">вектор узловых точек</param>
/// <returns>вектор указателей на элементы</returns>
template<class Elem>
vector<boundary_element*> make_elements(const vector<double>& nodes) {
	std::vector<boundary_element*> elements(nodes.size() - 1);
	for (size_t i = 0; i < elements.size(); i++)
	{
		elements[i] = new Elem(nodes[i], nodes[i + 1]);
	}
	return elements;
}


std::vector<std::complex<double>> get_solution(layer& l, double a, size_t size,
	double mass = 1.0, double tt = 1.0) {
	const auto points = make_nodes(size, a);
	const auto elements = make_elements<constant_element>(points);
	l.set_elements(elements);
	const auto weight = [a](double x) {return 1.0 / sqrt(a * a - x * x); };
	const auto m = l.matrix(weight);
	const auto matrix = square_matrix<complex<double>>(m.size(), m);
	auto sol = matrix.linSolve(vector<complex<double>>(m.size(), { 1,0 }));
	const auto lengths = l.lengths();
	const auto nodes = get_nodal_points(elements);
	for (size_t i = 0; i < nodes.size(); i++)
	{
		sol[i] *= weight(nodes[i]);
	}
	const auto integral = l.evaluate_integral(sol);
	const auto kappa = l.get_kappa();
	const auto gamma = tt / (integral - mass * kappa * kappa);
	std::ranges::transform(sol.begin(), sol.end(), sol.begin(),
		[gamma](auto x) {return x * gamma; });
	return sol;
}


std::complex<double> get_displacement(layer& l, double a, size_t size,
	double mass = 1.0, double tt = 1.0) {
	const auto points = make_nodes(size, a);
	const auto elements = make_elements<constant_element>(points);
	l.set_elements(elements);
	const auto weight = [a](double x) {return 1.0 / sqrt(a * a - x * x); };
	const auto m = l.matrix(weight);
	const auto matrix = square_matrix<complex<double>>(m.size(), m);
	auto sol = matrix.linSolve(vector<complex<double>>(m.size(), { 1,0 }));
	const auto lengths = l.lengths();
	const auto nodes = get_nodal_points(elements);
	for (size_t i = 0; i < nodes.size(); i++)
	{
		sol[i] *= weight(nodes[i]);
	}
	const auto integral = l.evaluate_integral(sol);
	const auto kappa = l.get_kappa();
	return tt / (integral - mass * kappa * kappa);
}

void solve_with_constant_elems(layer& l, double a, size_t size,
	const std::string& file_name, double mass = 1.0, double tt = 1.0, const std::string& colour = "red")
{
	const auto points = make_nodes(size, a);
	const auto elements = make_elements<constant_element>(points);
	const auto sol = get_solution(l, a, size, mass, tt);
	plot_the_wave_field(get_nodal_points(elements), { {colour, sol} }, file_name);
}

const std::vector<double> T2 = { 0.17 / 9 / pi * 1e-3, 3.25 / 9 / pi * 1e-3, 7.0 / 9 / pi * 1e-3 };
const std::vector<std::string> colours = { "black", "red", "blue" };

int main() {
	
	{
		const auto c_44 = [](double x) {return mu_1_dimentionless(x); };
		const auto c_55 = [](double x) {return mu_1_dimentionless(x); }; //
		const auto density = [](double x) {return 1.0; };
		// дисперсионные кривые
		// auto ds = dispersional_set<std::function<double(double)>>(9, density, c_44, c_55);
		// plot_the_dispersional_curves(ds, "ds_law1.txt");
		layer l = { c_44, density, 6.0 };
		for (size_t i = 0; i < T2.size(); i++)
		{
			solve_with_constant_elems(l, 1.0, 20, "law" + std::to_string(i) + ".txt", 1.0, T2[i]);
		}
	}
/*
	{
		const auto c_44 = [](double x) {return mu_2_dimentionless(x); };
		const auto c_55 = [](double x) {return mu_2_dimentionless(x); }; //
		const auto density = [](double x) {return 1.0; };
		//	дисперсионные кривые
		// auto ds = dispersional_set<std::function<double(double)>>(9, density, c_44, c_55);
		// plot_the_dispersional_curves(ds, "ds_law2.txt");
		layer l = { c_44, density, 4.0 };
		for (size_t i = 0; i < T2.size(); i++)
		{
			solve_with_constant_elems(l, 1.0, 20, "law" + std::to_string(i) + ".txt", 1.0, T2[i]);
		}

	}*/
	/*
	{
		const auto c_44 = [](double x) {return mu_3_dimentionless(x); };
		const auto c_55 = [](double x) {return mu_3_dimentionless(x); }; //
		const auto density = [](double x) {return 1.0; };
		layer l = { c_44, density, 7.0 };
		for (size_t i = 0; i < T2.size(); i++)
		{
			solve_with_constant_elems(l, 1.0, 20, "law" + std::to_string(i) + ".txt", 1.0, T2[i], colours[i]);
		}

		//	дисперсионные кривые
		//auto ds = dispersional_set<std::function<double(double)>>(9, density, c_44, c_55);
		//plot_the_dispersional_curves(ds, "ds_law3.txt");
	}//*/
}