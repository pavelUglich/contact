#include "layer.hpp"
#include <numeric>
#include"plots.h"
#include "Integral.h"
#include "constant_element.h"
#include "linear_element.h"



const double pi = acos(-1.0);

template<class Fun>
map<double, vector<complex<double>>> dispersional_set(
	double the_highest_frequency, Fun rho, Fun mu)
{
	map<double, vector<complex<double>>> result;
	double kappa = 0.1;
	const double step = 0.05;
	layer l(mu, rho, kappa, 5);
	auto roots = l.get_roots();
	while (kappa < the_highest_frequency)
	{
		cout << kappa << endl;
		result[kappa] = roots;
		kappa += step;
		l = { mu, rho, kappa, 7 };
		roots = l.get_roots();
	}
	return result;
}

void plot_the_dispersional_curves(map<double, 
	vector<complex<double>>>& dispersion_set, const string& file_name)
{
	ofstream os(file_name);
	os << "\\begin{tikzpicture}[scale=1.5]\n";
	os << "\\begin{axis}[grid]\n";
	const auto number_of_curves = dispersion_set.rbegin()->second.size();
	for (auto [x, y] : dispersion_set) {
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


void surface_plot(const std::vector<double>& x_1, const std::vector<double>& x_2,
	const std::vector<std::complex<double>>& solution, const layer& l,
	const std::string& file_name)
{
	const auto matrix = l.displacement(x_1, x_2, solution);
	std::ofstream ofs(file_name);
	std::ranges::for_each(x_1.begin(), x_1.end(), [&](auto x) {ofs << x << " "; });
	ofs << std::endl;
	std::ranges::for_each(x_2.begin(), x_2.end(), [&](auto x) {ofs << x << " "; });
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

void matrix_t_out(const std::vector<std::vector<double>>& u, std::ofstream& ofs)
{
	for (size_t i = 0; i < u.front().size(); i++)
	{
		std::ranges::for_each(u.begin(), u.end(),
			[&ofs, i](auto x) {ofs << x[i] << " "; });
		ofs << std::endl;
	}
}

void umoff_pointing1(const std::vector<double>& x_1, const std::vector<double>& x_2,
	const std::vector<std::complex<double>>& solution, const layer& l,
	const std::string& file_name)
{
	std::ofstream ofs(file_name);
	std::ranges::for_each(x_1.begin(), x_1.end(), [&ofs](auto x) {ofs << x << " "; });
	ofs << std::endl;
	std::ranges::for_each(x_2.begin(), x_2.end(), [&ofs](auto x) {ofs << x << " "; });
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
			auto norm = abs(up[0]) * abs(up[0]) + abs(up[1]) * abs(up[1]);
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


void umov_pointing(const std::vector<double>& x_1, double x_2,
	const std::vector<std::complex<double>>& solution, const layer& l,
	const std::string& file_name) {
	std::ofstream ofs(file_name);
	ofs << "\\begin{tikzpicture} [scale = 2.1] \\tiny\n";
	ofs << "\\begin{axis}[grid]\\n";
	for (auto xx : x_1)
	{
		auto x = xx;
		auto y = x_2;
		auto hh = 0.025;
		const size_t limit = 600;
		size_t counter = 0;
		ofs << "\\addplot[line width = 0.25mm, smooth, red] plot coordinates{\n";
		ofs << "(" << x << ", " << y << ") ";
		while (abs(x) < 2.0)
		{

			auto uu = l.umoff_pointing(x, y, solution);
			const auto norm = sqrt(uu[0] * uu[0] + uu[1] * uu[1]);
			uu = uu / norm;
			auto new_y = y + uu[1] * hh;
			while (new_y < 0 || new_y > 1)
			{
				hh /= 2;
				new_y = y + uu[1] * hh;
			}
			x += uu[0] * hh;
			y += uu[1] * hh;
			if (counter++ > limit)
			{
				break;
			}
			ofs << "(" << x << ", " << y << ") ";
		}
		ofs << "};\n";
	}
	ofs << "\\end{axis}\n";
	ofs << "\\end{tikzpicture}";
	ofs.close();
}

std::vector<double> make_nodes(size_t quantity, double a) {
	const auto fun = [=](double x) {return a * a - x * x; };
	std::map<size_t, double> lengths;
	const double h = 2.0 * a / quantity;
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
	}
	return nodes;
}

template<class Elem>
vector<boundary_element*> make_elements(const vector<double>& nodes) {
	std::vector<boundary_element*> elements(nodes.size() - 1);
	for (size_t i = 0; i < elements.size(); i++)
	{
		elements[i] = new Elem(nodes[i], nodes[i + 1]);
	}
	return elements;
}



int main() {
	const auto mu = [](double x) {return 1.0 + x; };
	const auto rho = [](double x) {return 1.0 + exp(-x); };
	const auto mass = 1.0;// масса штампа
	const auto kappa = 3.0; // волновое число
	const auto tt = 1.0; // интенсивность нагрузки	
	const double y = 1.0;
	size_t size = 20;// 0;
	std::vector<std::complex<double>> wavefield(size);
	std::vector<std::complex<double>> wavefield1(size + 1);
	const double a = 1.0;
	const auto h = 2 * a / (double)size;
	double left = -a;
	vector<complex<double>> sol1;
	const auto nodes = make_nodes(size, a);
	const auto elements = make_elements<linear_element>(nodes);
	const auto elements_ = make_elements<constant_element>(nodes);
	/*
		{
			std::vector<std::complex<double>> fr;
			for (size_t i = 0; i < 80; i++)
			{
				const double kappa = (i + 1) * 0.1;
				std::cout << kappa << std::endl;
				for (size_t i = 0; i < size; i++)
				{
					elements[i] = new linear_element(nodes[i], nodes[i + 1]);
				}
				layer l(mu, rho, kappa, elements, 40, 0.1e-6);
				const auto mmm = l.matrix();
				auto matrix = SquareMatrix<complex<double>>(size + 1, mmm);
				auto sol = matrix.linSolve(vector<complex<double>>(size + 1, { 1,0 }));
				const auto lengths = l.lengths();
				complex<double> integral_ = 0.5 * (sol.front() * lengths.front() + sol.back() * lengths.back());
				for (size_t i = 1; i < lengths.size() - 1; i++)
				{
					integral_ += 0.5 * sol[i] * (lengths[i - 1] + lengths[i]);
				}
				auto gamma = tt / (integral_ - mass * kappa * kappa);
				fr.push_back(gamma);
			}
			plotTheWaveField(0.1, { {"black", fr} }, "fr.txt", 0.1);
		}//*/

	{
		layer l(mu, rho, kappa, elements_, 80, 0.1e-6); //*/ (mu, rho, "10.txt", elements_); // = { mu, rho, "roots.txt", elements };
		const auto lengths = l.lengths();
		const auto weight = [=](double x) {return 1.0 / sqrt(a * a - x * x); };
		auto mmm = l.matrix();//(weight);
		auto matrix = SquareMatrix<complex<double>>(mmm.size(), mmm);
		auto sol = matrix.linSolve(vector<complex<double>>(mmm.size(), { 1,0 }));
		const auto nodes = get_nodal_points(elements_);
		for (size_t i = 0; i < sol.size(); i++)
		{
			sol[i] *= weight(nodes[i]);
		}
		complex<double> integral_ = 0.5 * (sol.front() * lengths.front() + sol.back() * lengths.back());
		for (size_t i = 1; i < lengths.size() - 1; i++)
		{
			integral_ += 0.5 * sol[i] * (lengths[i - 1] + lengths[i]);
		}
		auto gamma = tt / (integral_ - mass * kappa * kappa);
		std::ranges::transform(sol.begin(), sol.end(), sol.begin(), [=](auto x) {return x * gamma; });
		plotTheWaveField(nodes, { {"red", sol} }, "sol.txt");
		/*
		layer l1/*(mu, rho, kappa, elements, 50, 0.1e-6);//(mu, rho, "10.txt", elements_); // = { mu, rho, "roots.txt", elements };
		mmm = l1.matrix();
		matrix = SquareMatrix<complex<double>>(mmm.size(), mmm);
		sol1 = matrix.linSolve(vector<complex<double>>(mmm.size(), { 1,0 }));
		integral_ = 0;// 0.5 * (sol.front() * lengths.front() + sol.back() * lengths.back());
		for (size_t i = 0; i < lengths.size(); i++)
		{
			integral_ += sol1[i] * lengths[i];
		}
		gamma = tt / (integral_ - mass * kappa * kappa);//(integral_ - tt) / mass / kappa / kappa;
		std::ranges::transform(sol1.begin(), sol1.end(), sol1.begin(), [=](auto x) {return x * gamma; });
		plotTheWaveField(get_nodal_points(elements), { {"red", sol} }, "sol.txt");
		//*/

		/*
		std::vector<double> points_x(80);
		const auto h_x = 2 * a / points_x.size();
		for (size_t i = 1; i <= points_x.size(); i++)
		{
			points_x[i - 1] = (i - 0.5) * h_x;
		}

		std::vector<double> points_y(40);
		const auto h_y = 1.0 / points_y.size();
		for (size_t i = 1; i <= points_y.size(); i++)
		{
			points_y[i - 1] = (i - 0.5) * h_y;
		}*/
		//const auto displacement = l.displacement(points_x, points_y, sol1);
		//surface_plot(points_x, points_y, sol, l, "sp.txt");
		//umov_pointing(points_x, 0.99, sol, l, "up.txt");
		//umoff_pointing1(points_x, points_y, sol, l, "up.txt");//*/

	}//*/
	//plotTheWaveField(h, { {"red", wavefield}, {"blue", wavefield1} }, "wavefield.txt", -a + h / 2);//*/
	system("pause");
}