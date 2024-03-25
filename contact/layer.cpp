#include "OdeSolver.h"
#include "layer.hpp"
#include <numeric>
#include <fstream>
#include <iomanip>
#include <set>
#include <map>


vector<complex<double>> layer::function_and_derivative(complex<double> alpha) const {
	const auto sol1 = bigger_system(alpha).solve(0, 1, { 0,1,0,0 });
	return  { sol1[1], sol1[3] };
}


complex<double> layer::newton_method(complex<double> begin_param) const
{
	complex<double> xx;
	do {
		const auto fd = function_and_derivative(begin_param);
		xx = fd[0] / fd[1];
		begin_param -= xx;
	} while (abs(xx) > _tolerance);
	return begin_param;
}

double layer::equation(complex<double> alpha) const {
	return smaller_system(alpha).solve(0, 1, { 0,1 })[1].real();
}

void layer::roots(size_t roots_number) {
	_roots.resize(roots_number);
	_residues.resize(roots_number);
	double a = 1.5 * _kappa;
	double step = a / 10;
	auto b = a - step;
	size_t counter = 0;
	for (size_t i = 0; i < 10; i++) {
		const auto f_a = equation({ a, 0 });
		const auto f_b = equation({ b, 0 });
		if (f_a * f_b < 0) {
			_roots[counter] = newton_method({ std::midpoint(a, b),0 });
			const auto fad = function_and_derivative(_roots[counter]);
			_residues[counter] = fad.back();
			++counter;
		}
		a = b;
		b -= step;
	}

	a = 0;
	b = step;
	while (true) {
		const auto f_a = equation({ 0, a });
		const auto f_b = equation({ 0, b });
		if (f_a * f_b < 0) {
			_roots[counter] = newton_method({ 0, std::midpoint(a, b) });
			const auto fad = function_and_derivative(_roots[counter]);
			_residues[counter] = fad.back();
			++counter;
		}
		a = b;
		b += step;
		if (counter >= _roots.size())
			break;
	}
}

OdeSolver<complex<double>> layer::smaller_system(complex<double> alpha) const
{
	return
	{
		{
			[this](double x, const vector<complex<double>>& y) {return y[1] / _mu(x); },
			[this, alpha](double x, const vector<complex<double>>& y)
			{
				return (alpha * alpha * _mu(x) - _rho(x) * _kappa * _kappa) * y[0];
			}
		},
		_tolerance,
		RUNGE_KUTTA_FELDBERG
	};
}

OdeSolver<complex<double>> layer::bigger_system(complex<double> alpha) const
{
	return
	{
		{
			[this](double x, const vector<complex<double>>& y) {return y[1] / _mu(x); },
			[this, alpha](double x, const vector<complex<double>>& y)
			{
				return (alpha * alpha * _mu(x) - _rho(x) * _kappa * _kappa) * y[0];
			},
			[this](double x, const vector<complex<double>>& y) {return y[3] / _mu(x); },
			[this, alpha](double x, const vector<complex<double>>& y)
			{
				return  2.0 * alpha * _mu(x) * y[0]
					+ (alpha * alpha * _mu(x) - _rho(x) * _kappa * _kappa) * y[2];
			}
		},
		_tolerance,
		RUNGE_KUTTA_FELDBERG
	};
}

vector<complex<double>> layer::get_residues() const
{
	return this->_residues;
}

double layer::get_tolerance() const
{
	return this->_tolerance;
}


vector<complex<double>> layer::row(double node, size_t size, double y) const
{
	vector<complex<double>> result(size);
	for (size_t i = 0; i < _roots.size(); i++)
	{
		const auto multiplier = smaller_system(_roots[i]).solve(0, y, { 0,1 })
			.front() / _residues[i];
		const auto r = row(node, _roots[i]);
		result = result + I * multiplier * r;
	}
	return result;
}

vector<vector<complex<double>>> layer::row(double node, size_t size,
	const vector<double>& y) const
{
	auto points = y;
	process_points_vector(points);
	vector<vector<complex<double>>> result(points.size(),
		vector<complex<double>>(size));
	for (size_t i = 0; i < _roots.size(); i++)
	{
		const auto solution = smaller_system(_roots[i]).solve(points, { 0,1 });
		vector<complex<double>> multiplier(solution.size());
		for (size_t ii = 0; ii < solution.size(); ii++)
		{
			multiplier[ii] = I * solution[ii].front() / _residues[i];
		}
		const auto r = row(node, _roots[i]);
		const auto term = multiplier % r;
		result = result + term;
	}
	return result;
}

vector<vector<vector<complex<double>>>> layer::row(const vector<double>& nodes,
	size_t size, const vector<double>& y) const
{
	auto points = y;
	process_points_vector(points);
	vector<vector<vector<complex<double>>>> result(nodes.size(),
		vector<vector<complex<double>>>(points.size(), vector<complex<double>>(size)));
	for (size_t i = 0; i < _roots.size(); i++)
	{
		const auto solution = smaller_system(_roots[i]).solve(points, { 0,1 });
		vector<complex<double>> multiplier(solution.size());
		for (size_t ii = 0; ii < solution.size(); ii++)
		{
			multiplier[ii] = I * solution[ii].front() / _residues[i];
		}
		const auto r = row(nodes, _roots[i]);
		for (size_t ii = 0; ii < nodes.size(); ii++)
		{
			const auto term = multiplier % r[ii];
			result[ii] = result[ii] + term;
		}
	}
	return result;
}

vector<complex<double>> layer::row(double node,
	complex<double> root, bool odd) const
{
	std::map<double, complex<double>> result;
	for (const auto& element : this->_elements)
	{
		const auto val = element->evaluate(root, node, odd);
		for (const auto& x : val)
		{
			result[x.first] += x.second;
		}
	}
	std::vector<complex<double>> v(result.size());
	std::transform(result.begin(), result.end(), v.begin(), [](auto x) {return x.second; });
	return v;
}

vector<vector<complex<double>>> layer::row(const vector<double>& node,
	complex<double> root) const
{
	vector<vector<complex<double>>> result(node.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = row(node[i], root);
	}
	return result;
}

vector<complex<double>> layer::get_roots() const
{
	return this->_roots;
}



complex<double> layer::right_part(double xi_1, double xi_2) const {
	complex<double> result = 0;
	complex<double> term;
	complex<double> alpha;
	size_t n = 0;
	do {
		alpha = _roots[n++];
		auto bigger = bigger_system(alpha);
		const auto solution = bigger.solve({ 0, xi_2, 1.0 }, { 0, 1, 0, 0 });
		term = solution[1][0] / solution.back().back()
			* exp(I * alpha * abs(xi_1));
		result += term;
	} while (abs(term) > _tolerance && n < _roots.size());
	return I * result;
}

vector<vector<complex<double>>> layer::matrix() const
{
	const auto nodes = get_nodal_points(this->_elements);
	const size_t matrix_size = nodes.size();
	vector<vector<complex<double>>> result(matrix_size);
	for (size_t i = 0; i < matrix_size; i++)
	{
		result[i] = row(nodes[i], matrix_size);
	}
	return result;
}

std::vector<double> layer::umoff_pointing(double x_1, double x_2,
	const std::vector<std::complex<double>>& solution) const
{
	const auto size = solution.size();
	vector<complex<double>> u(size);
	vector<complex<double>> sigma_13(size);
	vector<complex<double>> sigma_23(size);
	for (size_t i = 0; i < 1;/*_roots.size();*/ i++)
	{
		const auto multiplier = smaller_system(_roots[i]).solve(0, x_2, { {0.0, 0.0},{1.0, 0.0} })
			/ _residues[i];
		const auto r = row(x_1, _roots[i]);
		const auto r_odd = row(x_1, _roots[i], true);
		u = u + I * multiplier[0] * r;
		sigma_13 = sigma_13 + _roots[i] * multiplier[0] * r_odd;
		sigma_23 = sigma_23 + I * multiplier[1] * r;
	}
	const auto U = u * solution;
	const auto Sigma_13 = sigma_13 * solution; // !!!!
	const auto Sigma_23 = sigma_23 * solution;
	return { (I * U * Sigma_13).real(), (I * U * Sigma_23).real() };
}

std::complex<double> layer::displacement(double x_1, double x_2,
	const std::vector<std::complex<double>>& solution) const
{
	const auto r = this->row(x_1, solution.size(), x_2);
	return r * solution;
}

std::vector<std::complex<double>> layer::displacement(double x_1,
	const std::vector<double>& x_2,
	const std::vector<std::complex<double>>& solution) const
{
	const auto r = this->row(x_1, solution.size(), x_2);
	return r * solution;
}

std::vector<std::vector<std::complex<double>>> layer::displacement(
	const std::vector<double>& x_1, const std::vector<double>& x_2,
	const std::vector<std::complex<double>>& solution) const
{
	const auto r = this->row(x_1, solution.size(), x_2);
	std::vector<std::vector<std::complex<double>>> result(x_1.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = r[i] * solution;
	}
	return result;
}

void layer::process_points_vector(std::vector<double>& points) const
{
	if (abs(points.front()) > DBL_EPSILON)
	{
		points.insert(points.begin(), 0.0);
	}
	if (abs(points.back() - 1) > DBL_EPSILON)
	{
		points.push_back(1.0);
	}
}

void layer::save_to_file(const std::string& filename) const
{
	std::ofstream ofs(filename);
	ofs << setprecision(30);
	ofs << this->_kappa << "\n";
	for (size_t i = 0; i < this->_roots.size(); i++)
	{
		ofs << _roots[i].real() << " " << _roots[i].imag() << " ";
		ofs << _residues[i].real() << " " << _residues[i].imag() << "\n";
	}
	ofs.close();
}

void layer::load_from_file(const std::string& filename)
{
	std::ifstream ifs(filename);
	ifs >> this->_kappa;
	_roots.clear();
	_residues.clear();
	double re, im, re1, im1;
	while (ifs >> re >> im >> re1 >> im1)
	{
		_roots.emplace_back(re, im);
		_residues.emplace_back(re1, im1);
	}
	ifs.close();
}

inline layer::~layer() {
	for (auto& element : this->_elements)
		delete element;
}

layer::layer(const function<double(double)>& mu,
	const function<double(double)>& rho, double kappa, size_t roots_number,
	double tolerance) :_mu(mu), _rho(rho), _kappa(kappa), _tolerance(tolerance)
{
	roots(roots_number);
}

layer::layer(const function<double(double)>& mu,
	const function<double(double)>& rho, const std::string& filename,
	double tolerance) :_mu(mu), _rho(rho), _tolerance(tolerance)
{
	this->load_from_file(filename);
}

layer::layer(const function<double(double)>& mu,
	const function<double(double)>& rho, double kappa,
	const std::vector<boundary_element*>& elements, size_t roots_number,
	double tolerance) :layer(mu, rho, kappa, roots_number, tolerance)
{
	this->_elements = elements;
}

layer::layer(const function<double(double)>& mu,
	const function<double(double)>& rho, const std::string& filename,
	const std::vector<boundary_element*>& elements, double tolerance)
	:layer(mu, rho, filename, tolerance)
{
	this->_elements = elements;
}

std::vector<double> get_nodal_points(
	const std::vector<boundary_element*>& elements)
{
	std::set<double> result;
	for (const auto& element : elements)
	{
		const auto& points = element->get_nodal_points();
		result.insert(points.begin(), points.end());
	}
	return std::vector<double>(result.begin(), result.end());
}
