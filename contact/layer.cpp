#include "OdeSolver.h"
#include "layer.hpp"
#include <numeric>
#include <fstream>
#include <iomanip>


vector<complex<double>> layer::function_and_derivative(complex<double> alpha) const {
	const auto sol1 = bigger_system(alpha).solve(0, 1, { 0,1,0,0 });//!!!
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

std::complex<double> layer::displacement(double x_1, double x_2, 
	const std::vector<std::complex<double>>& solution, double a) const
{
	const auto roots = get_roots();
	const auto residues = get_residues();
	std::complex<double> result;
	for (size_t i = 0; i < roots.size(); i++)
	{
		const auto numerator = smaller_system(roots[i]).solve(0, x_2, { 0,1 }).front();
		result += numerator / residues[i] * evaluate_exp(solution, a, x_1, roots[i]);
	}
	return I * result;
}

std::complex<double> layer::displacement_auxiliary(double x_2,
	std::complex<double> alpha) const
{
	std::vector<double> points = { x_2 };
	process_points_vector(points);
	const auto solution = this->bigger_system(alpha).solve(points, { 0,1,0,0 });
	return solution[1][0] / solution.back().back();
}

std::vector<std::complex<double>> layer::displacement_auxiliary(
	const std::vector<double>& x_2, std::complex<double> alpha) const
{
	auto points = x_2;
	process_points_vector(points);
	auto solution = this->bigger_system(alpha).solve(points, { 0,1,0,0 });
	const auto denumerator = solution.back().back();
	std::vector<std::complex<double>> result(solution.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i] = solution[i][0] / denumerator;
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
	for (size_t i = 0; i < this->_roots.size(); i++)
	{
		ofs << setprecision(30);
		ofs << this->_kappa << "\n";
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
	while (!ifs.eof())
	{
		double re, im;
		ifs >> re >> im;
		_roots.emplace_back(re, im);
		ifs >> re >> im;
		_residues.emplace_back(re, im);
	}
	ifs.close();
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
