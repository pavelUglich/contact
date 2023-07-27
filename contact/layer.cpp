//
//  layer.cpp
//  diploma
//
//  Created by Максим Целигоров on 26.04.2023.
//
#include "OdeSolver.h"
#include "layer.hpp"

const complex<double> I = { 0,1 };


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
			_roots[counter] = newton_method({ (a + b) / 2,0 });
			_residues[counter] = function_and_derivative(_roots[counter]).back();
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
			_roots[counter] = newton_method({ 0, (a + b) / 2 });
			_residues[counter] = function_and_derivative(_roots[counter]).back();
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
			[=](double x, const vector<complex<double>>& y) {return y[1] / _mu(x); },
			[=](double x, const vector<complex<double>>& y)
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

complex<double> layer::num1(OdeSolver<complex<double>>& smaller, double x2,
	double xi_2) const
{
	const auto sol1 = smaller.solve(0, 1, { 1, 0 });
	const auto sol2 = smaller.solve({ 0, x2, xi_2, 1 }, { 0, 1 });
	const auto& sigma_1_h = sol1.back();
	const auto& u_2_xi2 = sol2[2][0];
	const auto& u_2_x2 = sol2[1][0];
	const auto num = I * sigma_1_h * u_2_xi2 * u_2_x2;
	const auto& denum = sol2.back().back();
	return num;// / denum;
}

complex<double> layer::kernel_transformant(complex<double> alpha) const
{
	const auto sol = this->smaller_system(alpha).solve(0, 1, { 0,1 });
	return sol[0] / sol[1];
}

vector<complex<double>> layer::get_roots() const {
	return _roots;
}

complex<double> layer::fundamental_solution(double x1, double x2, double xi_1,
	double xi_2) const
{
	complex <double> result = 0;
	complex <double> term;
	size_t n = 0;
	do {
		const complex<double> alpha = _roots[n];
		auto smaller = smaller_system(alpha);
		const auto num = x2 < xi_2
			? num1(smaller, x2, xi_2)
			: num1(smaller, xi_2, x2);
		term = num * exp(I * alpha * abs(xi_1 - x1)) / _residues[n++];
		result += term;
	} while (abs(term) > _tolerance && n < _roots.size());
	return result;
}

complex<double> layer::right_part(double xi_1, double xi_2) const {
	complex <double> result = 0;
	complex <double> term;
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

complex<double> layer::b_constant(double uu) const
{
	complex<double> result = { 0,0 };
	for (size_t i = 0; i < _roots.size(); i++)
	{
		const complex<double> alpha = _roots[i];
		const auto u = smaller_system(alpha).solve(0, 1, { 0,1 });
		result += u.front() * exp(I * alpha * uu) / alpha / _residues[i];
	}
	return result;
}

layer::layer(const function<double(double)>& mu,
	const function<double(double)>& rho, double kappa, size_t roots_number,
	double tolerance) :_mu(mu), _rho(rho), _kappa(kappa), _tolerance(tolerance)
{
	roots(roots_number);
}
