#include "linear_element.h"
#include "Integral.h"
#include <iostream>

const std::complex<double> I = { 0,1 };

std::map<double, std::complex<double>> linear_element::evaluate(
	std::complex<double> alpha, double xi,
	const std::function<double(double)>& weight, bool odd) const
{
	const auto& nodes = this->get_nodal_points();
	const auto sign = odd ? -1.0 : 1.0;
	const auto phi_left = [this](double x) {return (x - _right) / (_left - _right); };
	const auto phi_right = [this](double x) {return (x - _left) / (_right - _left); };
	if (xi <= _left || abs(xi - _left) < DBL_EPSILON)
	{
		const auto left = integral_gauss(_left, _right, [=](auto x) {return exp(I * alpha * (x - xi)) * phi_left(x) * weight(x); });
		const auto right = integral_gauss(_left, _right, [=](auto x) {return exp(I * alpha * (x - xi)) * phi_right(x) * weight(x); });
		return { {nodes[0],left}, {nodes[1],right} };
	}
	if (xi >= _right || abs(xi - _right) < DBL_EPSILON)
	{
		const auto left = integral_gauss(_left, _right, [=](auto x) {return exp(-I * alpha * (x - xi)) * phi_left(x) * weight(x); });
		const auto right = integral_gauss(_left, _right, [=](auto x) {return exp(-I * alpha * (x - xi)) * phi_right(x) * weight(x); });
		return { {nodes[0],left}, {nodes[1],right} };
	}
	const auto left0 = integral_gauss(_left, xi, [=](auto x) {return exp(-I * alpha * (x - xi)) * phi_left(x) * weight(x); });
	const auto right0 = integral_gauss(_left, xi, [=](auto x) {return exp(-I * alpha * (x - xi)) * phi_right(x) * weight(x); });
	const auto left1 = integral_gauss(xi, _right, [=](auto x) {return exp(I * alpha * (x - xi)) * phi_left(x) * weight(x); });
	const auto right1 = integral_gauss(xi, _right, [=](auto x) {return exp(I * alpha * (x - xi)) * phi_right(x) * weight(x); });
	return { {nodes[0],(left0 + right0) * sign}, {nodes[1],left1 + right1} };
}

std::complex<double> linear_element::evaluate_integral(
	const std::vector<std::complex<double>>& solution, 
	const std::function<double(double)>& weight) const
{
	const auto phi_left = [this](double x) {return (x - _right) / (_left - _right); };
	const auto phi_right = [this](double x) {return (x - _left) / (_right - _left); };
	const auto left = integral_gauss(_left, _right, [=](auto x) {return phi_left(x) * weight(x); });
	const auto right = integral_gauss(_left, _right, [=](auto x) {return phi_right(x) * weight(x); });
	return solution[0] * left + solution[1] * right;
}

