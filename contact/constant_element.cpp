#include "constant_element.h"
#include "Integral.h"
#include <iostream>

const std::complex<double> I = { 0,1 };

std::map<double, std::complex<double>> constant_element::evaluate(std::complex<double> alpha,
	double xi, const std::function<double(double)>& weight, bool odd) const
{
	auto node = this->get_nodal_points().front();
	if (xi <= _left)
	{
		const auto result = integral_gauss(_left, _right, [=](auto x) {
			return exp(I * alpha * (x - xi)) * weight(x); });
		return { { node, result } };
	}
	if (xi >= _right)
	{
		const auto result = integral_gauss(_left, _right, [=](auto x) {
			return exp(-I * alpha * (x - xi)) * weight(x); });
		return { { node, result } };
	}
	const auto left = integral_gauss(_left, xi, [=](auto x) {
		return exp(-I * alpha * (x - xi)) * weight(x); });
	const auto right = integral_gauss(xi, _right, [=](auto x) {
		return exp(I * alpha * (x - xi)) * weight(x); });
	return { { node, left + right } };
}

std::complex<double> constant_element::evaluate_integral(
	const std::vector<std::complex<double>>& solution, 
	const std::function<double(double)>& weight) const
{	
	return solution.front() * integral_gauss(this->_left, this->_right, weight);
}
