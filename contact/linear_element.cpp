#include "linear_element.h"
#include "Integral.h"
#include <iostream>

const std::complex<double> I = { 0,1 };

/*
std::map<double, std::complex<double>> linear_element::evaluate(
	std::complex<double> alpha, double xi, bool odd) const
{
	//xi = _left + (_right - _left) * rand() / RAND_MAX;
	const auto& nodes = this->get_nodal_points();
	//const auto mm = 1.0 / alpha / alpha / (_left - _right);
	const auto sign = odd ? -1.0 : 1.0;
	const auto phi_left = [this](double x) {return (x - _right) / (_left - _right); };
	const auto phi_right = [this](double x) {return (x - _left) / (_right - _left); };
	if (xi <= _left || abs(xi - _left) < DBL_EPSILON)
	{
		const auto left = integral_gauss(_left, _right, [=](auto x) {return exp(I * alpha * (x - xi)) * phi_left(x); });
		const auto right = integral_gauss(_left, _right, [=](auto x) {return exp(I * alpha * (x - xi)) * phi_right(x); });
		/*		const auto er = exp(I * alpha * (_right - xi));
				const auto el = exp(I * alpha * (_left - xi));
				const auto left = mm * (er - el + I * alpha * (_left - _right) * el);
				const auto right = mm * (el - er - I * alpha * (_left - _right) * er);//*/
/*		return {{nodes[0],left}, {nodes[1],right}};
	}
	if (xi >= _right || abs(xi - _right) < DBL_EPSILON)
	{
		const auto left = integral_gauss(_left, _right, [=](auto x) {return exp(-I * alpha * (x - xi)) * phi_left(x); });
		const auto right = integral_gauss(_left, _right, [=](auto x) {return exp(-I * alpha * (x - xi)) * phi_right(x); });
		/*const auto er = exp(I * alpha * (xi - _right)) * sign;
		const auto el = exp(I * alpha * (xi - _left)) * sign;
		const auto left = mm * (er - el - I * alpha * (_left - _right) * el);
		const auto right = mm * (el - er + I * alpha * (_left - _right) * er);//*/
/*		return {{nodes[0],left}, {nodes[1],right}};
	}
	const auto left0 = integral_gauss(_left, xi, [=](auto x) {return exp(-I * alpha * (x - xi)) * phi_left(x); });
	const auto right0 = integral_gauss(_left, xi, [=](auto x) {return exp(-I * alpha * (x - xi)) * phi_right(x); });
	const auto left1 = integral_gauss(xi, _right, [=](auto x) {return exp(I * alpha * (x - xi)) * phi_left(x); });
	const auto right1 = integral_gauss(xi, _right, [=](auto x) {return exp(I * alpha * (x - xi)) * phi_right(x); });
	/*const auto erx = exp(I * alpha * (xi - _left));
	const auto exl = exp(I * alpha * (_right - xi));
	const auto left_0 = mm * (alpha * I * (xi - _right) - 1.0 + erx);
	const auto left_1 = -mm * (exl + I * alpha * (_left - _right) * exl - 1.0 +
		I * alpha * (_right - xi));
	const auto right_0 = -mm * (-1.0 + I * alpha * (xi - _left) + erx +
		I * alpha * (_left - _right) * erx);
	const auto right_1 = -mm * (-exl + 1.0 + I * alpha * (xi - _left));
	std::cout << left1 + right1 << " " << left_1 + right_1 << std::endl;*/
	//	return { {nodes[0],(left_0 + right_0) * sign}, {nodes[1],left_1 + right_1} };
/*	return {{nodes[0],(left0 + right0) * sign}, {nodes[1],left1 + right1}};
}
*/


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

