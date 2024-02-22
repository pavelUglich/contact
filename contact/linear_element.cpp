#include "linear_element.h"

const std::complex<double> I = { 0,1 };

std::vector<std::complex<double>> linear_element::evaluate(
	std::complex<double> alpha, double xi) const
{
	const auto mm = 1.0 / alpha / alpha / (_left - _right);
	if (xi <= _left)
	{
		const auto left = mm * (exp(I * alpha * (_right - xi)) -
			exp(I * alpha * (_left - xi))
			+ I * alpha * (_left - _right) * exp(I * alpha * (_left - xi)));
		const auto right = -mm *
			(exp(I * alpha * (_right - xi)) - exp(I * alpha * (_left - xi))
				+ I * alpha * (_left - _right) * exp(I * alpha * (_right - xi)));
		return { left, right };
	}
	if (xi >= _right)
	{
		const auto left = mm *
			(exp(I * alpha * (xi - _right)) - exp(I * alpha * (xi - _left))
				- I * alpha * (_left - _right) * exp(I * alpha * (xi - _left)));
		const auto right = mm *
			(exp(I * alpha * (xi - _left)) - exp(I * alpha * (xi - _right))
				+ I * alpha * (_left - _right) * exp(I * alpha * (xi - _right)));
		return { left, right };
	}
	const auto left_0 = mm * (alpha * I * (xi - _right)
		- 1.0 + exp(I * alpha * (_right - xi)));
	const auto left_1 = -mm *
		(exp(I * alpha * (xi - _left)) + I * alpha * (_left - _right)
			* exp(I * alpha * (xi - _left)) - 1.0 + I * alpha * (_right - xi));
	const auto right_0 = -mm * (-1.0 + I * alpha
		* (xi - _left) + exp(I * alpha * (_right - xi)) + I * alpha * (_left - _right)
		* exp(I * alpha * (_right - xi)));
	const auto right_1 = -mm * (-exp(I * alpha * (xi - _left))
		+ 1.0 + I * alpha * (xi - _left));
	return { left_0 + right_0,left_1 + right_1 };
}

/*
std::complex<double> layer_linear::int_auxiliary(std::complex<double> tau_left,
	std::complex<double> tau_right, std::complex<double> alpha, double x_1,
	double a, double b) const
{
	if (x_1 <= a)
	{
		const auto left = 1.0 / alpha / alpha / (a - b) *
			(exp(I * alpha * (b - x_1)) - exp(I * alpha * (a - x_1))
				+ I * alpha * (a - b) * exp(I * alpha * (a - x_1)));
		const auto right = -1.0 / alpha / alpha / (a - b) *
			(exp(I * alpha * (b - x_1)) - exp(I * alpha * (a - x_1))
				+ I * alpha * (a - b) * exp(I * alpha * (b - x_1)));
		return tau_left * left + tau_right * right;
	}
	if (x_1 >= b)
	{
		const auto left = 1.0 / alpha / alpha / (a - b) *
			(exp(I * alpha * (x_1 - b)) - exp(I * alpha * (x_1 - a))
				- I * alpha * (a - b) * exp(I * alpha * (x_1 - a)));
		const auto right = 1.0 / alpha / alpha / (a - b) *
			(exp(I * alpha * (x_1 - a)) - exp(I * alpha * (x_1 - b))
				+ I * alpha * (a - b) * exp(I * alpha * (x_1 - b)));
		return tau_left * left + tau_right * right;
	}
	const auto left_0 = 1.0 / (a - b) / alpha / alpha * (alpha * I * (x_1 - b)
		- 1.0 + exp(I * alpha * (b - x_1)));
	const auto left_1 = -1.0 / (a - b) / alpha / alpha *
		(exp(I * alpha * (x_1 - a)) + I * alpha * (a - b)
			* exp(I * alpha * (x_1 - a)) - 1.0 + I * alpha * (b - x_1));
	const auto right_0 = -1.0 / (a - b) / alpha / alpha * (-1.0 + I * alpha
		* (x_1 - a) + exp(I * alpha * (b - x_1)) + I * alpha * (a - b)
			* exp(I * alpha * (b - x_1)));
	const auto right_1 = -1.0 / (a - b) / alpha / alpha * (-exp(I * alpha * (x_1 - a))
		+ 1.0 + I * alpha * (x_1 - a));
	return tau_left * (left_0 + right_0) + tau_right * (left_1 + right_1);
}

*/