#include "linear_element.h"

const std::complex<double> I = { 0,1 };

std::vector<std::complex<double>> linear_element::evaluate(
	std::complex<double> alpha, double xi) const
{
	const auto mm = 1.0 / alpha / alpha / (_left - _right);
	if (xi <= _left)
	{
		const auto er = exp(I * alpha * (_right - xi));
		const auto el = exp(I * alpha * (_left - xi));
		const auto left = mm * (er - el + I * alpha * (_left - _right) * el);
		const auto right = -mm * (er - el + I * alpha * (_left - _right) * er);
		return { left, right };
	}
	if (xi >= _right)
	{
		const auto er = exp(I * alpha * (xi - _right));
		const auto el = exp(I * alpha * (xi - _left));
		const auto left = mm * (er - el - I * alpha * (_left - _right) * el);
		const auto right = mm * (el - er + I * alpha * (_left - _right) * er);
		return { left, right };
	}
	const auto erx = exp(I * alpha * (_right - xi));
	const auto exl = exp(I * alpha * (_right - xi));
	const auto left_0 = mm * (alpha * I * (xi - _right) - 1.0 + erx);
	const auto left_1 = -mm * (exl + I * alpha * (_left - _right) * exl - 1.0 +
		I * alpha * (_right - xi));
	const auto right_0 = -mm * (-1.0 + I * alpha * (xi - _left) + erx +
		I * alpha * (_left - _right) * erx);
	const auto right_1 = -mm * (-exl + 1.0 + I * alpha * (xi - _left));
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