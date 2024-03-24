#include "linear_element.h"

const std::complex<double> I = { 0,1 };

std::map<double, std::complex<double>> linear_element::evaluate(
	std::complex<double> alpha, double xi, bool odd) const
{
	const auto& nodes = this->get_nodal_points();
	const auto mm = 1.0 / alpha / alpha / (_left - _right);
	const auto sign = odd ? -1.0 : 1.0;
	if (xi <= _left || abs(xi - _left) < DBL_EPSILON)
	{
		const auto er = exp(I * alpha * (_right - xi));
		const auto el = exp(I * alpha * (_left - xi));
		const auto left = mm * (er - el + I * alpha * (_left - _right) * el);
		const auto right = mm * (el - er - I * alpha * (_left - _right) * er);
		return { {nodes[0],left}, {nodes[1],right} };
	}
	if (xi >= _right || abs(xi - _right) < DBL_EPSILON)
	{
		const auto er = exp(I * alpha * (xi - _right)) * sign;
		const auto el = exp(I * alpha * (xi - _left)) * sign;
		const auto left = mm * (er - el - I * alpha * (_left - _right) * el);
		const auto right = mm * (el - er + I * alpha * (_left - _right) * er);
		return { {nodes[0],left}, {nodes[1],right} };
	}
	const auto erx = exp(I * alpha * (xi - _left));
	const auto exl = exp(I * alpha * (_right - xi));
	const auto left_0 = mm * (alpha * I * (xi - _right) - 1.0 + erx);
	const auto left_1 = -mm * (exl + I * alpha * (_left - _right) * exl - 1.0 +
		I * alpha * (_right - xi));
	const auto right_0 = -mm * (-1.0 + I * alpha * (xi - _left) + erx +
		I * alpha * (_left - _right) * erx);
	const auto right_1 = -mm * (-exl + 1.0 + I * alpha * (xi - _left));
	return { {nodes[0],(left_0 + right_0) * sign}, {nodes[1],left_1 + right_1} };
}

