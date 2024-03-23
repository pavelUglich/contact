#include "constant_element.h"

const std::complex<double> I = { 0,1 };

std::map<double, std::complex<double>> constant_element::evaluate(
	std::complex<double> alpha, double xi, bool odd) const
{
	auto node = this->get_nodal_points().front();
	const auto sign = odd ? -1.0 : 1.0;
	if (xi <= _left)
	{
		const auto diff = exp(I * alpha * _left) - exp(I * alpha * _right);
		const auto result = exp(-I * alpha * xi) * diff * I / alpha;
		return { { node, result } };
	}
	if (xi >= _right)
	{
		const auto diff = exp(-I * alpha * _right) - exp(-I * alpha * _left) * sign;
		const auto result = exp(I * alpha * xi) * diff * I / alpha * sign;
		return { { node, result } };
	}
	const auto right = -I * (exp(I * alpha * (_right - xi)) - 1.0) / alpha;
	const auto left = -I * (exp(I * alpha * (xi - _left)) - 1.0) / alpha * sign;
	return { { node, left + right } };
}