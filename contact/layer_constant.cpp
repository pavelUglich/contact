#include "layer_constant.h"

std::complex<double> layer_constant::evaluate_exp(
	const std::vector<std::complex<double>>& tau, double a, double x_1,
	std::complex<double> alpha) const
{
	const auto size = tau.size();
	std::vector<double> nodes(size); // срединные точки отрезков
	std::vector<double> edges(size + 1); // концы отрезков
	const auto h = 2 * a / size;
	edges[0] = -a;
	for (size_t i = 0; i < size; i++)
	{
		nodes[i] = (2 * i + 1) * h / 2 - a;
		edges[i + 1] = (i + 1) * h - a;
	}
	std::complex<double> result = 0;
	for (size_t i = 0; i < size; i++)
	{
		result += int_auxiliary(tau[i], alpha, x_1, edges[i], edges[i + 1]);
	}
	return result;
}

std::complex<double> layer_constant::int_auxiliary(std::complex<double> tau,
	std::complex<double> alpha, double x_1, double a, double b) const
{
	if (x_1 <= a)
	{
		const auto diff = exp(I * alpha * a) - exp(I * alpha * b);
		return tau * exp(-I * alpha * x_1) * diff * I / alpha;
	}
	if (x_1 >= b)
	{
		const auto diff = exp(-I * alpha * b) - exp(-I * alpha * a);
		return tau * exp(I * alpha * x_1) * diff * I / alpha;
	}
	const auto left = int_auxiliary(tau, alpha, x_1, a, x_1);
	const auto right = int_auxiliary(tau, alpha, x_1, x_1, b);
	return tau * (left + right);
}

complex<double> layer_constant::b_constant(double uu) const
{
	complex<double> result = { 0,0 };
	const auto _roots = this->get_roots();
	const auto _residues = this->get_residues();
	const auto _tolerance = this->get_tolerance();
	for (size_t i = 0; i < _roots.size(); i++)
	{
		const complex<double> alpha = _roots[i];
		const auto u = smaller_system(alpha).solve(0, 1, { 0,1 });
		const auto term = u.front() * exp(I * alpha * uu) / _residues[i] / alpha;
		result += term;
		if (abs(term) < _tolerance)
		{
			break;
		}
	}
	return result;
}

layer_constant::layer_constant(const function<double(double)>& mu,
	const function<double(double)>& rho, double kappa, size_t roots_number,
	double tolerance) :layer(mu, rho, kappa, roots_number, tolerance)
{
}

layer_constant::layer_constant(const function<double(double)>& mu,
	const function<double(double)>& rho, const std::string& filename,
	double tolerance) :layer(mu, rho, filename, tolerance)
{
}

vector<vector<complex<double>>>
layer_constant::matrix(double a, size_t size) const
{
	const double h = 2 * a / (double)size;
	vector<complex<double>> bb(size + 1);
	bb[0] = b_constant(0);
	double x = h / 2; // !!!
	for (size_t i = 0; i < size; i++)
	{
		bb[i + 1] = b_constant(x);
		x += h;
	}
	vector<vector<complex<double>>> result(size, vector<complex<double>>(size));
	for (size_t i = 0; i < size; i++)
	{
		for (size_t ii = i; ii < size; ii++)
		{
			result[i][ii] = bb[ii - i + 1] - bb[ii - i];
			result[ii][i] = result[i][ii];
		}
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i][i] *= 2;
	}
	return result;
}

