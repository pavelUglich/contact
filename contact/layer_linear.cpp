#include "layer_linear.h"

complex<double> layer_linear::b_constant(double uu, double h) const
{
	complex<double> result = { 0,0 };
	const auto _roots = this->get_roots();
	const auto _residues = this->get_residues();
	const auto _tolerance = this->get_tolerance();
	for (size_t i = 0; i < _roots.size(); i++)
	{
		const complex<double> alpha = _roots[i];
		const auto u = smaller_system(alpha).solve(0, 1, { 0,1 }); // !!!
		const auto term = 2.0 * I * u.front() * exp(I * alpha * uu) / alpha
			/ alpha * (1.0 - cos(alpha * h)) / _residues[i] / h;
		result += term;
		if (abs(term) < _tolerance)
		{
			break;
		}
	}
	return result;
}

complex<double> layer_linear::b_diagonal(double h) const
{
	complex<double> result = { 0,0 };
	const auto _roots = this->get_roots();
	const auto _residues = this->get_residues();
	const auto _tolerance = this->get_tolerance();
	for (size_t i = 0; i < _roots.size(); i++)
	{
		const complex<double> alpha = _roots[i];
		const auto u = smaller_system(alpha).solve(0, 1, { 0,1 }); // !!!
		const auto term = 2.0 * I * u.front() / alpha / alpha
			* (1.0 + I * alpha * h - exp(I * alpha * h)) / _residues[i] / h;
		result += term;
		if (abs(term) < _tolerance)
		{
			break;
		}
	}
	return result;
}

complex<double> layer_linear::b_left(double uu, double h) const
{
	complex<double> result = { 0,0 };
	const auto _roots = this->get_roots();
	const auto _residues = this->get_residues();
	const auto _tolerance = this->get_tolerance();
	for (size_t i = 0; i < _roots.size(); i++)
	{
		const complex<double> alpha = _roots[i];
		const auto u = smaller_system(alpha).solve(0, 1, { 0,1 }); // !!!
		const auto term = I * u.front() * exp(I * alpha * abs(uu)) / alpha
			/ alpha * (1.0 - exp(I * alpha * h) + I * alpha * h) / _residues[i] / h;
		result += term;
		if (abs(term) < _tolerance)
		{
			break;
		}
	}
	return result;
}

std::complex<double> layer_linear::evaluate_exp(
	const std::vector<std::complex<double>>& tau, double a, double x_1,
	std::complex<double> alpha) const
{
	const auto size = tau.size();
	std::vector<double> nodes(size + 1); // срединные точки отрезков
	std::vector<double> edges(size + 1); // концы отрезков
	const auto h = 2 * a / size;
	for (size_t i = 0; i <= size; i++)
	{
		nodes[i] = i * h;
		edges[i] = i * h;
	}
	std::complex<double> result = 0;
	for (size_t i = 0; i < size; i++)
	{
		result += int_auxiliary(tau[i], tau[i + 1], alpha, x_1, edges[i - 1],
			edges[i]);
	}
	return result;
}

std::complex<double> layer_linear::int_auxiliary(std::complex<double> tau_left,
	std::complex<double> tau_right, std::complex<double> alpha, double x_1,
	double a, double b) const
{
	if (x_1 < a)
	{
		const auto left = 1.0 / alpha / alpha / (a - b) *
			(exp(I * alpha * (b - x_1)) - exp(I * alpha * (a - x_1))
				+ I * alpha * (a - b) * exp(I * alpha * (a - x_1)));
		const auto right = 1.0 / alpha / alpha / (a - b) *
			(exp(I * alpha * (b - x_1)) - exp(I * alpha * (a - x_1))
				- I * alpha * (a - b) * exp(I * alpha * (b - x_1)));
		return tau_left * left + tau_right * right;
	}
	if (x_1 > b)
	{
		const auto left = 1.0 / alpha / alpha / (a - b) *
			(exp(I * alpha * (x_1 - b)) - exp(I * alpha * (x_1 - a))
				- I * alpha * (a - b) * exp(I * alpha * (x_1 - a)));
		const auto right = 1.0 / alpha / alpha / (a - b) *
			(exp(I * alpha * (x_1 - a)) - exp(I * alpha * (x_1 - b))
				+ I * alpha * (a - b) * exp(I * alpha * (x_1 - b)));
		return tau_left * left + tau_right * right;
	}
	const auto left_0 = 1.0 / (a - b) / alpha / alpha * (alpha * I * (a - b)
		* exp(I * alpha * (a - x_1)) - exp(I * alpha * (a - x_1))
		+ exp(I * alpha * (a - x_1)));
	const auto left_1 = 1.0 / (a - b) / alpha / alpha * (-alpha * I * (a - b)
		* exp(I * alpha * (b - x_1)) + exp(I * alpha * (a - x_1))
		+ exp(I * alpha * (b - x_1)));
	const auto right_0 = - 1.0 / (a - b) / alpha / alpha * (alpha * I * (a - b)
		* exp(I * alpha * (x_1 - a)) + exp(I * alpha * (x_1 - a))
		- exp(I * alpha * (x_1 - b)));
	const auto right_1 = - 1.0 / (a - b) / alpha / alpha * (-alpha * I * (a - b)
		* exp(I * alpha * (x_1 - b)) + exp(I * alpha * (x_1 - b))
		- exp(I * alpha * (x_1 - a)));
	return tau_left * (left_0 + right_0) + tau_right * (left_1 + right_1);
}

layer_linear::layer_linear(const function<double(double)>& mu,
	const function<double(double)>& rho, double kappa, size_t roots_number,
	double tolerance) :layer(mu, rho, kappa, roots_number, tolerance)
{
}

vector<vector<complex<double>>> layer_linear::matrix(double a, size_t size) const
{
	const double h = 2 * a / size;
	vector<vector<complex<double>>> result(size + 1, vector<complex<double>>(size + 1));
	// 1. çàïîëíÿåì ëåâûé è ïðàâûé ñòîëáöû
	vector<complex<double>> bb(size + 1);
	double x = 0; // !!!
	for (size_t i = 0; i < bb.size(); i++)
	{
		bb[i] = b_left(x, h);
		x += h;
	}
	// 2. çàïîëíÿåì äèàãîíàëü
	const auto diag = b_diagonal(h);
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i][i] = diag;
	}
	vector<complex<double>> b(size);
	x = h;
	for (size_t i = 0; i < b.size(); i++)
	{
		b[i] = b_constant(x, h);
		x += h;
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t ii = i + 1; ii < result.size(); ii++)
		{
			result[i][ii] = b[ii - i - 1];
			result[ii][i] = b[ii - i - 1];
		}
	}
	for (size_t i = 0; i < result.size(); i++)
	{
		result[i][0] = bb[i];
		result[i].back() = bb[bb.size() - i - 1];
	}
	return result;
}

