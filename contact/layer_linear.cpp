#include "layer_linear.h"

layer_linear::layer_linear(const function<double(double)>& mu,
	const function<double(double)>& rho, double kappa, size_t roots_number,
	double tolerance) :layer(mu, rho, kappa, roots_number, tolerance)
{
}
