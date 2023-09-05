#include "Integral.h"
#include<functional>
using namespace std;
const double pi = 3.1415926538;

complex<double> Integral::gauss(double a, double b, double x,
	const function<complex<double>(double)>& contour,
	const function<complex<double>(double)>& derivative) const
{
	const double m = 0.5 * (a + b);
	const double n = 0.5 * (b - a);
	complex<double> result;
	for (size_t i = 0; i < 16; i++)
	{
		double y = m + n * xi[i];
		complex<double> alpha = contour(y);
		auto symbol = _integrand(alpha);
		symbol *= sin(alpha * x);
		complex<double> da = derivative(y);
		symbol *= da;
		result += eta[i] * symbol;
	}
	return result;
}

complex<double> Integral::gauss(double a, double b, double x) const
{
	const double m = 0.5 * (a + b);
	const double n = 0.5 * (b - a);
	complex<double> result;
	for (size_t i = 0; i < 16; i++)
	{
		double y = m + n * xi[i];
		auto symbol = _integrand(y);
		symbol *= sin(y * x);
		result += eta[i] * symbol;
	}
	return result;
}


complex<double> Integral::evaluate(double x) const {
	complex<double> result(0, 0);
	double delta = _kappa / 10;
	for (size_t j = 0; j < 8; j++) {
		result += gauss(0.25 * _kappa * j, 0.25 * _kappa * (j + 1), x,
			[delta, this](double y) {return complex<double>(y, -1 * delta * (1 - (y - _kappa) * (y - _kappa) / (_kappa * _kappa))); }, 
			[delta, this](double y) {return complex<double>(1, 2 * delta * (y - _kappa) / (_kappa * _kappa)); });
	}
	result *= 0.125 * _kappa;
	double segment = 2 * pi / x;
	double x0 = 2 * _kappa;
	complex<double> term;
	do
	{
		term = gauss(x0, x0 + segment, x);
		result += term * 0.5 * segment;
		x0 += segment;
	} while (abs(term) > 0.1e-3);//*/
	return result / pi;
}

Integral::~Integral() = default;
