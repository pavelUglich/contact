#include "layer.hpp"
#include <numeric>
#include"plots.h"
#include "layer_constant.h"
#include "layer_linear.h"
#include "Integral.h"

const double pi = acos(-1.0);
enum class element_type { constant, linear };

layer* make_layer(element_type et, const function<double(double)>& mu,
	const function<double(double)>& rho, double kappa,
	size_t roots_number = MAX_ROOTS_NUMBER, double tolerance = 0.1e-5)
{
	switch (et)
	{
	case element_type::constant:
		return new layer_constant(mu, rho, kappa, roots_number, tolerance);
	case element_type::linear:
		return new layer_linear(mu, rho, kappa, roots_number, tolerance);
	default:
		break;
	}
	throw std::exception();
}

complex<double> transformant(complex<double> alpha, double _kappa)
{
	//return 1.0 / alpha;
	/*const auto _mu = [](double) {return 1.0; };
	const auto _rho = [](double) {return 1.0; }; 
	OdeSolver<complex<double>> os =
	{
		{
			[=](double x, const vector<complex<double>>& y) {return y[1] / _mu(x); },
			[=](double x, const vector<complex<double>>& y)
			{
				return (alpha * alpha * _mu(x) - _rho(x) * _kappa * _kappa) * y[0];
			}
		},
		0.1e-6,
		RUNGE_KUTTA_FELDBERG
	};
	const auto solution = os.solve(0, 1, { 0,exp(-abs(alpha)) });
	return solution.front() / solution.back() / alpha;//*/
	const auto gamma = sqrt(alpha * alpha - _kappa * _kappa);
	return sinh(gamma) / gamma / cosh(gamma);
}

complex<double> analytical(double x, double _kappa, double eps = 0.1e-5)
{
	complex<double> result = 0;
	size_t n = 0;
	complex<double> term;
	do
	{
		const complex<double> frac = (2 * n++ + 1) / 2.0 * pi;
		const complex<double> alpha = sqrt(_kappa * _kappa - (frac * frac));
		term = exp(I * alpha * x) / alpha;
		result += term / alpha;
	} while (abs(term)>eps);
	return result;
}



int main() {
	const auto mu = [](double x) {return 1.0; };
	const auto rho = [](double) {return 1.0; };
	const auto mass = 1.0;// масса штампа
	const auto kappa = 1.0; // волновое число
	const auto tt = 1.0; // интенсивность нагрузки	
	//layer_constant l = { mu, rho, kappa };
	//std::vector<complex<double>> field(20);
	double h = 1.0 / 20;
	double x = h / 2;
/*
	Integral ii([=](auto x) {return transformant(x, kappa); }, kappa);
	layer_constant l = { mu, rho, kappa };
	std::vector<complex<double>> ff(20);
	std::vector<complex<double>> ff1(20);
	for (size_t i = 0; i < 20; i++)
	{
		cout << x << "\n";
		ff1[i] = analytical(x, kappa);//I*l.b_constant(x);
		ff[i] = l.b_constant(x); //ii.evaluate(x);
		x += h;
	}

	plotTheWaveField(h, { {"red", ff}, {"black", ff1}}, "2.txt", h / 2);
	*/
	
	{
		layer_constant l = { mu, rho, kappa };
		size_t size = 20;
		const double a = 1.0;
		const auto h = 2 * a / (double)size;
		const auto mm = l.matrix(a, size);
		auto matrix = SquareMatrix<complex<double>>(size, mm);
		auto sol = matrix.linSolve(vector<complex<double>>(size, { 1,0 }));
		//const complex<double> integral_ = h * std::accumulate(sol.begin(), sol.end(), complex<double>());
		//const auto gamma = (integral_ - tt) / mass / kappa / kappa;
		//std::transform(sol.begin(), sol.end(), sol.begin(), [=](auto x) {return x / gamma; });
		plotTheWaveField(h, { {"red", sol} }, "1.txt", -a  + h / 2); 
	}
	
	{
		layer_linear l = { mu, rho, kappa };
		size_t size = 20;
		const double a = 1.0;
		const auto h = 2 * a / (double)size;
		const auto mm = l.matrix(a, size);
		auto matrix = SquareMatrix<complex<double>>(size + 1, mm);
		auto sol = matrix.linSolve(vector<complex<double>>(size + 1, { 1,0 }));
		//const complex<double> integral_ = h * std::accumulate(sol.begin(), sol.end(), complex<double>());
		//const auto gamma = (integral_ - tt) / mass / kappa / kappa;
		//std::transform(sol.begin(), sol.end(), sol.begin(), [=](auto x) {return x / gamma; });
		plotTheWaveField(h, { {"black", sol} }, "2.txt", -a);
	}//*/
	system("pause");
}