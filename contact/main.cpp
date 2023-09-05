#include "layer.hpp"
#include <numeric>
#include"plots.h"
#include "layer_constant.h"
#include "layer_linear.h"
#include "Integral.h"

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
	const auto _mu = [](double x) {return 1.0; };
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
	const auto solution = os.solve(0, 1, { 0,exp(-abs(alpha))});
	return solution.front() / solution.back();
}


int main() {
	const auto mu = [](double x) {return 1.0; };
	const auto rho = [](double) {return 1.0; };
	const auto mass = 1.0;// масса штампа
	const auto kappa = 1.0; // волновое число
	const auto tt = 1.0; // интенсивность нагрузки	
	//layer_constant l = { mu, rho, kappa };
	//std::vector<complex<double>> field(20);
	double x = 0;
	double h = 1.0 / 20;

	/*
	for (size_t i = 0; i < 20; i++)
	{
		field[i] = l.b_constant(x);
		x += h;
	}
	plotTheWaveField(h, { {"red", field} }, "1.txt");
	*/
	Integral ii([](auto x) {return transformant(x, 1.0); }, 1.0);
	x = h;
	//double h = 1.0 / 20;
	std::vector<complex<double>> ff(20);
	for (size_t i = 0; i < 20; i++)
	{
		ff[i] = ii.evaluate(x);
		x += h;
	}
	plotTheWaveField(h, { {"red", ff} }, "2.txt");
	/*
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
	}*/
	system("pause");
}