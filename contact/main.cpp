#include "layer.hpp"
#include <numeric>
#include"plots.h"


int main() {
	const auto mu = [](double x) {return 1.0 + x; };
	const auto rho = [](double) {return 1.0; };
	const auto mass = 1.0;// масса штампа
	const auto kappa = 2.0; // волновое число
	const auto T = 1.0; // интенсивность нагрузки
	layer l = { mu, rho, kappa };
	size_t size = 20;
	const double a = 1.0;
	const auto h = 2 * a / size;
	auto matrix = SquareMatrix<complex<double>>(size, l.matrix(a, size));
	auto sol = matrix.linSolve(vector<complex<double>>(size, { 1,0 }));
	const complex<double> integral_ = h * std::accumulate(sol.begin(), sol.end(), complex<double>());
	const auto gamma = (integral_ - T) / mass / kappa / kappa;
	std::transform(sol.begin(), sol.end(), sol.begin(), [=](auto x) {return x / gamma; });
	plotTheWaveField(h, { {"red", sol} }, "1.txt", -a + h / 2);
	system("pause");

}