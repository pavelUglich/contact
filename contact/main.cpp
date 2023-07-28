#include "layer.hpp"



int main() {
	const auto mu = [](double) {return 1.0; };
	const auto rho = [](double) {return 1.0; };
	layer l = { mu, rho, 10.0 };
	auto matrix = SquareMatrix<complex<double>>(20, l.matrix(1.0, 20));
	auto sol = matrix.linSolve(vector<complex<double>>(20, { 1,0 }));
	system("pause");

}