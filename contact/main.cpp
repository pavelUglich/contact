#include "layer.hpp"



int main() {
	const auto mu = [](double) {return 1.0; };
	const auto rho = [](double) {return 1.0; };
	layer l = { mu, rho, 1.0, 0 };

}