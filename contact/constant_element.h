#pragma once
#include"boundary_element.h"

class constant_element : public boundary_element
{
public:
	constant_element(double left, double right) :boundary_element(left, right, { 0.5 * (left + right) }) {}
	// Унаследовано через boundary_element
	std::vector<std::complex<double>> evaluate(std::complex<double> alpha, double xi) const override;
};

