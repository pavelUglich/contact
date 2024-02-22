#pragma once
#include "boundary_element.h"
class linear_element :
    public boundary_element
{
public:

    linear_element(double left, double right) :boundary_element(left, right, { left, right }) {}

    // Унаследовано через boundary_element
    std::vector<std::complex<double>> evaluate(std::complex<double> alpha, double xi) const override;
};

