#include "boundary_element.h"
#include"vector_ops.h"

std::complex<double> boundary_element::calculate(std::complex<double> alpha, 
	double xi, const std::vector<std::complex<double>>& nodal_values) const
{
	if (this->_nodal_points.size() != nodal_values.size())
	{
		throw std::out_of_range("");
	}
    return nodal_values * evaluate(alpha, xi);
}
