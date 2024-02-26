#include "boundary_element.h"
#include"vector_ops.h"

std::vector<double> boundary_element::get_nodal_points() const
{
	return _nodal_points;
}
