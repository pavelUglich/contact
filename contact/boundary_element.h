#pragma once
#include <vector>
#include <complex>
class boundary_element
{
	// концы отрезка
	double _left, _right;
	std::vector<double> _nodal_points;
public:
	boundary_element(double left, double right,
		const std::vector<double> nodal_points) :_left(left), _right(right), 
		_nodal_points(nodal_points) {}

	/// <summary>
	/// вспомогательная функция для отыскания интегралов
	/// </summary>
	/// <param name="alpha">показатель экспоненты</param>
	/// <param name="xi">горизонтальная координата</param>
	/// <returns></returns>
	std::complex<double> evaluate(std::complex<double> alpha, double xi) const;

};

