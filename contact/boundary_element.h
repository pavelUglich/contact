#pragma once
#include <vector>
#include <map>
#include <complex>

/// <summary>
/// виртуальный класс
/// </summary>
class boundary_element
{
	// концы отрезка
protected:
	double _left;
	double _right;
	std::vector<double> _nodal_points;
public:
	boundary_element(double left, double right,
		const std::vector<double>& nodal_points) :_left(left), _right(right), 
		_nodal_points(nodal_points) {}

	/// <summary>
	/// вспомогательная функция для отыскания интегралов
	/// </summary>
	/// <param name="alpha">показатель экспоненты</param>
	/// <param name="xi">горизонтальная координата</param>
	/// <returns></returns>
	virtual std::map<double, std::complex<double>> evaluate(std::complex<double> alpha,
		double xi, bool odd = false) const = 0;
	std::vector<double> get_nodal_points() const;
	virtual ~boundary_element() = default;
};

