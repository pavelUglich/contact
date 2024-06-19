#pragma once
#include"boundary_element.h"
#include <map>

class constant_element : public boundary_element
{
	
public:
	constant_element(double left, double right) :boundary_element(left, right, { 0.5 * (left + right) }) {}
	//std::map<double, std::complex<double>> evaluate(std::complex<double> alpha, 
	//	double xi, bool odd) const override;
	std::map<double, std::complex<double>> evaluate(std::complex<double> alpha, 
		double xi, const std::function<double(double)>& weight, bool odd) const override;
	/*
		virtual std::map<double, std::complex<double>> evaluate(std::complex<double> alpha,
		double xi, const std::function<double(double)>& weight, 
		bool odd = false) const = 0;

	*/
};

