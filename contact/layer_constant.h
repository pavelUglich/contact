#pragma once

#include"layer.hpp"
class layer_constant :
	public layer
{
	std::complex<double> evaluate_exp(
		const std::vector<std::complex<double>>& tau, double a, double x_1, 
		size_t size, std::complex<double> alpha) const;

	// int_auxiliary(tau[i], alpha, x_1, edges[i - 1], edges[i]);
	std::complex<double> int_auxiliary(std::complex<double> tau, 
		std::complex<double> alpha, double x_1, double a, double b) const;

public:
	/// <summary>
	/// вспомогательная функция B для расчёта ядра
	/// </summary>
	/// <param name="uu">аргумент B</param>
	/// <returns></returns>
	complex<double> b_constant(double uu) const;

	/// <summary>
	/// конструктор
	/// </summary>
	/// <param name="mu"> распределение модуля сдвига </param>
	/// <param name="rho"> распределение плотности </param>
	/// <param name="kappa"> безразмерная частота колебаний </param>
	/// <param name="roots_number"> максимальное количество корней </param>
	/// <param name="tolerance"> точность расчётов </param>
	layer_constant(const function<double(double)>& mu,
		const function<double(double)>& rho, double kappa,
		size_t roots_number = MAX_ROOTS_NUMBER, double tolerance = 0.1e-5);

	vector<vector<complex<double>>> matrix(double a, size_t size) const 
		override;

	std::complex<double> displacement(double x_1, double x_2,
		std::vector<std::complex<double>> solution, double a) const override;
};

