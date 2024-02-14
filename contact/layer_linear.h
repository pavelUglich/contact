#pragma once
#include"layer.hpp"

class layer_linear :
	public layer
{
	/// <summary>
	/// вспомогательная функция B для расчёта ядра
	/// </summary>
	/// <param name="uu">аргумент B</param>
	/// <returns></returns>
	complex<double> b_constant(double uu, double h) const;
	complex<double> b_diagonal(double h) const;
	complex<double> b_left(double uu, double h) const;

	std::complex<double> evaluate_exp(
		const std::vector<std::complex<double>>& tau, double a, double x_1, 
		std::complex<double> alpha) const override;

	std::complex<double> int_auxiliary(std::complex<double> tau_left, 
		std::complex<double> tau_right,	std::complex<double> alpha, double x_1, 
		double a, double b) const;

public:
	/// <summary>
	/// конструктор
	/// </summary>
	/// <param name="mu"> распределение модуля сдвига </param>
	/// <param name="rho"> распределение плотности </param>
	/// <param name="kappa"> безразмерная частота колебаний </param>
	/// <param name="roots_number"> максимальное количество корней </param>
	/// <param name="tolerance"> точность расчётов </param>
	layer_linear(const function<double(double)>& mu,
		const function<double(double)>& rho, double kappa,
		size_t roots_number = MAX_ROOTS_NUMBER, double tolerance = 0.1e-5);

	layer_linear(const function<double(double)>& mu,
		const function<double(double)>& rho, const std::string& filename,
		double tolerance = 0.1e-5);



	vector<vector<complex<double>>> matrix(double a, size_t size) const
		override;


};

