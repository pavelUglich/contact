//
//  layer.hpp
//  diploma
//
//  Created by Максим Целигоров on 26.04.2023.
//

#ifndef layer_hpp
#define layer_hpp
#include<vector>
#include <complex>
#include <stdio.h>
#include <functional>
#include"OdeSolver.h"

using namespace std;
const size_t MAX_ROOTS_NUMBER = 100;
const complex<double> I = { 0, 1 };

class layer {
	// распределения модуля сдвига и плотности
	function<double(double)> _mu, _rho;

	// частота колебаний (безразмерная)
	double _kappa;

	// количество корней дисперсионного уравнения
	vector<complex<double>> _roots;
	vector<complex<double>> _residues;



	// точность вычислений
	double _tolerance;

	/// <summary>
	/// вспомогательная функция для отыскания значения трансформанты напряжения
	/// вместе с её производной по параметру \alpha
	/// </summary>
	/// <param name="alpha"> параметр преобразования Фурье</param>
	/// <returns> два комплексных числа</returns>
	vector<complex<double>> function_and_derivative(complex<double> alpha) const;

	/// <summary>
	/// метод Ньютона
	/// </summary>
	/// <param name="beginParam"></param>
	/// <returns></returns>
	complex<double> newton_method(complex<double> begin_param) const;

	/// <summary>
	/// вещественная часть трансформанты напряжения на верхней поверхности слоя
	///  (используется для решения дисперсионного уравнения)
	/// </summary>
	/// <param name="alpha"> параметр преобразования Фурье </param>
	/// <returns> вещественное число </returns>
	double equation(complex<double> alpha) const;

	/// <summary>
	/// заполнение вектора корней
	/// </summary>
	/// <param name="roots_number">количество корней</param>
	void roots(size_t roots_number);

	/// <summary>
	/// система ОДУ в трансформантах (используем, если производная по \alpha не нужна)
	/// </summary>
	/// <param name="alpha">  параметр преобразования Фурье  </param>
	/// <returns></returns>
	OdeSolver<complex<double>> smaller_system(complex<double> alpha) const;

	/// <summary>
	/// система ОДУ в трансформантах (используем, если производная по \alpha нужна)
	/// </summary>
	/// <param name="alpha">  параметр преобразования Фурье  </param>
	/// <returns></returns>
	OdeSolver<complex<double>> bigger_system(complex<double> alpha) const;


	complex<double> num1(OdeSolver<complex<double>>& smaller, double x2, double xi_2) const;
public:

	complex<double> kernel_transformant(complex<double> alpha) const;

	/// <summary>
	/// конструктор
	/// </summary>
	/// <param name="mu"> распределение модуля сдвига </param>
	/// <param name="rho"> распределение плотности </param>
	/// <param name="kappa"> безразмерная частота колебаний </param>
	/// <param name="roots_number"> максимальное количество корней </param>
	/// <param name="tolerance"> точность расчётов </param>
	layer(const function<double(double)>& mu,
		const function<double(double)>& rho, double kappa,
		size_t roots_number = MAX_ROOTS_NUMBER, double tolerance = 0.1e-5);


	vector<complex<double>> get_roots() const;

	/// <summary>
	/// фундаментальное решение
	/// </summary>
	/// <param name="x1"></param>
	/// <param name="x2"></param>
	/// <param name="xi_1"></param>
	/// <param name="xi_2"></param>
	/// <returns></returns>
	complex<double> fundamental_solution(double x1, double x2,
		double xi_1, double xi_2) const;


	complex<double> right_part(double xi_1, double xi_2) const;


	complex<double> b_constant(double uu) const;

};
#endif /* layer_hpp */

