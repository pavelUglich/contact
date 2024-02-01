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



const size_t MAX_ROOTS_NUMBER = 50;
const complex<double> I = { 0, 1 };
/// <summary>
/// класс - абстрактный. Метод matrix - чисто абстрактный, реализация зависит 
/// от типа граничного элемента
/// </summary>
class layer {
	// распределения модуля сдвига и плотности
	function<double(double)> _mu;
	function<double(double)> _rho;

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
	/// интегральный множитель в формуле (27)
	/// </summary>
	/// <param name="tau"></param>
	/// <param name="a"></param>
	/// <param name="size"></param>
	/// <param name="alpha"></param>
	/// <returns></returns>
	virtual std::complex<double> evaluate_exp(
		const std::vector<std::complex<double>>& tau, double a, double x_1, 
		size_t size, std::complex<double> alpha) const = 0;

protected:
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

	/// <summary>
	/// знаметатели вычетов в подынтегральных выражениях
	/// </summary>
	/// <returns>вектор комплексных чисел</returns>
	vector<complex<double>> get_residues() const;

	/// <summary>
	/// точность расчёта
	/// </summary>
	/// <returns>точность</returns>
	double get_tolerance() const;

public:

	/// <summary>
	/// найти корни (метод открытый, чтобы осталась возможность построения дисперсионных кривых)
	/// </summary>
	/// <returns>вектор комплексных чисел</returns>
	vector<complex<double>> get_roots() const;

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


	/// <summary>
	/// правая часть интегрального уравнения
	/// </summary>
	/// <param name="xi_1"></param>
	/// <param name="xi_2"></param>
	/// <returns></returns>
	complex<double> right_part(double xi_1, double xi_2) const;


	/// <summary>
	/// матрица для решения контактной задачи
	/// </summary>
	/// <param name="a"></param>
	/// <returns></returns>
	virtual vector<vector<complex<double>>> matrix(double a, size_t size) 
		const = 0;
		
	/// <summary>
	/// перемещение 
	/// </summary>
	/// <param name="x_1">координата x_1</param>
	/// <param name="x_2">координата x_2</param>
	/// <param name="solution">решение контактной задачи</param>
	/// <param name="a">полуширина зоны контакта</param>
	/// <returns>перемещение</returns>
	virtual std::complex<double> displacement(double x_1, double x_2, 
		std::vector<std::complex<double>> solution, double a) const = 0;

	/// <summary>
	/// вспомогательный метод для отыскания поля перемещений (первый сомножитель в
	/// формуле (27), вычисляется отдельно для каждого корня)
	/// </summary>
	/// <param name="x_2"></param>
	/// <param name="alpha"></param>
	/// <returns></returns>
	std::complex<double> displacement_auxiliary(double x_2, 
		std::complex<double> alpha) const;

	std::vector<std::complex<double>> displacement_auxiliary(const std::vector<double>& x_2,
		std::complex<double> alpha) const;

	void process_points_vector(std::vector<double>& points) const;


	virtual ~layer() = default;
};
#endif /* layer_hpp */

