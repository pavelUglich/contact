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
#include "boundary_element.h"
using namespace std;


const size_t MAX_ROOTS_NUMBER = 50;
const complex<double> I = { 0, 1 };

std::vector<double> get_nodal_points(const std::vector<boundary_element*>& elements);

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

	// набор элементов
	std::vector<boundary_element*> _elements;

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
	vector<complex<double>> row(double node, complex<double> root) const;
	vector<vector<complex<double>>> row(const vector<double>& node, complex<double> root) const;
	vector<complex<double>> row(double node, size_t size, double y = 1.0) const;
	vector<vector<complex<double>>> row(double node, size_t size, 
		const vector<double>& y) const;
	vector<vector<vector<complex<double>>>> row(const vector<double>& node, 
		size_t size, const vector<double>& y) const;

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

	layer(const function<double(double)>& mu,
		const function<double(double)>& rho, const std::string& filename, 
		double tolerance = 0.1e-5);

	layer(const function<double(double)>& mu,
		const function<double(double)>& rho, double kappa,
		const std::vector<boundary_element* >& elements,
		size_t roots_number = MAX_ROOTS_NUMBER, double tolerance = 0.1e-5);

	layer(const function<double(double)>& mu,
		const function<double(double)>& rho, const std::string& filename,
		const std::vector<boundary_element* >& elements,
		double tolerance = 0.1e-5);

	/// <summary>
	/// правая часть интегрального уравнения
	/// </summary>
	/// <param name="xi_1"></param>
	/// <param name="xi_2"></param>
	/// <returns></returns>
	complex<double> right_part(double xi_1, double xi_2) const;

	vector<vector<complex<double>>> matrix() const;
		
	/// <summary>
	/// перемещение 
	/// </summary>
	/// <param name="x_1">координата x_1</param>
	/// <param name="x_2">координата x_2</param>
	/// <param name="solution">решение контактной задачи</param>
	/// <param name="a">полуширина зоны контакта</param>
	/// <returns>перемещение</returns>
	std::complex<double> displacement(double x_1, double x_2, 
		const std::vector<std::complex<double>>& solution) const;

	std::vector<std::complex<double>> displacement(double x_1, 
		const std::vector<double>& x_2,	
		const std::vector<std::complex<double>>& solution) const;

	std::vector<std::vector<std::complex<double>>> displacement(const std::vector<double>& x_1,
		const std::vector<double>& x_2,	
		const std::vector<std::complex<double>>& solution) const;
	void process_points_vector(std::vector<double>& points) const;
	void save_to_file(const std::string& filename) const;
	void load_from_file(const std::string& filename);

	virtual ~layer();;
};
#endif /* layer_hpp */

