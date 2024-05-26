#pragma once
#include <complex>
#include <fstream>
#include <functional>
#include <vector>
#include <map>

void showVector(std::vector<double>& complexValuedVector,
                std::ostream& stream);


void showVector(const std::vector<double>& complexValuedVector,
                std::ostream& stream, double step, double min = 0);




/**
 * \brief добавить кривую на график
 * \param step шаг
 * \param complexValuedVector вектор
 * \param mapping отображение
 * \param stream файловый поток
 * \param a левая граница отрезка
 */
void showVector(double step,
    const std::vector<std::complex<double>>& complexValuedVector,
    const std::function<double(std::complex<double>)>& mapping, 
    std::ostream& stream, double a = 0);

void showVector(const std::vector<double>& nodes,
    const std::vector<std::complex<double>>& complexValuedVector,
    const std::function<double(std::complex<double>)>& mapping,
    std::ostream& stream);


void addTheCurve(const std::vector<double>& fieldHomogeneous, std::ofstream& stream,
                 const std::string& color, double step, double min = 0);



void plotTheWaveField(const std::map<std::string, std::vector<double>>& vectors, const std::string& fileName,
                      double step, double min = 0);


/**
 * \brief  добавить кривую на график
 * \param step шаг
 * \param fieldHomogeneous вектор
 * \param stream файловый поток
 * \param color цвет
 * \param a левая граница
 */
void addTheCurve(double step, 
    const std::vector<std::complex<double>>& fieldHomogeneous, 
    std::ofstream& stream, const std::string& color, double a = 0);

void addTheCurve(const std::vector<double>& nodes,
    const std::vector<std::complex<double>>& fieldHomogeneous,
    std::ofstream& stream, const std::string& color);


/**
 * \brief построить диаграмму
 * \param step шаг
 * \param vectors набор векторов
 * \param fileName имя файла
 * \param a левая граница
 */
void plotTheWaveField(double step, const std::map<std::string, std::vector<std::complex<double>>>& vectors,
                      const std::string& fileName, double a = 0);

/**
 * \brief построить диаграмму
 * \param step шаг
 * \param vectors набор векторов
 * \param fileName имя файла
 * \param a левая граница
 */
void plotTheWaveField(const std::vector<double>& nodes, 
    const std::map<std::string, std::vector<std::complex<double>>>& vectors, 
    const std::string& fileName);
