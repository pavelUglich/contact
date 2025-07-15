#include "plots.h"
#include <string>
#include <fstream>
#include <iostream>

void showVector(std::vector<double>& complexValuedVector, std::ostream& stream)
{
	const auto points = complexValuedVector.size();
	for (size_t i = 0; i < points; i++)
	{
		stream << "(" << (i + 0.5) / points << "," << complexValuedVector[i] << ") ";
	}
	std::cout << std::endl;
}

void showVector(const std::vector<double>& complexValuedVector, std::ostream& stream,
	double step, double min)
{
	for (size_t i = 0; i < complexValuedVector.size(); i++)
	{
		const auto x = min + (i + 0.5) * step;
		stream << "(" << x << ", " << abs(complexValuedVector[i]) << ") ";
	}
}

void addTheCurve(const std::vector<double>& fieldHomogeneous, std::ofstream& stream, const std::string& color,
	double step, double min)
{
	stream << "	\\addplot[line width = 0.25mm, smooth, ";
	stream << color;
	stream << "] plot coordinates{\n";
	showVector(fieldHomogeneous, stream, step, min);
	stream << "					};\n";
}

void plot_the_wave_field(const std::map<std::string, std::vector<double>>& vectors, const std::string& fileName,
	double step, double min)
{
	std::ofstream stream(fileName);
	stream << "\\begin{tikzpicture}[scale=1.5]\n";
	stream << "\\begin{axis}[grid]\n";
	for (const auto& item : vectors)
	{
		addTheCurve(item.second, stream, item.first, step, min);
	}
	stream << "\\end{axis}\n";
	stream << "\\end{tikzpicture}\n";
	stream.close();
}

void addTheCurve(double step, 
	const std::vector<std::complex<double>>& fieldHomogeneous, 
	std::ofstream& stream, const std::string& color, double a)
{
	stream << " \\addplot[line width = 0.4mm, smooth, ";
	stream << color;
	stream << "] plot coordinates{\n";
	showVector(step, fieldHomogeneous, [](auto x) { return x.real(); }, stream, a);
	stream << " };\n";
	stream << " \\addplot[line width = 0.4mm, smooth, dashed,";
	stream << color;
	stream << "] plot coordinates{\n";
	showVector(step, fieldHomogeneous, [](auto x) { return x.imag(); }, stream, a);
	stream << " };\n";
}

void addTheCurve(const std::vector<double>& nodes, 
	const std::vector<std::complex<double>>& fieldHomogeneous, 
	std::ofstream& stream, const std::string& color)
{
	stream << " \\addplot[line width = 0.25mm, smooth, ";
	stream << color;
	stream << "] plot coordinates{\n";
	showVector(nodes, fieldHomogeneous, [](auto x) { return x.real(); }, stream);
	stream << " };\n";
	stream << " \\addplot[smooth, dashed,";
	stream << color;
	stream << "] plot coordinates{\n";
	showVector(nodes, fieldHomogeneous, [](auto x) { return x.imag(); }, stream);
	stream << " };\n";
}

void plot_the_wave_field(double step,
	const std::map<std::string, std::vector<std::complex<double>>>& vectors,
	const std::string& fileName, double a)
{
	std::ofstream stream(fileName);
	stream << "\\begin{tikzpicture}[scale=1.1]\\footnotesize\n";
	stream << "\\begin{axis}[grid]\n";
	for (const auto& item : vectors)
	{
		addTheCurve(step, item.second, stream, item.first, a);
	}
	stream << "\\end{axis}\n";
	stream << "\\end{tikzpicture}\n";
	stream.close();
}

void plot_the_wave_field(const std::vector<double>& nodes, 
	const std::map<std::string, std::vector<std::complex<double>>>& vectors, 
	const std::string& fileName)
{
	std::ofstream stream(fileName);
	stream << "\\begin{tikzpicture}[scale=1.1]\\tiny\n";
	stream << "\\begin{axis}[grid]\n";
	for (const auto& item : vectors)
	{
		addTheCurve(nodes, item.second, stream, item.first);
	}
	stream << "\\end{axis}\n";
	stream << "\\end{tikzpicture}\n";
	stream.close();
}

void showVector(double step, const std::vector<std::complex<double>>& complexValuedVector,
	const std::function<double(std::complex<double>)>& mapping, std::ostream& stream, double a)
{
	for (size_t i = 0; i < complexValuedVector.size(); i++)
	{
		const auto x = a + i * step;
		stream << "(" << x << ", " << mapping(complexValuedVector[i]) << ") ";
	}
}

void showVector(const std::vector<double>& nodes, 
	const std::vector<std::complex<double>>& complexValuedVector, 
	const std::function<double(std::complex<double>)>& mapping, 
	std::ostream& stream)
{
	for (size_t i = 0; i < complexValuedVector.size(); i++)
	{
		stream << "(" << nodes[i] << ", " << mapping(complexValuedVector[i]) << ") ";
	}
}

