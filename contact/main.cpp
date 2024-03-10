#include "layer.hpp"
#include <numeric>
#include"plots.h"
#include "Integral.h"
#include "constant_element.h"
#include "linear_element.h"

const double pi = acos(-1.0);

complex<double> transformant(complex<double> alpha, double _kappa)
{
	//return 1.0 / alpha;
	/*const auto _mu = [](double) {return 1.0; };
	const auto _rho = [](double) {return 1.0; };
	OdeSolver<complex<double>> os =
	{
		{
			[=](double x, const vector<complex<double>>& y) {return y[1] / _mu(x); },
			[=](double x, const vector<complex<double>>& y)
			{
				return (alpha * alpha * _mu(x) - _rho(x) * _kappa * _kappa) * y[0];
			}
		},
		0.1e-6,
		RUNGE_KUTTA_FELDBERG
	};
	const auto solution = os.solve(0, 1, { 0,exp(-abs(alpha)) });
	return solution.front() / solution.back() / alpha;//*/
	const auto gamma = sqrt(alpha * alpha - _kappa * _kappa);
	return sinh(gamma) / gamma / cosh(gamma);
}

complex<double> analytical(double x, double _kappa, double eps = 0.1e-5)
{
	complex<double> result = 0;
	size_t n = 0;
	complex<double> term;
	do
	{
		const complex<double> frac = (2 * n++ + 1) / 2.0 * pi;
		const complex<double> alpha = sqrt(_kappa * _kappa - (frac * frac));
		term = exp(I * alpha * x) / alpha;
		result += term / alpha;
	} while (abs(term) > eps);
	return result;
}

template<class Fun>
map<double, vector<complex<double>>> dispersionalSet(double theHighestFrequency, Fun rho, Fun mu)
{
	map<double, vector<complex<double>>> result;
	double kappa = 0.1;
	const double step = 0.05;

	layer l(mu, rho, kappa, 5);
	auto roots = l.get_roots();
	while (kappa < theHighestFrequency)
	{
		cout << kappa << endl;
		result[kappa] = roots;
		kappa += step;
		l = { mu, rho, kappa, 7 };
		roots = l.get_roots();
	}
	return result;
}



void plotTheDispersionalCurves(
	map<double, vector<complex<double>>>& dispersionSet,
	const string& fileName)
{
	ofstream os(fileName);
	os << "\\begin{tikzpicture}[scale=1.5]\n";
	os << "\\begin{axis}[grid]\n";
	const auto numberOfCurves = dispersionSet.rbegin()->second.size();
	for (auto [x, y] : dispersionSet) {
		for (auto& i : y) {
			if (i.real() < DBL_EPSILON)
				i = { -i.imag(),0 };
		}
		std::ranges::sort(y.begin(), y.end(),
			[](auto x, auto y) { return x.real() < y.real(); });
	}
	//const auto maxFreq = dispersionSet.rbegin()->first;
	for (size_t i = 0; i < numberOfCurves; i++) {
		os << "\\addplot[smooth, black] plot coordinates{\n";
		for (auto [x, y] : dispersionSet) {
			if (i < y.size())
			{
				os << "(" << y[i].real() << ", " << x << ") ";
			}
		}
		os << "};\n";
	}
	os << "\\end{axis}\n";
	os << "\\end{tikzpicture}\n";
	os.close();
}


void surface_plot(const std::vector<double>& x_1, const std::vector<double>& x_2, 
	const std::vector<std::complex<double>>& solution, const layer& l, 
	const std::string& file_name) 
{
	const auto matrix = l.displacement(x_1, x_2, solution);
	std::ofstream ofs(file_name);
	for (size_t i = 0; i < x_2.size(); i++)
	{
		for (size_t ii = 0; ii < x_1.size(); ii++)
		{
			ofs <<"(" << x_1[ii] << ", " << x_2[i] << ", " << abs(matrix[i][ii]) << ") ";
		}
	}
	ofs.close();
}

int main() {
	const auto mu = [](double x) {return 1.0 + x; };
	const auto rho = [](double x) {return 1.0 + exp(-x); };
	const auto mass = 1.0;// масса штампа
	const auto kappa = 2.0; // волновое число
	const auto tt = 1.0; // интенсивность нагрузки	
	const double y = 1.0;
	size_t size = 20;
	std::vector<std::complex<double>> wavefield(size);
	std::vector<std::complex<double>> wavefield1(size + 1);
	const double a = 1.0;
	const auto h = 2 * a / (double)size;
	std::vector<boundary_element*> elements(size);
	double left = -a;
	for (size_t i = 0; i < size; i++)
	{
		double right = left + h;
		elements[i] = new constant_element(left, right);
		left = right;
	}
	vector<complex<double>> sol, sol1;
	{
		layer l = { mu, rho, "roots.txt", elements };
		const auto mmm = l.matrix();
		auto matrix = SquareMatrix<complex<double>>(size, mmm);
		sol = matrix.linSolve(vector<complex<double>>(size, { 1,0 }));
		const complex<double> integral_ = h * std::accumulate(sol.begin(), sol.end(), complex<double>());
		const auto gamma = (integral_ - tt) / mass / kappa / kappa;
		std::ranges::transform(sol.begin(), sol.end(), sol.begin(), [=](auto x) {return x / gamma; });
		for (size_t i = 0; i < size; i++)
		{
			const double x = -a + (i + 0.5) * h;
			wavefield[i] = l.displacement(x, y, sol);
		}
	}
	std::vector<double> nodes(size + 1);
	for (size_t i = 0; i < nodes.size(); i++)
	{
		nodes[i] = -a + i * h;
	}
	left = -a;
	for (size_t i = 0; i < size; i++)
	{
		double right = left + h;
		elements[i] = new linear_element(nodes[i], nodes[i + 1]);
		left = right;
	}

	{
		layer l = { mu, rho, "roots.txt", elements };
		const auto mmm = l.matrix();
		auto matrix = SquareMatrix<complex<double>>(size + 1, mmm);
		sol1 = matrix.linSolve(vector<complex<double>>(size + 1, { 1,0 }));
		const complex<double> integral_ = h * (std::accumulate(sol1.begin() + 1, sol1.end() - 1, complex<double>())
			+ 0.5 * (sol1.front() + sol1.back()));
		const auto gamma = (integral_ - tt) / mass / kappa / kappa;
		std::ranges::transform(sol1.begin(), sol1.end(), sol1.begin(), [=](auto x) {return x / gamma; });
		const auto dd = l.displacement({ -1.0,0.0,1.0 }, { 0.25, 0.5, 1.0 }, sol1);
		surface_plot({ -1.0,0.0,1.0 }, { 0.25, 0.5, 1.0 }, sol1, l, "1.dat");
	}//*/
	//plotTheWaveField(h, { {"red", sol}, {"blue", sol1} }, "wf.txt", -a + h / 2);
	plotTheWaveField(h, { {"red", wavefield}, {"blue", wavefield1} }, "wavefield.txt", -a + h / 2);
	system("pause");
}