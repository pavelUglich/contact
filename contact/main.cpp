#include "layer.hpp"
#include <numeric>
#include"plots.h"
#include "layer_constant.h"
#include "layer_linear.h"
#include "Integral.h"

const double pi = acos(-1.0);
enum class element_type { constant, linear };

layer* make_layer(element_type et, const function<double(double)>& mu,
	const function<double(double)>& rho, double kappa,
	size_t roots_number = MAX_ROOTS_NUMBER, double tolerance = 0.1e-5)
{
	switch (et)
	{
	case element_type::constant:
		return new layer_constant(mu, rho, kappa, roots_number, tolerance);
	case element_type::linear:
		return new layer_linear(mu, rho, kappa, roots_number, tolerance);
	default:
		break;
	}
	throw std::exception();
}

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

	layer_constant l(mu, rho, kappa, 5);
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


int main() {
	const auto mu = [](double x) {return 1.0 + x; };
	const auto rho = [](double x) {return 1.0 + exp(-x); };
	const auto mass = 1.0;// масса штампа
	const auto kappa = 5.0; // волновое число
	const auto tt = 1.0; // интенсивность нагрузки	
	const double y = .001;
	size_t size = 20;
	//const std::vector<std::complex<double>> sol(size, { 1.0, 0 });
	//const std::vector<std::complex<double>> sol1(size + 1, { 1.0, 0 });
	std::vector<std::complex<double>> wavefield(size);
	std::vector<std::complex<double>> wavefield1(size + 1);
	const double a = 1.0;
	const auto h = 2 * a / (double)size;
	{
		layer_constant l = { mu, rho, "roots.txt"/*kappa*/ };
		const auto mm = l.matrix(a, size);
		auto matrix = SquareMatrix<complex<double>>(size, mm);
		auto sol = matrix.linSolve(vector<complex<double>>(size, { 1,0 }));
		const complex<double> integral_ = h * std::accumulate(sol.begin(), sol.end(), complex<double>());
		const auto gamma = (integral_ - tt) / mass / kappa / kappa;
		std::ranges::transform(sol.begin(), sol.end(), sol.begin(), [=](auto x) {return x / gamma; });
		/*for (size_t i = 0; i < size; i++)
		{
			const double x = -a + (i + 0.5) * h;
			wavefield[i] = l.displacement(x, y, sol, a);
		}*/
		const auto rp = matrix * sol;
		cout << rp.front() << endl;
		cout << matrix.getElements().front().front() << endl;
		std::cout << l.displacement(-a + 0.5 * h, 1.0, sol, a) << std::endl;
	}
	/*
	{
		layer_linear l = { mu, rho, "roots.txt" };
		const auto mm = l.matrix(a, size);
		auto matrix = SquareMatrix<complex<double>>(size + 1, mm);
		auto sol = matrix.linSolve(vector<complex<double>>(size + 1, { 1,0 }));
		const complex<double> integral_ = h * (std::accumulate(sol.begin() + 1, sol.end() - 1, complex<double>())
			+ 0.5 * (sol.front() + sol.back()));
		const auto gamma = (integral_ - tt) / mass / kappa / kappa;
		std::ranges::transform(sol.begin(), sol.end(), sol.begin(), [=](auto x) {return x / gamma; });
		for (size_t i = 0; i <= size; i++)
		{
			const double x = -a + i * h;
			wavefield1[i] = l.displacement(x, y, sol, a);
		}
	}//*/
	//plotTheWaveField(h, { {"red", wavefield}, {"blue", wavefield1} }, "wavefield.txt", -a + h / 2);
	system("pause");
}