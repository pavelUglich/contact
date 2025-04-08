#pragma once
#include "Matrix.h"

template<class T>
class square_matrix : public matrix<T>
{
	std::vector<T> transform_the_right_part(
		const std::vector<T>& rightPart, const std::vector<size_t>& p) const;
	std::vector<T> there_and_back_again(const std::vector<std::vector<T>>& lu,
		const std::vector<T>& rp) const;


public:
	std::vector<std::vector<T>> LUDecompose(std::vector<size_t>& p) const;

	/**
	 * \brief êîíñòðóêòîð çàäà¸ò ðàçìåðû êâàäðàòíîé ìàòðèöû (âñå êîíñòðóêòîðû
	 * äîëæíû îáðàùàòüñÿ ê êîíñòðóêòîðàì áàçîâîãî êëàññà)
	 * \param rows êîëè÷åñòâî ñòðîê
	 */
	square_matrix(size_t rows);

	/**
	 * \brief êîíñòðóêòîð çàäà¸ò ðàçìåðû êâàäðàòíîé ìàòðèöû è çàïîëíÿåò å¸
	 * çíà÷åíèì value
	 * \param rows êîëè÷åñòâî ñòðîê
	 * \param value çíà÷åíèå
	 */
	square_matrix(size_t rows, const T& value);

	/**
	 * \brief êîíñòðóêòîð çàäà¸ò ðàçìåðû êâàäðàòíîé ìàòðèöû è çàïîëíÿåò å¸
	 * çíà÷åíèÿìè èç âåêòîðà
	 * \param rows ðàçìåð ìàòðèöû
	 * \param elements âåêòîð çíà÷åíèé
	 */
	square_matrix(size_t rows, const std::vector<std::vector<T>>& elements);

	square_matrix(const matrix<T>& matrix);


	/**
	 * \brief îïðåäåëèòåëü ìàòðèöû
	 * \return îïðåäåëèòåëü ìàòðèöû
	 */
	T det() const;

	/**
	 * \brief ðåøåíèå ñèñòåìû ëèíåéíûõ àëãåáðàè÷åñêèõ óðàâíåíèé
	 * \param rightPart âåêòîð ïðàâîé ÷àñòè
	 * \return ðåøåíèå ÑËÀÓ
	 */
	std::vector<T> linSolve(const std::vector<T>& rightPart) const;

	/**
	 * \brief îáðàùåíèå ìàòðèöû
	 * \return îáðàòíàÿ ìàòðèöà
	 */
	square_matrix operator~() const;

	/**
	 * \brief äåëåíèå íà äðóãóþ ìàòðèöó (óìíîæåíèå ñïðàâà íà ìàòðèöó, îáðàòíóþ
	 * ê äðóãîé)
	 * \param squareMatrix
	 * \return
	 */
	square_matrix& operator/=(const square_matrix<T>& squareMatrix);
	square_matrix operator/(const square_matrix<T>& squareMatrix) const;
};

template <class T>
square_matrix<T>::square_matrix(size_t rows) : matrix<T>(rows, rows)
{
}

template <class T>
square_matrix<T>::square_matrix(size_t rows, const T& value) :
	matrix<T>(rows, rows, value)
{
}

template <class T>
square_matrix<T>::square_matrix(size_t rows,
	const std::vector<std::vector<T>>& elements) :
	matrix<T>(rows, rows, elements)
{
}

template<class T>
square_matrix<T>::square_matrix(const matrix<T>& matrix) :matrix<T>(
	matrix.getRows(), matrix.getColumns(), matrix.getElements())
{
	if (matrix.getRows() != matrix.getColumns())
	{
		throw std::invalid_argument("");
	}
}

template<class T>
T square_matrix<T>::det() const
{
	std::vector<size_t> p(this->getRows());
	auto lu = LUDecompose(p);
	T result = 1;
	for (size_t i = 0; i < p.size(); i++)
	{
		result *= lu[i][i];
	}
	int parity = 0;
	for (size_t i = 0; i < p.size(); i++)
	{
		for (size_t ii = i + 1; ii < p.size(); ii++)
		{
			if (p[i] > p[ii])
			{
				++parity;
			}
		}
	}
	return result * (parity % 2 == 0 ? 1 : -1);
}

template <class T>
std::vector<std::vector<T>> square_matrix<T>::LUDecompose(
	std::vector<size_t>& p) const
{
	auto elems = this->getElements();
	p.resize(this->getRows());
	for (size_t i = 0; i < this->getRows(); i++)
		p[i] = i;
	for (size_t i = 0; i < this->getRows(); i++) {
		auto max = abs(elems[i][i]);
		size_t imax = i;
		for (size_t k = i; k < this->getRows(); k++) {
			if (abs(elems[k][i]) > max) {
				max = abs(elems[k][i]);
				imax = k;
			}
		}
		if (max < DBL_EPSILON) throw std::invalid_argument("");
		if (imax != i) {
			std::swap(p[i], p[imax]);
			std::swap(elems[i], elems[imax]);
		}
		for (size_t j = i + 1; j < this->getRows(); j++) {
			elems[j][i] /= elems[i][i];
			for (size_t k = i + 1; k < this->getRows(); k++)
				elems[j][k] -= elems[j][i] * elems[i][k];
		}
	}
	return elems;
}

template <class T>
std::vector<T> square_matrix<T>::transform_the_right_part(
	const std::vector<T>& rightPart, const std::vector<size_t>& p) const
{
	std::vector<T> rp(this->getRows());
	for (size_t i = 0; i < this->getRows(); i++)
	{
		rp[i] = rightPart[p[i]];
	}
	return rp;
}

template <class T>
std::vector<T> square_matrix<T>::there_and_back_again(
	const std::vector<std::vector<T>>& lu,
	const std::vector<T>& rp) const
{
	std::vector<T> y(this->getRows());
	for (size_t i = 0; i < this->getRows(); i++)
	{
		T sum = 0;
		for (size_t j = 0; j < i; j++)
		{
			sum += y[j] * lu[i][j];
		}
		y[i] = rp[i] - sum;
	}
	for (int i = static_cast<int>(this->getRows() - 1); i >= 0; i--)
	{
		T sum = 0;
		for (size_t j = i + 1; j < this->getRows(); j++)
		{
			sum += y[j] * lu[i][j];
		}
		y[i] = (y[i] - sum) / lu[i][i];
	}
	return y;
}

template<class T>
std::vector<T>
square_matrix<T>::linSolve(const std::vector<T>& rightPart) const
{
	std::vector<size_t> p(this->getRows());
	auto lu = LUDecompose(p);
	auto rp = transform_the_right_part(rightPart, p);
	return there_and_back_again(lu, rp);
}

template<class T>
square_matrix<T> square_matrix<T>::operator~() const
{
	square_matrix<T> result(this->getRows());
	std::vector<size_t> p(this->getRows());
	auto lu = LUDecompose(p);
	for (size_t i = 0; i < this->getRows(); i++)
	{
		std::vector<T> right_part(this->getRows(), 0);
		right_part[i] = 1;
		auto rp = transform_the_right_part(right_part, p);
		auto column = there_and_back_again(lu, rp);
		for (size_t ii = 0; ii < this->getRows(); ii++)
		{
			result[ii][i] = column[ii];
		}
	}
	return result;
}

template<class T>
square_matrix<T>& square_matrix<T>::operator/=(
	const square_matrix<T>& squareMatrix)
{
	*this *= ~squareMatrix;
	return *this;
}

template<class T>
square_matrix<T> square_matrix<T>::operator/(
	const square_matrix<T>& squareMatrix) const
{
	auto temp(*this);
	return temp /= squareMatrix;
}