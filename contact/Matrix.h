#pragma once
#include <vector>
#include <stdexcept>
#include <iostream>


template<class T>
void gauss(std::vector<T>& first, std::vector<T>& second)
{
	size_t zeros = 0;
	while (zeros < first.size() && first[zeros] == 0 && second[zeros] == 0)
	{
		++zeros;
	}
	if (zeros == first.size())
	{
		return;
	}
	if (first[zeros] == 0)
	{
		std::swap(first, second);
	}
	if (second[zeros] == 0)
	{
		return;
	}
	auto a = first[zeros];
	auto b = second[zeros];
	for (size_t i = zeros; i < first.size(); i++)
	{
		second[i] *= a;
		second[i] -= b * first[i];
	}
}

template<class T>
bool equals_zero(const std::vector<T>& vec)
{
	/*for (auto x : vec)
	{
		if (x != 0)
		{
			return false;
		}
	}
	return true;*/
	return std::all_of(vec.begin(), vec.end(), [](T x) {return x != 0; });
}


template<class T>
class matrix
{
	size_t rows;
	size_t columns; // ÷èñëî ñòðîê è ñòîëáöîâ
	std::vector<std::vector<T>> elements; // âåêòîð ýëåìåíòîâ
	void setElements(const std::vector<std::vector<T>>& elements);

public:
	std::vector<std::vector<T>> getElements() const;

	size_t getRows() const;
	size_t getColumns() const;

	/**
	 * \brief êîíñòðóêòîð çàäà¸ò ðàçìåðû êâàäðàòíîé ìàòðèöû (âñå êîíñòðóêòîðû
	 * äîëæíû îáðàùàòüñÿ ê êîíñòðóêòîðàì áàçîâîãî êëàññà)
	 * \param rows êîëè÷åñòâî ñòðîê
	 * \param columns êîëè÷åñòâî ñòîëáöîâ
	 */
	matrix(size_t rows = 0, size_t columns = 0);
	matrix(size_t rows, size_t columns, const T& value);
	matrix(size_t rows, size_t columns,
		const std::vector<std::vector<T>>& elements);

	friend std::ostream& operator<<(std::ostream& os, const matrix& matrix)
	{
		for (auto element : matrix.elements)
		{
			os << "{ " << element[0];
			for (size_t i = 1; i < element.size(); i++)
			{
				os << ", " << element[i];
			}
			os << " },\n";
		}
		return os;
	}
	const std::vector<T>& operator[](size_t i) const;
	std::vector<T>& operator[](size_t i);

	matrix& operator+=(const matrix& that);
	matrix operator+(const matrix& that) const;

	matrix& operator-=(const matrix& that);
	matrix operator-(const matrix& that) const;

	matrix& operator*=(const matrix& that);
	matrix operator*(const matrix& that) const;

	matrix transpose() const;

	size_t rang() const;
};

template <class T>
size_t matrix<T>::getRows() const
{
	return rows;
}

template <class T>
size_t matrix<T>::getColumns() const
{
	return  columns;
}

template <class T>
matrix<T>::matrix(size_t rows, size_t columns) : rows(rows),
columns(columns)
{
	elements.resize(rows, std::vector<T>(columns));
}

template <class T>
matrix<T>::matrix(size_t rows, size_t columns, const T& value) :
	matrix(rows, columns)
{
	elements.resize(rows, std::vector<T>(columns, value));
}

template <class T>
matrix<T>::matrix(size_t rows, size_t columns,
	const std::vector<std::vector<T>>& elements) : matrix(rows, columns)
{
	setElements(elements);
}

template <class T>
const std::vector<T>& matrix<T>::operator[](size_t i) const
{
	return elements[i];
}

template <class T>
std::vector<T>& matrix<T>::operator[](size_t i)
{
	return elements[i];
}

template <class T>
matrix<T>& matrix<T>::operator+=(const matrix& that)
{
	if (this->rows != that.rows || this->columns != that.columns)
	{
		throw std::invalid_argument("");
	}
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < columns; j++)
		{
			elements[i][j] += that[i][j];
		}
	}
	return *this;
}

template <class T>
matrix<T> matrix<T>::operator+(const matrix& that) const
{
	auto temp(*this);
	return temp += that;
}

template <class T>
matrix<T>& matrix<T>::operator-=(const matrix& that)
{
	if (this->rows != that.rows || this->columns != that.columns)
	{
		throw std::invalid_argument("");
	}
	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < columns; j++)
		{
			elements[i][j] -= that[i][j];
		}
	}
	return *this;
}

template <class T>
matrix<T> matrix<T>::operator-(const matrix& that) const
{
	auto temp(*this);
	return temp -= that;
}

template <class T>
matrix<T>& matrix<T>::operator*=(const matrix& that)
{
	if (this->columns != that.rows)
	{
		throw std::invalid_argument("");
	}
	std::vector<std::vector<T>> temp(rows);
	for (size_t i = 0; i < rows; i++)
	{
		temp[i].resize(that.columns);
		for (size_t j = 0; j < that.columns; j++)
		{
			temp[i][j] = 0;
			for (size_t k = 0; k < columns; k++)
			{
				temp[i][j] += elements[i][k] * that[k][j];
			}
		}
	}
	elements = temp;
	return *this;
}

template <class T>
matrix<T> matrix<T>::operator*(const matrix& that) const
{
	auto temp(*this);
	return temp *= that;
}

template<class T>
matrix<T> matrix<T>::transpose() const
{
	matrix<T> result(this->getColumns(), this->getRows());
	for (size_t i = 0; i < this->getRows(); i++)
	{
		for (size_t j = 0; j < this->getColumns(); j++)
		{
			result[i][j] = this->getElements()[j][i];
		}
	}
	return result;
}

template <class T>
size_t matrix<T>::rang() const
{
	auto result = rows < columns ? rows : columns;
	auto temp = elements;
	for (size_t i = 0; i < rows - 1; i++)
	{
		for (size_t j = i + 1; j < rows; j++)
		{
			gauss(temp[i], temp[j]);
		}
	}
	while (equals_zero(temp[result - 1]))
	{
		--result;
	}
	return result;
}

template <class T>
void matrix<T>::setElements(const std::vector<std::vector<T>>& elements)
{
	if (elements.size() != rows)
	{
		throw std::invalid_argument("");
	}
	for (auto element : elements)
	{
		if (elements.size() != columns)
		{
			throw std::invalid_argument("");
		}
	}
	this->elements = elements;
}

template <class T>
std::vector<std::vector<T>> matrix<T>::getElements() const
{
	return elements;
}