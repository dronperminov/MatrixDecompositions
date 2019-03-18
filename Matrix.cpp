#include "Matrix.h"

// конструктор матрицы из размера
Matrix::Matrix(int n) {
	this->n = n; // сохраняем размер

	// выделяем память под значения матрицы
	values = new double*[n];

	for (int i = 0; i < n; i++) {
		values[i] = new double[n];

		for (int j = 0; j < n; j++)
			values[i][j] = 0; // обнуляем все элементы
	}
}

// конструктор копирования
Matrix::Matrix(const Matrix& matrix) {
	n = matrix.n;

	// выделяем память под значения матрицы
	values = new double*[n];

	for (int i = 0; i < n; i++) {
		values[i] = new double[n];

		for (int j = 0; j < n; j++)
			values[i][j] = matrix.values[i][j]; // копируем значения матрицы
	}
}

// оператор присваивания
Matrix& Matrix::operator=(const Matrix& matrix) {
	if (this == &matrix)
		return *this;

	for (int i = 0; i < n; i++)
		delete[] values[i];

	delete[] values;

	n = matrix.n;

	// выделяем память под значения матрицы
	values = new double*[n];

	for (int i = 0; i < n; i++) {
		values[i] = new double[n];

		for (int j = 0; j < n; j++)
			values[i][j] = matrix.values[i][j]; // копируем значения матрицы
	}

	return *this;
}

// получение размера матрицы
int Matrix::size() const {
	return n;
}

// получение элемента по индексу
double Matrix::operator()(int i, int j) const {
	return values[i][j];
}

// получение элемента по индексу
double& Matrix::operator()(int i, int j) {
	return values[i][j];
}

// умножение матрицы на матрицу
Matrix Matrix::operator*(const Matrix &matrix) const {
	Matrix result(n);

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double sum = 0;

			for (int k = 0; k < n; k++)
				sum += values[i][k] * matrix.values[k][j];

			result.values[i][j] = sum;
		}
	}

	return result;
}

// заполнение матрицы с клавиатуры
void Matrix::Read() {
	for (int i = 0; i < n; i++) {
		std::cout << "Enter row " << (i + 1) << ": ";

		for (int j = 0; j < n; j++)
			std::cin >> values[i][j];
	}
}

// получение транспонированной матрицы
Matrix Matrix::Transpose() const {
	Matrix T(n); // создаём матрицу

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			T.values[i][j] = values[j][i]; // копируем транспонированные элементы

	return T; // возвращаем матрицу
}

// деструктор (освобождение памяти)
Matrix::~Matrix() {
	for (int i = 0; i < n; i++)
		delete[] values[i];

	delete[] values;
}

// вывод в поток
std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
	for (int i = 0; i < matrix.n; i++) {
		for (int j = 0; j < matrix.n; j++) {
			os << std::setw(6); 
			os << (fabs(matrix.values[i][j]) > 1e-15 ? matrix.values[i][j] : 0);
			os << " ";
		}

		os << std::endl;
	}

	return os;
}