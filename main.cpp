#include <iostream>
#include <cmath>

#include "Matrix.h"

using namespace std;

// считывание матрицы
Matrix ReadMatrix() {
	int n; // размер матрицы
	cout << "Enter size of matrix: ";
	cin >> n; // считываем размер матрицы

	Matrix matrix(n); // создаём матрицу для поиска собственных значений

	cout << "Enter values of matrix:" << endl;
	matrix.Read(); // считываем значения матрицы

	return matrix; // возвращаем матрицу
}

// LU разложение
// A — квадратная обратимая матрица
// L — нижняя треугольная матрица
// U — верхняя треугольная матрица)
// A = L*U
void LUDecomposition(Matrix A) {
	int n = A.size();
	
	Matrix L(n);
	Matrix U(n);

	// выполняем LU разложение матрицы
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			double sum = 0;

			for (int k = 0; k < i; k++)
				sum += L(i, k) * U(k, j);

			U(i, j) = A(i, j) - sum;
		
			sum = 0;

			for (int k = 0; k < i; k++)
				sum += L(j, k) * U(k, i);

			L(j, i) = (A(j, i) - sum) / U(i, i);
		}
	}

	cout << "LU decomposition:" << endl;
	cout << "Matrix L:" << endl << L << endl;
	cout << "Matrix U:" << endl << U << endl;
	cout << "L * U:" << endl << (L * U);
	cout << endl;
}

// Разложение Холецкого
// A — симметричная положительно-определённая матрица
// L — нижняя треугольная матрица со строго положительными элементами на диагонали
// A = L*L^T
void CholeskyDecomposition(Matrix A) {
	int n = A.size();

	Matrix L(n);

	for (int i = 0; i < n; i++) {
		double sum = A(i, i);

		for (int k = 0; k < i; k++)
			sum -= L(i, k) * L(i, k);

		// корень из отрицательных чисел извлечь нельзя
		if (sum < 0) {
			cout << "Unable to find Cholesky decomposition" << endl;
			return;
		}

		L(i, i) = sqrt(sum);

		for (int j = i; j < n; j++) {
			sum = 0;

			for (int k = 0; k < i; k++)
				sum += L(i, k) * L(j, k);

			L(j, i) = (A(j, i) - sum) / L(i, i);
		}
	}

	cout << "Cholesky decomposition:" << endl;
	cout << "Matrix L:" << endl << L << endl;
	cout << "Matrix L^T:" << endl << L.Transpose() << endl;
	cout << "L * L^T:" << endl << (L * L.Transpose());
	cout << endl;
}

int main() {
	Matrix A = ReadMatrix(); // считываем матрицу

	cout << endl << "Entered matrix: " << endl << A << endl; // выводим введённую матрицу

	LUDecomposition(A); // находим LU разложение
	CholeskyDecomposition(A); // находим разложение Холецкого
}