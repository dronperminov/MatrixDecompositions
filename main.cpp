#include <iostream>

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
	for (int j = 0; j < n; j++) {
		U(0, j) = A(0, j);
		L(j, 0) = A(j, 0) / U(0, 0);
	}

	for (int i = 1; i < n; i++) {
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

int main() {
	Matrix A = ReadMatrix(); // считываем матрицу

	cout << endl << "Entered matrix: " << endl << A << endl; // выводим введённую матрицу

	LUDecomposition(A); // находим LU разложение
}