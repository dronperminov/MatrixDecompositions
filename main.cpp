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
			U(i, j) = A(i, j);
			L(j, i) = A(j, i);

			for (int k = 0; k < i; k++) {
				U(i, j) -= L(i, k) * U(k, j);
				L(j, i) -= L(j, k) * U(k, i);
			}

			L(j, i) /= U(i, i);
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
		L(i, i) = A(i, i);

		for (int k = 0; k < i; k++)
			L(i, i) -= L(i, k) * L(i, k);

		// корень из отрицательных чисел извлечь нельзя
		if (L(i, i) < 0) {
			cout << "Unable to find Cholesky decomposition" << endl;
			return;
		}

		L(i, i) = sqrt(L(i, i));

		for (int j = i + 1; j < n; j++) {
			L(j, i) = A(j, i);

			for (int k = 0; k < i; k++)
				L(j, i) -= L(i, k) * L(j, k);

			L(j, i) /= L(i, i);
		}
	}

	cout << "Cholesky decomposition:" << endl;
	cout << "Matrix L:" << endl << L << endl;
	cout << "Matrix L^T:" << endl << L.Transpose() << endl;
	cout << "L * L^T:" << endl << (L * L.Transpose());
	cout << endl;
}

// QR разложение
// A — квадратная невырожденная матрица
// Q — унитарная матрица
// R — верхнетреугольная матрица
// A = Q*R
void QRDecomposition(Matrix A) {
	int n = A.size();

	Matrix Q(n);
	Matrix R(n);

	for (int j = 0; j < n; j++) {
		for (int i = 0; i < j; i++)
			for (int k = 0; k < n; k++)
				R(i, j) += A(k, j) * Q(k, i);

		for (int i = 0; i < n; i++) {
			Q(i, j) = A(i, j);

			for (int k = 0; k < j; k++)
				Q(i, j) -= R(k, j) * Q(i, k);

			R(j, j) += Q(i, j) * Q(i, j);
		}

		R(j, j) = sqrt(R(j, j));

		for (int i = 0; i < n; i++)
			Q(i, j) /= R(j, j); // нормируем вектор
	}

	cout << "QR decomposition:" << endl;
	cout << "Matrix Q:" << endl << Q << endl;
	cout << "Matrix R:" << endl << R << endl;
	cout << "Q * R:" << endl << (Q * R);
	cout << endl;
}

// LDL разложение
// A — симметричная положительно-определённая матрица
// L — нижняя треугольная матрица с единицами на диагонали
// D — диагональная матрица
// A = L*D*L^T
void LDLDecomposition(Matrix A) {
	int n = A.size();

	Matrix L(n);
	Matrix D(n);

	for (int j = 0; j < n; j++) {
		D(j, j) = A(j, j);
		L(j, j) = 1;

		for (int k = 0; k < j; k++)
			D(j, j) -= L(j, k) * L(j, k) * D(k, k);

		for (int i = j + 1; i < n; i++) {
			L(i, j) = A(i, j);

			for (int k = 0; k < j; k++)
				L(i, j) -= L(i, k) * L(j, k) * D(k, k);

			L(i, j) /= D(j, j);
		}
	}

	cout << "LDL decomposition:" << endl;
	cout << "Matrix L:" << endl << L << endl;
	cout << "Matrix D:" << endl << D << endl;
	cout << "Matrix L^T:" << endl << L.Transpose() << endl;
	cout << "L * D * L^T:" << endl << (L * D * L.Transpose());
	cout << endl;
}

int main() {
	Matrix A = ReadMatrix(); // считываем матрицу

	cout << endl << "Entered matrix: " << endl << A << endl; // выводим введённую матрицу

	LUDecomposition(A); // находим LU разложение
	CholeskyDecomposition(A); // находим разложение Холецкого
	QRDecomposition(A); // находим QR разложение
	LDLDecomposition(A); // находим LDL разложение
}