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

// QR разложение методом вращений (Гивенса)
// A — квадратная невырожденная матрица
// Q — унитарная матрица
// R — верхнетреугольная матрица
// A = Q*R
void QRDecompositionGivens(Matrix A) {
	int n = A.size();

	Matrix Q(n);
	Matrix R(A);

	for (int i = 0; i < n; i++)
		Q(i, i) = 1;

	for (int j = 0; j < n; j++) {
		for (int i = n - 1; i > j; i--) {
			Matrix G(n);

			for (int k = 0; k < n; k++)
				G(k, k) = 1;

			double a = R(i - 1, j);
			double b = R(i, j);

			double c;
			double s;

			if (fabs(b) == 0 && fabs(a) == 0) {
				c = 1;
				s = 0;
			}
			else if (fabs(b) > fabs(a)) {
				double r = a / b;
				s = 1 / sqrt(1 + r * r);
				c = s * r;
			}
			else {
				double r = b / a;
				c = 1 / sqrt(1 + r * r);
				s = c * r;
			}

			G(i - 1, i - 1) = c;
			G(i - 1, i) = -s;
			G(i, i - 1) = s;
			G(i, i) = c;

			R = G.Transpose() * R;
			Q = Q * G;
		}
	}

	// для единственности QR разложения требуется положительность диагонали R
	for (int i = 0; i < n; i++) {
		if (R(i, i) < 0) {
			for (int j = 0; j < n; j++) {
				Q(j, i) = -Q(j, i);
				R(i, j) = -R(i, j);
			}
		}
	}

	cout << "QR decomposition (Givens):" << endl;
	cout << "Matrix Q:" << endl << Q << endl;
	cout << "Matrix R:" << endl << R << endl;
	cout << "Q * R:" << endl << (Q * R);
	cout << endl;
}

// QR разложение методом отражений (Хаусхолдера)
// A — квадратная невырожденная матрица
// Q — унитарная матрица
// R — верхнетреугольная матрица
// A = Q*R
void QRDecompositionHouseholder(Matrix A) {
	int n = A.size();

	Matrix Q(n);
	Matrix R(A);

	for (int i = 0; i < n; i++)
		Q(i, i) = 1;

	for (int k = 0; k < n; k++) {
		double *u = new double[n - k];
		double norm = 0;

		for (int i = k; i < n; i++) {
			u[i - k] = R(i, k);
			norm += R(i, k) * R(i, k);
		}

		double b = u[0] > 0 ? sqrt(norm) : -sqrt(norm);

		u[0] += b;

		Matrix H(n);

		for (int i = 0; i < n; i++)
			H(i, i) = 1;

		for (int i = 0; i < n - k; i++)
			for (int j = 0; j < n - k; j++)
				H(i + k, j + k) -= u[i] * u[j] / (R(k, k) * b + norm);

		R = H * R;
		Q = Q * H;

		delete[] u;
	}

	// для единственности QR разложения требуется положительность диагонали R
	for (int i = 0; i < n; i++) {
		if (R(i, i) < 0) {
			for (int j = 0; j < n; j++) {
				Q(j, i) = -Q(j, i);
				R(i, j) = -R(i, j);
			}
		}
	}

	cout << "QR decomposition (Householder):" << endl;
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

// LR разложение
// A — квадратная невырожденная матрица
// L — единичная нижняя треугольная матрица
// R — верхняя треугольная
// A = L*R
void LRDecomposition(Matrix A) {
	int n = A.size();

	Matrix L(n);
	Matrix R(n);

	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			R(i, j) = A(i, j);
			L(j, i) = A(j, i);

			for (int k = 0; k < i; k++) {
				R(i, j) -= L(i, k) * R(k, j);
				L(j, i) -= L(j, k) * R(k, i);
			}
			
			L(j, i) /= R(i, i);
		}
	}

	cout << "LR decomposition:" << endl;
	cout << "Matrix L:" << endl << L << endl;
	cout << "Matrix R:" << endl << R << endl;
	cout << "L * R:" << endl << (L * R);
	cout << endl;
}

int main() {
	Matrix A = ReadMatrix(); // считываем матрицу

	cout << endl << "Entered matrix: " << endl << A << endl; // выводим введённую матрицу

	LUDecomposition(A); // находим LU разложение
	CholeskyDecomposition(A); // находим разложение Холецкого
	QRDecomposition(A); // находим QR разложение
	QRDecompositionGivens(A); // находим QR разложение методом вращений (Гивенса)
	QRDecompositionHouseholder(A); // находим QR разложение методом отражений (Хаусхолдера)
	LDLDecomposition(A); // находим LDL разложение
	LRDecomposition(A); // находим LR разложение
}