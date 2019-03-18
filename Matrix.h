#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>

// класс матрицы
class Matrix {
	double **values;
	int n;

public:
	Matrix(int n); // конструктор матрицы из размера

	Matrix(const Matrix& matrix); // конструктор копирования

	Matrix& operator=(const Matrix& matrix); // оператор присваивания

	int size() const; // получение размера строк

	double operator()(int i, int j) const; // получение элемента по индексу
	double& operator()(int i, int j); // получение элемента по индексу

	Matrix operator*(const Matrix &matrix) const; // умножение матрицы на матрицу
	
	void Read(); // заполнение матрицы с клавиатуры
	
	Matrix Transpose() const; // получение транспонированной матрицы

	~Matrix(); // деструктор (освобождение памяти)

	friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix); // вывод в поток
};