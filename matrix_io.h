#ifndef MATRIX_IO_H_INCLUDED
#define MATRIX_IO_H_INCLUDED

#include <string>

using namespace std;

// Инициализация матрицы A (n x n) по формуле k (i,j = 1..n)
void init_matrix_formula(int n, int k, double* A);

// Чтение матрицы из файла. Возвращает true при успехе.
bool read_matrix_file(const string& filename, int n, double* A);

// Построение правой части b из матрицы A
void build_rhs(int n, const double* A, double* b);

// Вывод матрицы (l x n) не более чем m строк и столбцов
void print_matrix(int l, int n, const double* mat, int m);

#endif // MATRIX_IO_H_INCLUDED
