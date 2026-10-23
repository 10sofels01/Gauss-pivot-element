#ifndef GAUSS_SOLVER_H_INCLUDED
#define GAUSS_SOLVER_H_INCLUDED

// Решение системы Ax = b методом Гаусса с выбором главного элемента по всей матрице.
// A и b модифицируются. Решение помещается в x.
// Возвращает true при успехе, false если матрица вырождена.
bool gauss_solve(int n, double* A, double* b, double* x);


#endif // GAUSS_SOLVER_H_INCLUDED
