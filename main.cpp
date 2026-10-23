#include <iostream>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdio>
#include "matrix_io.h"
#include "gauss_solver.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " n m k [filename]" << endl;
        return 1;
    }

    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    int k = atoi(argv[3]);

    if (n <= 0) {
        cerr << "Error: n must be positive." << n << endl;
        return 1;
    }
    if (m <= 0) {
        cerr << "Error: m must be positive." << m << endl;
        return 1;
    }

    if (k < 0 || k > 4) {
        cerr << "Error: k must be 0, 1, 2, 3 or 4." << k << endl;
        return 1;
    }

    string filename;
    if (k == 0) {
        if (argc < 5) {
            cerr << "Error: filename required when k=0" << endl;
            return 1;
        }
        filename = argv[4];
    }

    vector<double> A;
    vector<double> b;
    vector<double> x;
    try {
        A.resize(n * n);
        b.resize(n);
        x.resize(n);
    } catch (const bad_alloc& e) {
        cerr << "Error: cannot allocate memory for n = " << n << endl;
        return 1;
    }

    // Инициализация матрицы
    if (k == 0) {
        if (!read_matrix_file(filename, n, A.data())) {
            return 1;
        }
    } else {
        init_matrix_formula(n, k, A.data());
    }

    // Построение правой части
    build_rhs(n, A.data(), b.data());

    // Вывод исходной матрицы A (не более m строк/столбцов)
    cout << "Initial matrix A:" << endl;
    print_matrix(n, n, A.data(), m);

    // Вывод правой части b (как столбец, не более m элементов)
    cout << "Right-hand side b:" << endl;
    for (int i = 0; i < min(n, m); ++i) {
        printf(" %10.3e\n", b[i]);
    }

    // Решение системы
    auto start = chrono::high_resolution_clock::now(); // время до
    bool success = gauss_solve(n, A.data(), b.data(), x.data());
    auto end = chrono::high_resolution_clock::now(); // после
    chrono::duration<double> elapsed = end - start;

    if (!success) {
        cerr << "Error: matrix is singular or nearly singular." << endl;
        return 1;
    }

    // Вывод решения x как 1 x n матрицы
    cout << "Solution x:" << endl;
    print_matrix(1, n, x.data(), m);

    // Повторная инициализация A и b для вычисления невязки
    if (k == 0) {
        read_matrix_file(filename, n, A.data());
    } else {
        init_matrix_formula(n, k, A.data());
    }
    build_rhs(n, A.data(), b.data());

    // Вычисление невязки ||Ax - b|| / ||b||
    vector<double> Ax(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Ax[i] += A[i * n + j] * x[j];
        }
    }
    double norm_res = 0.0, norm_b = 0.0;
    for (int i = 0; i < n; ++i) {
        norm_res += (Ax[i] - b[i]) * (Ax[i] - b[i]);
        norm_b += b[i] * b[i];
    }
    norm_res = sqrt(norm_res);
    norm_b = sqrt(norm_b);
    double rel_res = (norm_b > 0) ? norm_res / norm_b : norm_res;
    printf("Relative residual norm ||Ax - b|| / ||b|| = %e\n", rel_res);

    // Вычисление нормы погрешности ||x - x_exact||, x_exact = (1,0,1,0,...)
    double err_norm = 0.0;
    for (int i = 0; i < n; ++i) {
        double exact = (i % 2 == 0) ? 1.0 : 0.0;
        err_norm += (x[i] - exact) * (x[i] - exact);
    }
    err_norm = sqrt(err_norm);
    printf("Error norm ||x - x_exact|| = %e\n", err_norm);

    // Вывод времени счёта
    printf("Elapsed time: %.2f seconds\n", elapsed.count());

    return 0;
}
