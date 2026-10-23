#include "gauss_solver.h"
#include <cmath>
#include <limits>
#include <algorithm>

using namespace std;

bool gauss_solve(int n, double* A, double* b, double* x) {

    // Находим максимальный по модулю элемент во всей матрице.
    // Он нужен, чтобы задать порог tol, ниже которого элемент считается нулём.
    double max_elem = 0.0;
    for (int i = 0; i < n; ++i) {
        const double* Ai = A + i * n;
        for (int j = 0; j < n; ++j) {
            double v = Ai[j];
            if (v < 0) v = -v;
            if (v > max_elem) max_elem = v;
        }
    }
    if (max_elem < 1e-15) {
        return false;
    }

    double eps = numeric_limits<double>::epsilon();
    double tol = eps * n * max_elem;

    // Массив перестановок столбцов, который после всех перестановок оказался на позиции i.
    int* col_perm = new int[n];
    for (int i = 0; i < n; ++i) col_perm[i] = i;   // изначально порядок не тронут

    // Кэш построчных максимумов. Для каждой строки i храним:
    //   row_max[i] — максимальный по модулю элемент среди ещё не обработанных столбцов
    //   row_arg[i] — номер столбца, где этот максимум достигается
    double* row_max = new double[n];
    int*    row_arg = new int[n];

    // Временный массив для обратного хода
    double* y = new double[n];

    // Для каждой строки ищем максимум по всем столбцам
    for (int i = 0; i < n; ++i) {
        const double* Ai = A + i * n;
        double mx = 0.0;
        int arg = 0;
        for (int j = 0; j < n; ++j) {
            double v = Ai[j];
            if (v < 0) v = -v;
            if (v > mx) { mx = v; arg = j; }
        }
        row_max[i] = mx;
        row_arg[i] = arg;
    }

    // Прямой ход: обнуляем поддиагональные элементы, шаг за шагом
    for (int k = 0; k < n - 1; ++k) {

        // --- Поиск главного элемента по всей подматрице A[i][j] (i>=k, j>=k) ---

        double best = -1.0;
        int pivot_row = k;
        for (int i = k; i < n; ++i) {
            if (row_max[i] > best) {
                best = row_max[i];
                pivot_row = i;
            }
        }

        // 2) Теперь в найденной строке ищем столбец с максимумом по j >= k.
        double max_val = 0.0;
        int pivot_col = k;
        {
            const double* Ap = A + pivot_row * n;
            for (int j = k; j < n; ++j) {
                double v = Ap[j];
                if (v < 0) v = -v;
                if (v > max_val) { max_val = v; pivot_col = j; }
            }
        }

        // Если максимальный элемент подматрицы меньше tol — матрица вырождена
        if (max_val < tol) {
            delete[] col_perm;
            delete[] row_max;
            delete[] row_arg;
            delete[] y;
            return false;
        }

        // --- Перестановка строк: ставим строку с главным элементом на место k ---
        if (pivot_row != k) {
            double* Ak = A + k * n;
            double* Ap = A + pivot_row * n;
            for (int j = 0; j < n; ++j) {           // меняем местами две строки поэлементно
                double t = Ak[j]; Ak[j] = Ap[j]; Ap[j] = t;
            }
            double tb = b[k]; b[k] = b[pivot_row]; b[pivot_row] = tb;   // и правую часть
            // Кэш построчных максимумов едет вместе со строками
            double tm = row_max[k]; row_max[k] = row_max[pivot_row]; row_max[pivot_row] = tm;
            int    ta = row_arg[k]; row_arg[k] = row_arg[pivot_row]; row_arg[pivot_row] = ta;
        }

        // --- Перестановка столбцов: ставим столбец с главным элементом на место k ---
        if (pivot_col != k) {
            for (int i = 0; i < n; ++i) {           // меняем местами два столбца поэлементно
                double* Ai = A + i * n;
                double t = Ai[k]; Ai[k] = Ai[pivot_col]; Ai[pivot_col] = t;
            }
            // Запоминаем перестановку столбцов для восстановления ответа в конце
            int tc = col_perm[k]; col_perm[k] = col_perm[pivot_col]; col_perm[pivot_col] = tc;

            // row_arg ссылается на номера столбцов, а мы только что два столбца
            // поменяли местами. Значит, все row_arg[i], равные k или pivot_col,
            // надо тоже поменять местами, иначе кэш будет указывать не туда.
            for (int i = 0; i < n; ++i) {
                if      (row_arg[i] == k)         row_arg[i] = pivot_col;
                else if (row_arg[i] == pivot_col) row_arg[i] = k;
            }
        }

        // --- Исключение ниже диагонали: обнуляем A[i][k] для i > k ---
        double* Ak = A + k * n;
        double akk = Ak[k];                         // ведущий элемент (диагональный после перестановок)
        double bk  = b[k];                          // соответствующий элемент правой части
        for (int i = k + 1; i < n; ++i) {
            double* Ai = A + i * n;
            double factor = Ai[k] / akk;            // множитель: aik/akk
            Ai[k] = 0.0;                            // явно зануляем (в цикле j начинаем с k+1)
            for (int j = k + 1; j < n; ++j) {       // j с k+1, т.к. столбец k уже обработан
                Ai[j] -= factor * Ak[j];            // aij -= factor * akj
            }
            b[i] -= factor * bk;                    // и правую часть тоже
        }


        // Пересчитываем строку только если её максимум стоял в столбце <= k, то есть в уже исключённом столбце.
        for (int i = k + 1; i < n; ++i) {
            if (row_arg[i] <= k) {
                const double* Ai = A + i * n;
                double mx = 0.0;
                int arg = k + 1;
                for (int j = k + 1; j < n; ++j) {
                    double v = Ai[j];
                    if (v < 0) v = -v;
                    if (v > mx) { mx = v; arg = j; }
                }
                row_max[i] = mx;
                row_arg[i] = arg;
            }
        }
    }

    // --- Обратный ход: находим y (решение в переставленном порядке столбцов) ---
    for (int i = n - 1; i >= 0; --i) {
        const double* Ai = A + i * n;
        double sum = b[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= Ai[j] * y[j];                    // вычитаем уже известные слагаемые
        }
        double diag = Ai[i];
        if (diag < 0) diag = -diag;                 // модуль диагонального элемента
        if (diag < tol) {
            // Диагональный элемент практически ноль — матрица вырождена
            delete[] col_perm;
            delete[] row_max;
            delete[] row_arg;
            delete[] y;
            return false;
        }
        y[i] = sum / Ai[i];                         // делим на диагональный элемент
    }

    // Восстановление исходного порядка переменных:
    // y[i] соответствует столбцу col_perm[i] исходной матрицы
    for (int i = 0; i < n; ++i) {
        x[col_perm[i]] = y[i];
    }

    delete[] col_perm;
    delete[] row_max;
    delete[] row_arg;
    delete[] y;

    return true;
}
