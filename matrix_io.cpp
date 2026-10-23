#include "matrix_io.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

void init_matrix_formula(int n, int k, double* A) {
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            double val = 0.0;
            switch (k) {
                case 1: val = n - max(i, j) + 1; break;
                case 2: val = max(i, j); break;
                case 3: val = abs(i - j); break;
                case 4: val = 1.0 / (i + j - 1); break;
                default: val = 0.0; break;
            }
            A[(i - 1) * n + (j - 1)] = val;
        }
    }
}

bool read_matrix_file(const string& filename, int n, double* A) {
    ifstream fin(filename);
    if (!fin) {
        cerr << "Error: cannot open file " << filename << endl;
        return false;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!(fin >> A[i * n + j])) {
                cerr << "Error: insufficient data or wrong format in file " << filename << endl;
                return false;
            }
        }
    }
    return true;
}

void build_rhs(int n, const double* A, double* b) {
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        const double* row = A + i * n;
        for (int j = 0; j < n; j += 2) {
            sum += row[j];
        }
        b[i] = sum;
    }
}

void print_matrix(int l, int n, const double* mat, int m) {
    int rows = min(l, m);
    int cols = min(n, m);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf(" %10.3e", mat[i * n + j]);
        }
        printf("\n");
    }
}
