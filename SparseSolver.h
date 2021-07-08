
#pragma once
#include "CSRMatrix.h"
#include <vector>
#include <memory>

template <class T>
class SparseSolver
{
public:
    CSRMatrix<T> A;

    std::vector<T> b{};

    SparseSolver(CSRMatrix<T> &A, std::vector<T> &b);

    ~SparseSolver();

    void stationaryIterative(std::vector<T> &x, double &tol, int &it_max, bool isGaussSeidel);

    T residualCalc(std::vector<T> &x, std::vector<T> &output_b);

    void conjugateGradient(std::vector<T> &x, double &tol, int &it_max);

    std::shared_ptr<CSRMatrix<T>> lu_decomp();
    void lu_solve(CSRMatrix<T> &LU, std::vector<int> &piv, std::vector<T> &x);

    std::shared_ptr<CSRMatrix<T>> cholesky_decomp();
    void cholesky_solve(CSRMatrix<T> &R, std::vector<T> &x);
};