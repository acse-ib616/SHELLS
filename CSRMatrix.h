#pragma once
#include "Matrix.h"
#include <vector>
#include <memory>

template <class T>
class CSRMatrix : public Matrix<T>
{
public:
    // default constructor
    CSRMatrix(void);

    // Constructor
    CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
    CSRMatrix(int size, double sparsity);

    //Constructor
    CSRMatrix(int rows, int cols, int nnzs, std::shared_ptr<T[]> values_ptr, std::shared_ptr<int[]> row_pos, std::shared_ptr<int[]> col_ind);

    // Copy constructor
    CSRMatrix(const CSRMatrix<T>& M2);

    CSRMatrix<T> &operator=(const CSRMatrix<T> &M2);

    ~CSRMatrix();

    virtual void printMatrix();

    virtual void print2DMatrix();

    void matVecMult(std::vector<T> &input, std::vector<T> &output);

    std::shared_ptr<CSRMatrix<T>> matMatMult(CSRMatrix<T> &mat_right);
    std::shared_ptr<CSRMatrix<T>> matMatMultSymbolic(CSRMatrix<T> &mat_right);

    CSRMatrix<T> cholesky();
    std::shared_ptr<CSRMatrix<T>> transpose();

    std::shared_ptr<int[]> row_position; //create nullpointer
    std::shared_ptr<int[]> col_index;    // create nullpointer

    // number of non-zeros
    int nnzs = -1;

    // we're inheriting the values pointer so we don't have to include it here
};