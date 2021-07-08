#pragma once
#include <iostream>
#include <vector>
#include <memory>

template <class T>
class Matrix
{
public:
    // default constructor
    Matrix(void);

    // constructor where we want to preallocate ourselves
    Matrix(int rows, int cols, bool preallocate);

    // constructor where we already have allocated memory outside
    Matrix(int rows, int cols, std::shared_ptr<T[]> values_ptr);

    // Copy constructor
    Matrix(const Matrix<T> &M2);

    // Overload assignment operator to deepcopy
    Matrix<T> &operator=(const Matrix<T> &M2);

    // destructor
    virtual ~Matrix();

    // Print out the values in our matrix
    void printValues();
    virtual void printMatrix();

    void matMatMult(Matrix<T> &mat_right, Matrix<T> &output);

    void matVecMult(std::vector<T> &vec, std::vector<T> &output);

    std::shared_ptr<T[]> values;
    int rows = -1;
    int cols = -1;

    int size_of_values = -1;
    bool preallocated = false;
};