#pragma once
#include <iostream>
#include <vector>
// This header file includes functions that do not fall
// under the scope of the Matrix or Solver classes

template <typename T>
void printVector(std::vector<T> vec)
{
    std::cout << "Vector is: \n";
    for (int i = 0; i < vec.size(); i++)
    {
        std::cout << vec[i] << " ";
    }
    std::cout << "\n";
}

template <typename T>
void checkDimensions(Matrix<T> &M1, std::vector<T> &vec)
{

    // Check if square matrix
    if (M1.cols != M1.rows)
    {
        throw std::invalid_argument("Only implemented for square matrix");
    }

    if (M1.cols != vec.size())
    {
        throw std::invalid_argument("Dimensions don't match");
    }
}

template <typename T>
double vecDotProduct(std::vector<T> v1, std::vector<T> v2)
{
    double result = 0;
    if (v1.size() == v2.size())
    {
        for (int i = 0; i < v1.size(); i++)
        {
            result += v1[i] * v2[i];
        }
    }
    return result;
}
