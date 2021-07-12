#include <iostream>
#include <math.h>
#include <chrono>
#include <vector>
#include <memory>
// #include "Matrix.h"
// #include "Matrix.cpp"
// #include "CSRMatrix.h"
// #include "CSRMatrix.cpp"
// #include "SparseSolver.h"
// #include "SparseSolver.cpp"
// #include "utilities.h"
#include <memory>
#include <tuple>
#include <algorithm>
#include <numeric>
#include "mex.hpp"        // for the C++ MEX API.
#include "mexAdapter.hpp" // for the implementation of MexFunction class
#include "MatlabDataArray.hpp"

class MexFunction : public matlab::mex::Function
{
public:
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {
        // checkArguments(outputs, inputs);
        double E = inputs[0][0];  // Young's Modulus,
        double nu = inputs[1][0]; // Poisson coefficient
        double t = inputs[2][0];  // Thickness of shell

        // Change inputs into vectors, which are easier to manipulate in C++ and can link with sparse solvers
        // and matrix classes. Values are stored in ROW-MAJOR order in the vectors (A[i][j] = v[i*columns +j])
        matlab::data::TypedArray<double> coords_Array = std::move(inputs[3]);
        std::vector<double> node_coords(coords_Array.begin(), coords_Array.end()); // Node coordinates
        matlab::data::TypedArray<uint64_t> els_Array = std::move(inputs[4]);
        std::vector<uint64_t> ELEMENTS(els_Array.begin(), els_Array.end()); // Nodal connectivity

        // Initialise sparse K matrix vectors to be passed by reference to the function that assembles K
        std::vector<uint64_t> col_index, row_position;
        std::vector<double> values;

        // Call function to assemble K
        CST_K(nu, E, t, node_coords, ELEMENTS, col_index, row_position, values);

        // Change row_position vector to include a row index for each NNZ value
        uint64_t r, r_next, row_nnzs;
        std::vector<uint64_t> rows;
        for (uint64_t i = 0; i < row_position.size() - 1; i++)
        {
            r = row_position[i];
            r_next = row_position[i + 1];
            row_nnzs = r_next - r;

            for (uint64_t n = 0; n < row_nnzs; n++)
            {
                rows.push_back(i);
            }
        }

        // Build sparse matrix object for output to MATLAB. This is taken from an example in the MATHWROKS documentation
        size_t nnz = values.size();

        std::cout << "Size: " << values.size() << std::endl;

        matlab::data::ArrayFactory factory;
        auto values_p = factory.createBuffer<double>(nnz);
        auto rows_p = factory.createBuffer<size_t>(nnz);
        auto cols_p = factory.createBuffer<size_t>(nnz);

        double *valuesPtr = values_p.get();
        size_t *rowsPtr = rows_p.get();
        size_t *colsPtr = cols_p.get();
        std::for_each(values.begin(), values.end(), [&](const double &e)
                      { *(valuesPtr++) = e; });
        std::for_each(rows.begin(), rows.end(), [&](const size_t &e)
                      { *(rowsPtr++) = e; });
        std::for_each(col_index.begin(), col_index.end(), [&](const size_t &e)
                      { *(colsPtr++) = e; });

        // Assign MEX function output
        outputs[0] = factory.createSparseArray<double>({node_coords.size(), node_coords.size()}, nnz, std::move(values_p), std::move(rows_p), std::move(cols_p));
    }

    void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;

        if (inputs.size() != 2)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("Two inputs required")}));
        }

        if (inputs[0].getNumberOfElements() != 1)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("Input multiplier must be a scalar")}));
        }

        if (inputs[0].getType() != matlab::data::ArrayType::DOUBLE ||
            inputs[0].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("Input multiplier must be a noncomplex scalar double")}));
        }

        if (inputs[3].getType() != matlab::data::ArrayType::DOUBLE ||
            inputs[3].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("Input matrix must be type double")}));
        }

        if (inputs[3].getDimensions().size() != 2)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("Input must be m-by-n dimension")}));
        }
    }

    // CST element stiffness matrix
    std::vector<double> k_CST(double nu, double E, double t, double x1, double x2, double x3, double y1, double y2, double y3)
    {
        std::vector<double> k(36, 0);
        // The following block of code has been generated by Maple's CodeGeneration
        k[0] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * (pow(y2 - y3, 0.2e1) + pow(-x2 + x3, 0.2e1) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1)) / 0.2e1;
        k[1] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((y2 - y3) * (double)nu * (-x2 + x3) + (-x2 + x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (y2 - y3)) / 0.2e1;
        k[2] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((y2 - y3) * (-y1 + y3) + (-x2 + x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (x1 - x3)) / 0.2e1;
        k[3] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((y2 - y3) * (double)nu * (x1 - x3) + (-x2 + x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-y1 + y3)) / 0.2e1;
        k[4] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((y2 - y3) * (y1 - y2) + (-x2 + x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-x1 + x2)) / 0.2e1;
        k[5] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((y2 - y3) * (double)nu * (-x1 + x2) + (-x2 + x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (y1 - y2)) / 0.2e1;
        k[1 * 6] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((y2 - y3) * (double)nu * (-x2 + x3) + (-x2 + x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (y2 - y3)) / 0.2e1;
        k[1 * 6 + 1] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * (pow(-x2 + x3, 0.2e1) + pow(y2 - y3, 0.2e1) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1)) / 0.2e1;
        k[1 * 6 + 2] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-x2 + x3) * (double)nu * (-y1 + y3) + (y2 - y3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (x1 - x3)) / 0.2e1;
        k[1 * 6 + 3] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-x2 + x3) * (x1 - x3) + (y2 - y3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-y1 + y3)) / 0.2e1;
        k[1 * 6 + 4] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-x2 + x3) * (double)nu * (y1 - y2) + (y2 - y3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-x1 + x2)) / 0.2e1;
        k[1 * 6 + 5] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-x2 + x3) * (-x1 + x2) + (y2 - y3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (y1 - y2)) / 0.2e1;
        k[2 * 6] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((y2 - y3) * (-y1 + y3) + (-x2 + x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (x1 - x3)) / 0.2e1;
        k[2 * 6 + 1] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-x2 + x3) * (double)nu * (-y1 + y3) + (y2 - y3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (x1 - x3)) / 0.2e1;
        k[2 * 6 + 2] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * (pow(-y1 + y3, 0.2e1) + pow(x1 - x3, 0.2e1) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1)) / 0.2e1;
        k[2 * 6 + 3] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-y1 + y3) * (double)nu * (x1 - x3) + (x1 - x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-y1 + y3)) / 0.2e1;
        k[2 * 6 + 4] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-y1 + y3) * (y1 - y2) + (x1 - x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-x1 + x2)) / 0.2e1;
        k[2 * 6 + 5] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-y1 + y3) * (double)nu * (-x1 + x2) + (x1 - x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (y1 - y2)) / 0.2e1;
        k[3 * 6] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((y2 - y3) * (double)nu * (x1 - x3) + (-x2 + x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-y1 + y3)) / 0.2e1;
        k[3 * 6 + 1] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-x2 + x3) * (x1 - x3) + (y2 - y3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-y1 + y3)) / 0.2e1;
        k[3 * 6 + 2] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-y1 + y3) * (double)nu * (x1 - x3) + (x1 - x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-y1 + y3)) / 0.2e1;
        k[3 * 6 + 3] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * (pow(x1 - x3, 0.2e1) + pow(-y1 + y3, 0.2e1) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1)) / 0.2e1;
        k[3 * 6 + 4] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((x1 - x3) * (double)nu * (y1 - y2) + (-y1 + y3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-x1 + x2)) / 0.2e1;
        k[3 * 6 + 5] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((x1 - x3) * (-x1 + x2) + (-y1 + y3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (y1 - y2)) / 0.2e1;
        k[4 * 6] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((y2 - y3) * (y1 - y2) + (-x2 + x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-x1 + x2)) / 0.2e1;
        k[4 * 6 + 1] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-x2 + x3) * (double)nu * (y1 - y2) + (y2 - y3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-x1 + x2)) / 0.2e1;
        k[4 * 6 + 2] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-y1 + y3) * (y1 - y2) + (x1 - x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-x1 + x2)) / 0.2e1;
        k[4 * 6 + 3] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((x1 - x3) * (double)nu * (y1 - y2) + (-y1 + y3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (-x1 + x2)) / 0.2e1;
        k[4 * 6 + 4] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * (pow(y1 - y2, 0.2e1) + pow(-x1 + x2, 0.2e1) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1)) / 0.2e1;
        k[4 * 6 + 5] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((y1 - y2) * (double)nu * (-x1 + x2) + (-x1 + x2) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (y1 - y2)) / 0.2e1;
        k[5 * 6] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((y2 - y3) * (double)nu * (-x1 + x2) + (-x2 + x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (y1 - y2)) / 0.2e1;
        k[5 * 6 + 1] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-x2 + x3) * (-x1 + x2) + (y2 - y3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (y1 - y2)) / 0.2e1;
        k[5 * 6 + 2] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((-y1 + y3) * (double)nu * (-x1 + x2) + (x1 - x3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (y1 - y2)) / 0.2e1;
        k[5 * 6 + 3] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((x1 - x3) * (-x1 + x2) + (-y1 + y3) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (y1 - y2)) / 0.2e1;
        k[5 * 6 + 4] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * ((y1 - y2) * (double)nu * (-x1 + x2) + (-x1 + x2) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1) * (y1 - y2)) / 0.2e1;
        k[5 * 6 + 5] = E * t / (double)(-nu * nu + 1) / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2) * (pow(-x1 + x2, 0.2e1) + pow(y1 - y2, 0.2e1) * (0.1e1 / 0.2e1 - (double)nu / 0.2e1)) / 0.2e1;

        return k;
    }

    // This function assembles the global stiffness matrix for CST elements
    void CST_K(double nu, double E, double t, std::vector<double> node_coords, std::vector<uint64_t> ELEMENTS, std::vector<uint64_t> &col_index, std::vector<uint64_t> &row_position, std::vector<double> &values)
    {
        uint64_t els = ELEMENTS.size() / 3;                            // Number of elements
        uint64_t nodes = node_coords.size() / 2;                       // Number of nodes
        uint64_t n_dofs = node_coords.size();                          // Number of DOFs
        double x1, x2, x3, y1, y2, y3;                                 // Element nodal coordinates
        uint64_t n1, n2, n3, dof11, dof12, dof21, dof22, dof31, dof32; // Element DOFs

        // Create vector of DOFs
        std::vector<uint64_t> dofs(n_dofs);
        std::iota(std::begin(dofs), std::end(dofs), 0);

        // Initialise vector of tuples containing row and column indices of NNZs
        std::vector<std::tuple<uint64_t, uint64_t>> v;

        // Determine sparsity pattern
        for (uint64_t i = 0; i < els; i++)
        {
            // Identify element node numbers. SUBTRACT 1 DUE TO MATLAB INDEXING STARTING AT 1
            n1 = ELEMENTS[i * 3] - 1;
            n2 = ELEMENTS[i * 3 + 1] - 1;
            n3 = ELEMENTS[i * 3 + 2] - 1;

            // element node 1 - dofs
            dof11 = dofs[n1 * 2];
            dof12 = dofs[n1 * 2 + 1];
            // element node 2 - dofs
            dof21 = dofs[n2 * 2];
            dof22 = dofs[n2 * 2 + 1];
            // element node 3 - dofs
            dof31 = dofs[n3 * 2];
            dof32 = dofs[n3 * 2 + 1];

            std::vector<uint64_t> e_dofs{};
            e_dofs.push_back(dof11);
            e_dofs.push_back(dof12);
            e_dofs.push_back(dof21);
            e_dofs.push_back(dof22);
            e_dofs.push_back(dof31);
            e_dofs.push_back(dof32);

            // Register all NNZ values. SOme will be duplicate
            for (uint64_t j = 0; j < e_dofs.size(); j++)
            {
                for (uint64_t k = 0; k < e_dofs.size(); k++)
                {
                    v.push_back(std::make_tuple(e_dofs[j], e_dofs[k]));
                }
            }
        }

        // Sort vector of tuples by row indices and remove duplicates
        sort(v.begin(), v.end());
        v.erase(unique(v.begin(), v.end()), v.end());

        values.resize(v.size(), 0);         // Non-zero values of stiffness matrix
        row_position.resize(n_dofs + 1, 0); // Indices in nnz values of stiffness matrix at start of each row

        // Construct row_position and col_index vectors
        uint64_t count = 0; // Counter for number of nnzs
        for (uint64_t i = 0; i < n_dofs; i++)
        {
            // Vector containing column indices for row i
            std::vector<uint64_t> col_nnzs{};
            while (count <= v.size())
            {
                if (std::get<0>(v[count]) == i)
                {
                    // Append column indices for row i
                    col_nnzs.push_back(std::get<1>(v[count]));
                    count++;
                }
                else if (std::get<0>(v[count]) != i || count == v.size())
                {
                    // Sort column indices
                    sort(col_nnzs.begin(), col_nnzs.end());
                    for (uint64_t k = 0; k < col_nnzs.size(); k++)
                    {
                        // Append column indices to full column index vector
                        col_index.push_back(col_nnzs[k]);
                    }
                    if (count == v.size())
                    {
                        count++;
                    }
                    // When all the NNZs in the row are registered, move on to the next row
                    break;
                }
            }
            row_position[i + 1] = count;
        }

        // Calculate global stiffness matrix coefficients with known sparsity pattern
        for (uint64_t i = 0; i < els; i++)
        {
            // Identify element node numbers
            n1 = ELEMENTS[i * 3] - 1;
            n2 = ELEMENTS[i * 3 + 1] - 1;
            n3 = ELEMENTS[i * 3 + 2] - 1; // element node 3 - x,y coordinates

            // element node 1 - x,y coordinates
            x1 = node_coords[n1 * 2];
            y1 = node_coords[n1 * 2 + 1];
            // element node 2 - x,y coordinates
            x2 = node_coords[n2 * 2];
            y2 = node_coords[n2 * 2 + 1];
            // element node 3 - x,y coordinates
            x3 = node_coords[n3 * 2];
            y3 = node_coords[n3 * 2 + 1];

            // element node 1 - dofs
            dof11 = dofs[n1 * 2];
            dof12 = dofs[n1 * 2 + 1];
            // element node 2 - dofs
            dof21 = dofs[n2 * 2];
            dof22 = dofs[n2 * 2 + 1];
            // element node 3 - dofs
            dof31 = dofs[n3 * 2];
            dof32 = dofs[n3 * 2 + 1];

            std::vector<uint64_t> e_dofs{};
            e_dofs.push_back(dof11);
            e_dofs.push_back(dof12);
            e_dofs.push_back(dof21);
            e_dofs.push_back(dof22);
            e_dofs.push_back(dof31);
            e_dofs.push_back(dof32);

            // Call function for element stiffness matrix
            std::vector<double> ke = k_CST(nu, E, t, x1, x2, x3, y1, y2, y3);
            uint64_t r, r_next;

            for (uint64_t j = 0; j < 6; j++)
            {
                r = row_position[e_dofs[j]];
                r_next = row_position[e_dofs[j] + 1];
                for (uint64_t k = 0; k < 6; k++)
                {
                    for (uint64_t c = r; c < r_next; c++)
                    {
                        // If the column index matches with the DOF, add contribution to global K
                        if (col_index[c] == e_dofs[k])
                        {
                            values[c] += ke[j * 6 + k];
                            break;
                        }
                    }
                }
            }
        }
    };
};