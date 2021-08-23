#include <iostream>
#include <math.h>
#include <chrono>
#include <vector>
#include <memory>
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
        checkArguments(outputs, inputs);

        double E = inputs[0][0]; // Young's Modulus Pa
        double A = inputs[1][0]; // Area of bar mm2

        // Change inputs into vectors, which are easier to manipulate in C++ and can link with sparse solvers
        // and matrix classes. Values are stored in ROW-MAJOR order in the vectors (A[i][j] = v[i*columns +j])
        matlab::data::TypedArray<double> coords_Array = std::move(inputs[2]);
        std::vector<double> node_coords(coords_Array.begin(), coords_Array.end()); // Node coordinates
        matlab::data::TypedArray<uint64_t> els_Array = std::move(inputs[3]);
        std::vector<uint64_t> ELEMENTS(els_Array.begin(), els_Array.end()); // Nodal connectivity

        // Initialise sparse K matrix vectors to be passed by reference to the function that assembles K
        std::vector<uint64_t> col_index, row_position;
        std::vector<double> values;

        // Call function to assemble K
        Truss_K(E, A, node_coords, ELEMENTS, col_index, row_position, values);

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

        // Build sparse matrix object for output to MATLAB. This is taken from an example in the MATHWORKS documentation
        size_t nnz = values.size();

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

    // CST element stiffness matrix
    std::vector<double> k_Truss(double E, double A, double x1, double x2, double y1, double y2)
    {
        std::vector<double> k(16, 0);

        double alpha = atan2(y2 - y1, x2 - x1); // angle of inclination relative to the POSITIVE direction of the x axis

        // angle parameters
        double c = cos(alpha);
        double s = sin(alpha);

        double Le = sqrt(pow((x2 - x1), 2.0) + pow((y2 - y1), 2.0)); // element length
        double ke = E * A / Le;                                      // element axial stiffness

        //Row 1 - element dof11
        k[0] = ke * c * c;  //Col 1 - element dof11
        k[1] = ke * c * s;  //Col 2 - element dof12
        k[2] = -ke * c * c; //Col 3 - element dof21
        k[3] = -ke * c * s; //Col 4 - element dof11

        //Row 2 - element dof12
        k[1 * 4] = ke * c * s;      //Col 1 - element dof11
        k[1 * 4 + 1] = ke * s * s;  //Col 2 - element dof12
        k[1 * 4 + 2] = -ke * c * s; //Col 3 - element dof21
        k[1 * 4 + 3] = -ke * s * s; //Col 4 - element dof11

        //Row 3 - element dof21
        k[2 * 4] = -ke * c * c;     //Col 1 - element dof11
        k[2 * 4 + 1] = -ke * c * s; //Col 2 - element dof12
        k[2 * 4 + 2] = ke * c * c;  //Col 3 - element dof21
        k[2 * 4 + 3] = ke * c * s;  //Col 4 - element dof11

        //Row 4 - element dof22
        k[3 * 4] = -ke * c * s;     //Col 1 - element dof11
        k[3 * 4 + 1] = -ke * s * s; //Col 2 - element dof12
        k[3 * 4 + 2] = ke * c * s;  //Col 3 - element dof21
        k[3 * 4 + 3] = ke * s * s;  //Col 4 - element dof11

        return k;
    }

    // This function assembles the global stiffness matrix for CST elements
    void Truss_K(double E, double A, std::vector<double> node_coords, std::vector<uint64_t> ELEMENTS, std::vector<uint64_t> &col_index, std::vector<uint64_t> &row_position, std::vector<double> &values)
    {
        uint64_t els = ELEMENTS.size() / 2;          // Number of elements
        uint64_t nodes = node_coords.size() / 2;     // Number of nodes
        uint64_t n_dofs = node_coords.size();        // Number of DOFs
        double x1, x2, y1, y2;                       // Element nodal coordinates
        uint64_t n1, n2, dof11, dof12, dof21, dof22; // Element DOFs

        // Create vector of DOFs
        std::vector<uint64_t> dofs(n_dofs);
        std::iota(std::begin(dofs), std::end(dofs), 0);

        // Initialise vector of tuples containing row and column indices of NNZs
        std::vector<std::tuple<uint64_t, uint64_t>> v;

        // Determine sparsity pattern
        for (uint64_t i = 0; i < els; i++)
        {
            // Identify element node numbers. SUBTRACT 1 DUE TO MATLAB INDEXING STARTING AT 1
            n1 = ELEMENTS[i * 2] - 1;
            n2 = ELEMENTS[i * 2 + 1] - 1;

            // element node 1 - dofs
            dof11 = dofs[n1 * 2];
            dof12 = dofs[n1 * 2 + 1];
            // element node 2 - dofs
            dof21 = dofs[n2 * 2];
            dof22 = dofs[n2 * 2 + 1];

            std::vector<uint64_t> e_dofs{};
            e_dofs.push_back(dof11);
            e_dofs.push_back(dof12);
            e_dofs.push_back(dof21);
            e_dofs.push_back(dof22);

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
            n1 = ELEMENTS[i * 2] - 1;
            n2 = ELEMENTS[i * 2 + 1] - 1;

            // element node 1 - x,y coordinates
            x1 = node_coords[n1 * 2];
            y1 = node_coords[n1 * 2 + 1];
            // element node 2 - x,y coordinates
            x2 = node_coords[n2 * 2];
            y2 = node_coords[n2 * 2 + 1];

            // element node 1 - dofs
            dof11 = dofs[n1 * 2];
            dof12 = dofs[n1 * 2 + 1];
            // element node 2 - dofs
            dof21 = dofs[n2 * 2];
            dof22 = dofs[n2 * 2 + 1];

            std::vector<uint64_t> e_dofs{};
            e_dofs.push_back(dof11);
            e_dofs.push_back(dof12);
            e_dofs.push_back(dof21);
            e_dofs.push_back(dof22);

            // Call function for element stiffness matrix
            std::vector<double> ke = k_Truss(E, A, x1, x2, y1, y2);
            uint64_t r, r_next;

            for (uint64_t j = 0; j < 4; j++)
            {
                r = row_position[e_dofs[j]];
                r_next = row_position[e_dofs[j] + 1];
                for (uint64_t k = 0; k < 4; k++)
                {
                    for (uint64_t c = r; c < r_next; c++)
                    {
                        // If the column index matches with the DOF, add contribution to global K
                        if (col_index[c] == e_dofs[k])
                        {
                            values[c] += ke[j * 4 + k];
                            break;
                        }
                    }
                }
            }
        }
    }

    void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
    {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        matlab::data::ArrayFactory factory;

        if (inputs.size() != 4)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("4 inputs required")}));
        }

        if (inputs[0].getNumberOfElements() != 1)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("1st input must be a scalar")}));
        }

        if (inputs[0].getType() != matlab::data::ArrayType::DOUBLE ||
            inputs[0].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("1st input must be a noncomplex scalar double")}));
        }
        if (inputs[1].getNumberOfElements() != 1)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("2nd input must be a scalar")}));
        }

        if (inputs[1].getType() != matlab::data::ArrayType::DOUBLE ||
            inputs[1].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("2nd input must be a noncomplex scalar double")}));
        }

        if (inputs[2].getType() != matlab::data::ArrayType::DOUBLE ||
            inputs[2].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("3rd input matrix must be type double")}));
        }

        if (inputs[2].getDimensions()[0] != 1 && inputs[2].getDimensions()[1] != 1)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("3rd input must be one-dimensional matrix")}));
        }

        if (inputs[3].getType() != matlab::data::ArrayType::UINT64 ||
            inputs[3].getType() == matlab::data::ArrayType::COMPLEX_UINT64)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("4th input matrix must be type uint64")}));
        }

        if (inputs[3].getDimensions()[0] != 1 && inputs[3].getDimensions()[1] != 1)
        {
            matlabPtr->feval(u"error",
                             0, std::vector<matlab::data::Array>({factory.createScalar("4th input must be one-dimensional matrix")}));
        }
    }
};