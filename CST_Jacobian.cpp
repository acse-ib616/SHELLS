#include <iostream>
#include <math.h>
#include <chrono>
#include <vector>
#include <memory>
#include <memory>
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
        double Ly = inputs[0][0]; // Height of domain mm

        // Change inputs into vectors, which are easier to manipulate in C++ and can link with sparse solvers
        // and matrix classes. Values are stored in ROW-MAJOR order in the vectors (A[i][j] = v[i*columns +j])
        matlab::data::TypedArray<double> q_Array = std::move(inputs[1]);
        std::vector<double> q(q_Array.begin(), q_Array.end()); // UDL values
        matlab::data::TypedArray<double> coords_Array = std::move(inputs[2]);
        std::vector<double> node_coords(coords_Array.begin(), coords_Array.end()); // Node coordinates
        std::cout << node_coords.size() << std::endl;
        matlab::data::TypedArray<uint64_t> els_Array = std::move(inputs[3]);
        std::vector<uint64_t> ELEMENTS(els_Array.begin(), els_Array.end()); // Nodal connectivity

        // Initialise sparse Jacobian matrix vectors to be passed by reference to the function that assembles the Jacobian
        std::vector<uint64_t> col_index, row_position;
        std::vector<double> values;

        std::cout << "Hello mate" << std::endl;

        // Call function to assemble K
        Jacobian(Ly, q, node_coords, ELEMENTS, col_index, row_position, values);

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
        outputs[0] = factory.createSparseArray<double>({node_coords.size(), q.size()}, nnz, std::move(values_p), std::move(rows_p), std::move(cols_p));
    }

    // This function assembles the vector the constitute the Jacobian matrix in sparse CSR format
    void Jacobian(double Ly, std::vector<double> q, std::vector<double> node_coords, std::vector<uint64_t> ELEMENTS, std::vector<uint64_t> &col_index, std::vector<uint64_t> &row_position, std::vector<double> &values)
    {

        uint64_t els = ELEMENTS.size() / 3;      // Number of elements
        uint64_t nodes = node_coords.size() / 2; // Number of nodes
        uint64_t n_dofs = node_coords.size();    // Number of DOFs
        double x1, x3, y1, y3, x31;              // Element nodal coordinates and element side (nodes 1 to 3)
        uint64_t n1, n3, dof12, dof32;           // Element DOFs
        uint64_t counter = 0;                    // Initialise UDL counter. Max. value will be length of q-1

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
            n3 = ELEMENTS[i * 3 + 2] - 1;

            // element node 1 - y coordinate
            y1 = node_coords[n1 * 2 + 1];
            // element node 3 - y coordinate
            y3 = node_coords[n3 * 2 + 1];

            dof12 = dofs[n1 * 2 + 1]; // element node 1 - dofs
            dof32 = dofs[n3 * 2 + 1];

            // If element is at the top edge of planar body
            if (y1 == Ly && y3 == Ly)
            {
                // Register all NNZ values. Some will be duplicate
                v.push_back(std::make_tuple(dof12, counter));
                v.push_back(std::make_tuple(dof32, counter));
                counter++;
            }
        }

        // Sort vector and remove duplicates
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

        uint64_t r, r_next;
        counter = 0;
        // Calculate global stiffness matrix coefficients with known sparsity pattern
        for (uint64_t i = 0; i < els; i++)
        {
            // Identify element node numbers
            n1 = ELEMENTS[i * 3] - 1;
            n3 = ELEMENTS[i * 3 + 2] - 1;

            // element node 1 - x,y coordinates
            x1 = node_coords[n1 * 2];
            y1 = node_coords[n1 * 2 + 1];
            // element node 3 - y coordinate
            x3 = node_coords[n3 * 2];
            y3 = node_coords[n3 * 2 + 1];

            // If element is at the top edge of planar body
            if (y1 == Ly && y3 == Ly)
            {
                x31 = x3 - x1;

                // element node 1 - dofs
                dof12 = dofs[n1 * 2 + 1];
                // element node 3 - dofs
                dof32 = dofs[n3 * 2 + 1];
                std::vector<uint64_t> e_dofs{};
                e_dofs.push_back(dof12);
                e_dofs.push_back(dof32);

                for (uint64_t j = 0; j < 2; j++)
                {
                    r = row_position[e_dofs[j]];
                    r_next = row_position[e_dofs[j] + 1];
                    for (uint64_t c = r; c < r_next; c++)
                    {
                        // If the column index matches with the UDL counter, add contribution to Jacobian
                        if (col_index[c] == counter)
                        {
                            values[c] += 0.5 * abs(x31) * 1e-3;
                            break;
                        }
                    }
                }
                counter++;
            }
        }
    };
};