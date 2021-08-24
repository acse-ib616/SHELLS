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
        matlab::data::TypedArray<uint32_t> els_Array = std::move(inputs[3]);
        std::vector<uint32_t> ELEMENTS(els_Array.begin(), els_Array.end()); // Nodal connectivity

        // Initialise sparse Jacobian matrix vectors to be passed by reference to the function that assembles the Jacobian
        std::vector<uint32_t> col_index, row_position;
        std::vector<double> values;

        // Call function to assemble Jacobian
        Jacobian(Ly, q, node_coords, ELEMENTS, col_index, row_position, values);

        // Change row_position vector to include a row index for each NNZ value
        uint32_t r, r_next, row_nnzs;
        std::vector<uint32_t> rows;
        for (uint32_t i = 0; i < row_position.size() - 1; i++)
        {
            r = row_position[i];
            r_next = row_position[i + 1];
            row_nnzs = r_next - r;

            for (uint32_t n = 0; n < row_nnzs; n++)
            {
                rows.push_back(i);
            }
        }

        matlab::data::ArrayFactory f_values;
        matlab::data::ArrayFactory f_rows;
        matlab::data::ArrayFactory f_cols;

        outputs[0] = f_values.createArray<double>({1, values.size()}, values.data(), values.data() + values.size());
        outputs[1] = f_rows.createArray<uint32_t>({1, rows.size()}, rows.data(), rows.data() + rows.size());
        outputs[2] = f_cols.createArray<uint32_t>({1, col_index.size()}, col_index.data(), col_index.data() + col_index.size());
    }

    // This function assembles the vector the constitute the Jacobian matrix in sparse CSR format
    void Jacobian(double Ly, std::vector<double> q, std::vector<double> node_coords, std::vector<uint32_t> ELEMENTS, std::vector<uint32_t> &col_index, std::vector<uint32_t> &row_position, std::vector<double> &values)
    {

        uint32_t els = ELEMENTS.size() / 6;       // Number of elements
        uint32_t nodes = node_coords.size() / 2;  // Number of nodes
        uint32_t n_dofs = node_coords.size();     // Number of DOFs
        double x1, x2, y1, y2, x4, y4, x41;       // Element nodal coordinates and element side (nodes 1 to 3)
        uint32_t n1, n2, n4, dof12, dof22, dof42; // Element DOFs
        uint32_t counter = 0;                     // Initialise UDL counter. Max. value will be length of q-1

        // Create vector of DOFs
        std::vector<uint32_t> dofs(n_dofs);
        std::iota(std::begin(dofs), std::end(dofs), 0);

        // Initialise vector of tuples containing row and column indices of NNZs
        std::vector<std::tuple<uint32_t, uint32_t>> v;

        // Determine sparsity pattern
        for (uint32_t i = 0; i < els; i++)
        {
            // Identify element node numbers. SUBTRACT 1 DUE TO MATLAB INDEXING STARTING AT 1
            n1 = ELEMENTS[i * 6] - 1;
            n2 = ELEMENTS[i * 6 + 1] - 1;
            n4 = ELEMENTS[i * 6 + 3] - 1;

            // element node 1 - y coordinate
            y1 = node_coords[n1 * 2 + 1];
            // element node 2 - y coordinate
            y2 = node_coords[n2 * 2 + 1];
            // element node 4 - y coordinate
            y4 = node_coords[n4 * 2 + 1];

            dof12 = dofs[n1 * 2 + 1]; // element node 1 - dofs
            dof22 = dofs[n2 * 2 + 1];
            dof42 = dofs[n4 * 2 + 1];

            // If element is at the top edge of planar body
            if (y1 == Ly && y2 == Ly && y4 == Ly)
            {
                // Register all NNZ values. Some will be duplicate
                v.push_back(std::make_tuple(dof12, counter));
                v.push_back(std::make_tuple(dof22, counter));
                v.push_back(std::make_tuple(dof42, counter));
                counter++;
            }
        }

        // Sort vector and remove duplicates
        sort(v.begin(), v.end());
        v.erase(unique(v.begin(), v.end()), v.end());

        values.resize(v.size(), 0);         // Non-zero values of stiffness matrix
        row_position.resize(n_dofs + 1, 0); // Indices in nnz values of stiffness matrix at start of each row

        // Construct row_position and col_index vectors
        uint32_t count = 0; // Counter for number of nnzs
        for (uint32_t i = 0; i < n_dofs; i++)
        {
            // Vector containing column indices for row i
            std::vector<uint32_t> col_nnzs{};
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
                    for (uint32_t k = 0; k < col_nnzs.size(); k++)
                    {
                        // Append column indices to full column index vector
                        col_index.push_back(col_nnzs[k]);
                    }

                    // When all the NNZs in the row are registered, move on to the next row
                    break;
                }
            }
            row_position[i + 1] = count;
        }

        uint32_t r, r_next;
        counter = 0;
        // Calculate global stiffness matrix coefficients with known sparsity pattern
        for (uint32_t i = 0; i < els; i++)
        {
            // Identify element node numbers
            n1 = ELEMENTS[i * 6] - 1;
            n2 = ELEMENTS[i * 6 + 1] - 1;
            n4 = ELEMENTS[i * 6 + 3] - 1;

            // element node 1 - x,y coordinates
            x1 = node_coords[n1 * 2];
            y1 = node_coords[n1 * 2 + 1];
            // element node 2 - y coordinate
            y2 = node_coords[n2 * 2 + 1];
            // element node 4 - x,y coordinates
            x4 = node_coords[n4 * 2];
            y4 = node_coords[n4 * 2 + 1];

            // If element is at the top edge of planar body
            if (y1 == Ly && y2 == Ly && y4 == Ly)
            {
                x41 = x4 - x1;

                // element node 1 - dofs
                dof12 = dofs[n1 * 2 + 1];
                // element node 2 - dofs
                dof22 = dofs[n2 * 2 + 1];
                // element node 4 - dofs
                dof42 = dofs[n4 * 2 + 1];
                std::vector<uint32_t> e_dofs{};
                e_dofs.push_back(dof12);
                e_dofs.push_back(dof22);
                e_dofs.push_back(dof42);

                for (uint32_t j = 0; j < e_dofs.size(); j++)
                {
                    r = row_position[e_dofs[j]];
                    r_next = row_position[e_dofs[j] + 1];
                    for (uint32_t c = r; c < r_next; c++)
                    {
                        // If the column index matches with the UDL counter, add contribution to Jacobian
                        if (col_index[c] == counter)
                        {
                            if (j == 2) // The midpoint along the edge recives double the contribution than the vertices
                            {
                                values[c] += abs(x41) * 1e-3;
                                break;
                            }
                            else
                            {
                                values[c] += 0.5 * abs(x41) * 1e-3;
                                break;
                            }
                        }
                    }
                }
                counter++;
            }
        }
    };
};