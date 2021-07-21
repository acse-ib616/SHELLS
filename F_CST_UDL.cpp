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
        double t = inputs[0][0];      // Thickness of shell mm
        double weight = inputs[1][0]; // Unit weight of material kN/m3
        // double dx = inputs[2][0];     // Distance in x-direction between nodes mm

        // Change inputs into vectors, which are easier to manipulate in C++ and can link with sparse solvers
        // and matrix classes. Values are stored in ROW-MAJOR order in the vectors (A[i][j] = v[i*columns +j])
        matlab::data::TypedArray<double> coords_Array = std::move(inputs[2]);
        std::vector<double> node_coords(coords_Array.begin(), coords_Array.end()); // Node coordinates
        matlab::data::TypedArray<uint64_t> els_Array = std::move(inputs[3]);
        std::vector<uint64_t> ELEMENTS(els_Array.begin(), els_Array.end()); // Nodal connectivity
        matlab::data::TypedArray<double> F_Array = std::move(inputs[4]);
        std::vector<double> F(F_Array.begin(), F_Array.end()); // Nodal forces vector

        uint64_t n_dofs = node_coords.size(); // Number of DOFs

        // Call function to assemble K
        UDL(t, weight, node_coords, ELEMENTS, F);

        matlab::data::ArrayFactory f;

        // Assign MEX function output
        outputs[0] = f.createArray<double>({F.size(), 1}, F.data(), F.data() + F.size());
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

    // This function assembles the global nodal force vector for CST elements
    void UDL(double t, double weight, std::vector<double> node_coords, std::vector<uint64_t> ELEMENTS, std::vector<double> &forces)
    {
        uint64_t els = ELEMENTS.size() / 3;                            // Number of elements
        uint64_t nodes = node_coords.size() / 2;                       // Number of nodes
        uint64_t n_dofs = node_coords.size();                          // Number of DOFs
        double x1, x2, x3, y1, y2, y3, A, x21, x31, y21, y31;          // Element nodal coordinates, area & sides
        uint64_t n1, n2, n3, dof11, dof12, dof21, dof22, dof31, dof32; // Element DOFs

        // Create vector of DOFs
        std::vector<uint64_t> dofs(n_dofs);
        std::iota(std::begin(dofs), std::end(dofs), 0);

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

            // Triangle sides
            x21 = x2 - x1;
            x31 = x3 - x1;

            y21 = y2 - y1;
            y31 = y3 - y1;
            // Area of element
            A = abs(x21 * y31 - x31 * y21) / 2;

            // Weight contribution
            forces[dof12] = forces[dof12] - A * weight * t * 1e-9 / 3;
            forces[dof22] = forces[dof22] - A * weight * t * 1e-9 / 3;
            forces[dof32] = forces[dof32] - A * weight * t * 1e-9 / 3;
        }
    };
};