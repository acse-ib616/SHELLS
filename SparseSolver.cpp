#include <iostream>
#include <math.h>
#include "SparseSolver.h"
#include <stdexcept>
#include <vector>
#include "utilities.h"
#include <memory>
#include <algorithm>

template <class T>
SparseSolver<T>::SparseSolver(CSRMatrix<T> &A, std::vector<T> &b) : A(A), b(b)
{
    // Check our dimensions match
    if (A.cols != b.size())
    {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
        return;
    }
}

// destructor
template <class T>
SparseSolver<T>::~SparseSolver()
{
}

template <class T>
T SparseSolver<T>::residualCalc(std::vector<T> &x, std::vector<T> &output_b)
{
    T residual = 0;
    // A x = b(estimate)
    A.matVecMult(x, output_b);

    // Find the norm between old value and new guess
    for (int i = 0; i < A.rows; i++)
    {
        residual += pow(output_b[i] - b[i], 2.0);
    }
    return sqrt(residual);
}

template <class T>
void SparseSolver<T>::stationaryIterative(std::vector<T> &x, double &tol, int &it_max, bool isGaussSeidel)
{
    double residual;
    std::vector<T> output_b(x.size(), 0);

    // vector for storing previous iteration if necessary
    std::vector<T> x_old;
    if (isGaussSeidel == false)
    {
        x_old = std::vector<T>(x.size(), 0);
    }

    // Check our dimensions match
    checkDimensions(A, b);
    checkDimensions(A, x);

    // Set values to zero before hand
    for (int i = 0; i < x.size(); i++)
    {
        x[i] = 0;
    }
    // loop up to a max number of iterations in case the solution doesn't converge
    int k;
    for (k = 0; k < it_max; k++)
    {
        // loop over rows
        for (int r = 0; r < A.rows; r++)
        {
            // Initialise sums of aij * xj
            T diagonal = 0;
            T sum = 0;

            // loop over non-zero values in row r
            for (int item_index = A.row_position[r]; item_index < A.row_position[r + 1]; item_index++)
            {
                int col_ind = A.col_index[item_index];
                if (r == col_ind)
                {
                    // this is a diagonal element
                    diagonal = A.values[item_index];
                }
                else if (!isGaussSeidel)
                {
                    sum += A.values[item_index] * x_old[col_ind];
                }
                else
                {
                    sum += A.values[item_index] * x[col_ind];
                }
            }

            x[r] = (1.0 / diagonal) * (b[r] - sum);
        }

        // Call residual calculation method
        residual = residualCalc(x, output_b);

        if (residual < tol)
        {
            break;
        }
        if (isGaussSeidel == false)
        {
            // Update the solution from previous iteration with
            // new estimate for Jacobi
            for (int i = 0; i < x.size(); i++)
            {
                x_old[i] = x[i];
            }
        }
    }
    std::cout << "k is :" << k << std::endl;
    std::cout << "residual is :" << residual << std::endl;
}

template <class T>
void SparseSolver<T>::conjugateGradient(std::vector<T> &x, double &tol, int &it_max)
{
    double residual;
    double alpha;
    double beta;
    std::vector<T> output_b(x.size(), 0);
    std::vector<T> residue_vec(x.size(), 0);
    std::vector<T> p(x.size(), 0);
    std::vector<T> Ap_product(x.size(), 0);
    std::vector<T> r_old(x.size(), 0);

    // Check our dimensions match
    checkDimensions(A, b);
    checkDimensions(A, x);

    // Set values to zero before hand
    for (int i = 0; i < x.size(); i++)
    {
        x[i] = 0;
    }

    // Find the norm between old value and new guess
    for (int i = 0; i < x.size(); i++)
    {
        r_old[i] = b[i];
        p[i] = r_old[i];
    }

    int k;
    for (k = 0; k < it_max; k++)
    {
        A.matVecMult(p, Ap_product);

        // Calculate alpha gradient
        alpha = vecDotProduct(r_old, r_old) / vecDotProduct(p, Ap_product);

        residual = 0.0;
        for (int i = 0; i < x.size(); i++)
        {
            x[i] += alpha * p[i];
            residue_vec[i] = r_old[i] - alpha * Ap_product[i];
            residual += pow(residue_vec[i], 2.0);
        }

        residual = sqrt(residual);

        if (residual < tol)
        {
            break;
        }

        // Calculate beta gradient
        beta = vecDotProduct(residue_vec, residue_vec) / vecDotProduct(r_old, r_old);
        for (int i = 0; i < x.size(); i++)
        {
            p[i] = residue_vec[i] + beta * p[i];

            // Update "old"(k) residue vector with new residue (k+1)
            // for next iteration
            r_old[i] = residue_vec[i];
        }
    }
    std::cout << "k is :" << k << std::endl;
    std::cout << "residual is :" << residual << std::endl;
}

// LU decomposition
template <class T>
std::shared_ptr<CSRMatrix<T>> SparseSolver<T>::lu_decomp()
/*
LU decomposition
Algorithm based on similar method as in 'Numerical recipes C++'.
The input matrix A is copied to LU, which is modified 'in place'.
Uses Crout's method by setting U_ii = 1.
Partial pivoting is implemented to ensure the stability of the method.
Implicit pivoting used to make it independent of scaling of equations.
An extensions would be to create a pivoting scheme that minimised the 
number of non-zero elements.
*/
{
    int n, max_ind, i, j, k, row_start, row_start2, row_len, col_indx;
    n = A.rows;
    T max, temp, piv_el;

    checkDimensions(A, b);

    std::vector<int> perm_indx(n);    // Store index of permutation
    std::vector<T> scaling(n);        // Store implicit scaling of each row
    std::vector<T> pivot_elements(n); // Store implicit scaling of each row

    // Implicit scaling, find max in each row and store scaling factor
    for (i = 0; i < n; i++)
    {
        max = 0.0;
        row_start = A.row_position[i];
        row_len = A.row_position[i + 1] - row_start;
        for (j = 0; j < row_len; j++)
        {
            temp = abs(A.values[row_start + j]);
            if (temp > max)
                max = temp;
        }
        if (max == 0)
            throw std::invalid_argument("Matrix is singular");
        scaling[i] = 1.0 / max;
    }

    bool matching = false;

    CSRMatrix<T> matrix_before = A;
    std::shared_ptr<CSRMatrix<T>> LU;
    while (!matching)
    {
        LU = matrix_before.matMatMultSymbolic(matrix_before);

        if (matrix_before.nnzs != LU->nnzs || matrix_before.rows != LU->rows)
        {
            // If they are not the same size
            matrix_before = *LU;
            continue;
        }

        bool allSame = true;

        for (int i = 0; i < LU->rows; i++)
        {
            if (matrix_before.row_position[i] != LU->row_position[i])
            {
                allSame = false;
                break;
            }
        }

        if (!allSame)
        {
            matrix_before = *LU;
            continue;
        }

        for (int i = 0; i < LU->nnzs; i++)
        {
            if (matrix_before.col_index[i] != LU->col_index[i])
            {
                allSame = false;
                break;
            }
        }

        if (!allSame)
        {
            matrix_before = *LU;
            continue;
        }

        matching = true;
    }

    // the lower and upper triangular matrix values are stored in the LU matrix
    // loop over rows of LU
    for (int row = 0; row < LU->rows; row++)
    {
        // loop over non-zero columns in that row
        for (int cii = LU->row_position[row]; cii < LU->row_position[row + 1]; cii++)
        {
            // col index of a non-zero
            int col = LU->col_index[cii];

            // first element is simply equal to A_00
            if (col == 0 && row == 0)
            {
                LU->values[0] = A.values[0];
                continue;
            }

            T a_ij = 0.0;
            // search in our original matrix for a_ij
            for (int a_row_pos = A.row_position[row]; a_row_pos < A.row_position[row + 1]; a_row_pos++)
            {
                if (A.col_index[a_row_pos] == col)
                {
                    a_ij = A.values[a_row_pos];
                    break;
                }
            }

            // number of products to sum over (k-limit)
            int cutoff = std::min(row, col) - 1;

            // sum over alpha_ik * beta_kj
            // look for values in same row first, if they exist check for col equivalents
            T valsum = 0.0;
            for (int col_vec_ind = LU->row_position[row]; col_vec_ind < LU->row_position[row + 1]; col_vec_ind++)
            {
                // check on this row, preceding the current value
                int k = LU->col_index[col_vec_ind];

                // if k is above cutoff, break the loop
                if (k > cutoff)
                {
                    break;
                }

                // if alpha_ik is non-zero
                if (LU->values[col_vec_ind])
                {
                    // check whether corresponding beta_kj also exists -> add to valsum
                    for (int tempcols = LU->row_position[k]; tempcols < LU->row_position[k + 1]; tempcols++)
                    {
                        if (LU->col_index[tempcols] == col)
                        {
                            valsum += (LU->values[col_vec_ind] * LU->values[tempcols]);
                        }
                    }
                }
            }

            if (col >= row)
            {
                // this means that we are in the upper triangle or the diagonal (U)
                LU->values[cii] = a_ij - valsum;
            }
            else if (col < row)
            {
                // we are in the lower triangle (L)
                // we need the beta value from the LU_jj above
                T b_jj = 0.0;

                for (int LU_vals_index = LU->row_position[col]; LU_vals_index < LU->row_position[col + 1]; LU_vals_index++)
                {
                    int temp_col_ind = LU->col_index[LU_vals_index];
                    if (temp_col_ind == col)
                    {
                        b_jj = LU->values[LU_vals_index];
                        break;
                    }
                }

                if (b_jj == 0)
                {
                    std::cerr << "Error: b_jj is zero" << std::endl;
                }

                LU->values[cii] = (1.0 / b_jj) * (a_ij - valsum);
            }
        }
    }

    return LU;
}

// Linear solver that uses LU decomposition
template <class T>
void SparseSolver<T>::lu_solve(CSRMatrix<T> &LU, std::vector<int> &perm_indx, std::vector<T> &x)
// Solve the equations L*y = b and U*x = y to find x.
{
    int n, ip, i, j, row_start, row_len, col_start, col_indx;
    n = LU.rows;
    T sum;

    checkDimensions(A, x);

    // The unknown x will be used as temporary storage for y.
    // The equations for forward and backward substitution have
    // been simplified by combining (b and sum) and (y and sum).
    for (i = 0; i < n; i++)
    {
        x[i] = b[i];
    }
    // Perform forward substitution to solve L*y = b.
    // Need to keep track of permutation of RHS as well
    for (i = 0; i < n; i++)
    {
        ip = perm_indx[i];
        sum = x[ip];
        x[ip] = x[i];
        // row_start = LU.row_position[i];
        // row_len = LU.row_position[i + 1] - row_start;
        for (j = LU.row_position[i]; j < LU.row_position[i + 1]; j++)
        {
            col_indx = LU.col_index[j];
            // check if valid and exits loop to avoid uneccesary checks
            if (col_indx >= i)
                break;
            sum -= LU.values[j] * x[col_indx];
        }
        x[i] = sum;
    }

    // Perform backward substitution to solve U*x = y
    // Here x = y before being updated.
    for (i = n - 1; i >= 0; i--)
    {
        sum = x[i];
        row_start = LU.row_position[i];
        row_len = LU.row_position[i + 1] - row_start;
        for (j = 0; j < row_len; j++)
        {
            col_indx = LU.col_index[row_start + j];
            if (col_indx >= i + 1)
                sum -= LU.values[row_start + j] * x[col_indx];
        }
        // Find diagonal element
        for (j = 0; j < row_len; j++)
        {
            col_indx = LU.col_index[row_start + j];
            if (col_indx == i)
            {
                x[i] = sum / LU.values[row_start + j];
                break;
            }
        }
    }
}

template <class T>
std::shared_ptr<CSRMatrix<T>> SparseSolver<T>::cholesky_decomp()
{
    // for now, we assume output has been preallocated
    std::vector<T> R_values{};
    std::vector<int> R_cols{};
    std::vector<int> R_row_position(A.rows + 1, 0);

    // loop over rows of A
    for (int i = 0; i < A.rows; i++)
    {
        // rows indices of left matrix
        int r_start = A.row_position[i];
        int r_end = A.row_position[i + 1];

        // loop over column indices for this row in A - equivalent to rows in B
        std::vector<int> cols_nnzs{};
        // cii - index of col_index of R array
        int ci = 0;
        std::vector<int> infills_left{};
        for (int cii = r_start; (cii < r_end && ci < i); cii++)
        {
            ci = A.col_index[cii];

            for (int k = 0; k < ci; k++)
            {
                if (cii != r_start && k == A.col_index[cii - 1])
                {
                    infills_left.push_back(A.col_index[cii - 1]);
                }
                for (int r = 0; r < k && k > 0; r++)
                {
                    int r_start_above = A.row_position[r];
                    int r_end_above = A.row_position[r + 1];
                    for (int a = r_start_above; a < r_end_above; a++)
                    {
                        if (infills_left.size() > 0 && A.col_index[a] == k)
                        {
                            cols_nnzs.push_back(k);
                        }
                    }
                }
            }
            if (ci < i)
            {
                cols_nnzs.push_back(ci);
            }
        }
        cols_nnzs.push_back(i);
        // sort non_zeros & remove duplicates
        sort(cols_nnzs.begin(), cols_nnzs.end());
        cols_nnzs.erase(unique(cols_nnzs.begin(), cols_nnzs.end()), cols_nnzs.end());

        // Concatenate
        for (int j = 0; j < cols_nnzs.size(); j++)
        {
            R_cols.push_back(cols_nnzs[j]);
        }
        R_row_position[i + 1] = R_row_position[i] + cols_nnzs.size();

        if (i == 0)
        {
            R_values.push_back(sqrt(A.values[0]));
        }

        // Entries in R on row i
        int R_r_start = R_row_position[i];
        int R_r_end = R_row_position[i + 1];
        if (i != 0)
        {
            for (int cii = R_r_start; cii < R_r_end; cii++)
            {
                ci = R_cols[cii];

                if (ci == 0)
                {

                    T A_i0;
                    for (int k = r_start; k < r_end; k++)
                    {
                        if (A.col_index[k] == 0)
                        {
                            A_i0 = A.values[k];
                            R_values.push_back(A_i0 / R_values[0]);
                        }
                    }
                }

                T sum_ij = 0;
                if (ci != 0 && ci != i)
                {

                    int r_start_above = R_row_position[ci];
                    int r_end_above = R_row_position[ci + 1];

                    for (int n = R_r_start; n < R_r_end; n++)
                    {
                        // Loop through columns[:j] in row ci above current entry row[i, ci]
                        for (int a = r_start_above; a < r_end_above; a++)
                        {
                            // Match cols in current and above rows
                            if (R_cols[n] == R_cols[a] && R_cols[n] != ci)
                            {
                                sum_ij += R_values[n] * R_values[a];
                            }
                        }
                    }

                    int diag_r_start = R_row_position[ci];
                    int diag_r_end = R_row_position[ci + 1];

                    // Loop through row that contains diagonal L entry on current column
                    for (int diag = diag_r_start; diag < diag_r_end; diag++)
                    {
                        if (R_cols[diag] == ci)
                        {
                            T A_ij = 0;
                            // Determine whether ij entry in L also appears in A, else we use a 0
                            for (int k = r_start; k < r_end; k++)
                            {
                                if (ci == A.col_index[k])
                                {
                                    A_ij = A.values[k];
                                }
                            }

                            R_values.push_back((A_ij - sum_ij) / R_values[diag]);
                        }
                    }
                }

                T sum_jj = 0;
                if (ci != 0 && ci == i)
                {
                    // Loop over Ljk entries in row j
                    for (int k = R_r_start; (k < R_r_end && R_cols[k] < i); k++)
                    {
                        sum_jj += pow(R_values[k], 2);
                    }

                    // Access jj value in A
                    T A_jj;
                    int diag;
                    for (diag = r_start; diag < r_end; diag++)
                    {
                        if (A.col_index[diag] == i)
                        {
                            A_jj = A.values[diag];
                        }
                    }

                    // Ljj = sqrt(Ajj - sum(Ljk[:j]^2)
                    R_values.push_back(sqrt(A_jj - sum_jj));
                }
            }
        }
    }

    std::shared_ptr<CSRMatrix<T>> sparse_mat_ptr(new CSRMatrix<T>(A.rows, A.cols, R_values.size(), true));

    for (int i = 0; i < R_values.size(); i++)
    {
        sparse_mat_ptr->col_index[i] = R_cols[i];
        sparse_mat_ptr->values[i] = R_values[i];
    }
    for (int i = 0; i <= A.rows; i++)
    {
        sparse_mat_ptr->row_position[i] = R_row_position[i];
    }

    return sparse_mat_ptr;
}

// Linear solver that uses LU decomposition
template <class T>
void SparseSolver<T>::cholesky_solve(CSRMatrix<T> &R, std::vector<T> &x)
// Solve the equations L*y = b and U*x = y to find x.
{
    int n, ip, i, j, row_start, row_len, col_start, col_indx;
    n = R.rows;
    T sum;

    checkDimensions(A, x);

    std::shared_ptr<CSRMatrix<T>> R_T = R.transpose();

    // The unknown x will be used as temporary storage for y.
    // The equations for forward and backward substitution have
    // been simplified by combining (b and sum) and (y and sum).
    for (i = 0; i < n; i++)
    {
        x[i] = b[i];
    }
    // Perform forward substitution to solve L*y = b.
    // Need to keep track of permutation of RHS as well
    for (i = 0; i < n; i++)
    {
        sum = x[i];

        for (j = R.row_position[i]; j < R.row_position[i + 1]; j++)
        {
            col_indx = R.col_index[j];
            // check if valid and exits loop to avoid uneccesary checks
            if (col_indx >= i)
                break;
            sum -= R.values[j] * x[col_indx];
        }
        x[i] = sum / R.values[j];
    }

    // Perform backward substitution to solve U*x = y
    // Here x = y before being updated.
    for (i = n - 1; i >= 0; i--)
    {
        sum = x[i];
        row_start = R_T->row_position[i];
        row_len = R_T->row_position[i + 1] - row_start;
        for (j = 0; j < row_len; j++)
        {
            col_indx = R_T->col_index[row_start + j];
            if (col_indx >= i + 1)
                sum -= R_T->values[row_start + j] * x[col_indx];
        }
        // Find diagonal element
        for (j = 0; j < row_len; j++)
        {
            col_indx = R_T->col_index[row_start + j];
            if (col_indx == i)
            {
                x[i] = sum / R_T->values[row_start + j];
                break;
            }
        }
    }
}