# Systems of Linear Equations


This project implements iterative methods for solving linear systems of equations in the form \(Ax = b\), where \(A\) is a matrix and \(b\) is a vector. The methods implemented include **Jacobi Method**, **Gauss-Seidel Method**, and **LU Decomposition Method**. The primary goal is to compare the performance and convergence of these methods for large systems.

## Main Structure

- **Matrix**: Contains the structure and methods for working with matrices, including matrix creation, initialization, and LU decomposition.
- **Vector**: Defines a vector and methods to initialize it, including creating a vector of size \(N\) with sine values.
- **SolutionVector**: Inherits from the `Vector` structure and adds additional attributes for tracking iterations, residuals, and time.
- **Methods**: Implements the iterative methods such as Jacobi, Gauss-Seidel, and LU Decomposition.
- **CSV Saving**: Functionality for saving results to CSV files for further analysis.

## Features

- **Jacobi Method**: An iterative method that updates each element of the solution vector based on the previous values of all other elements. It is generally slower but simple.
- **Gauss-Seidel Method**: Similar to Jacobi but updates the solution immediately after computing each element. It often converges faster than the Jacobi method.
- **LU Decomposition**: Solves the system by first decomposing the matrix into lower (L) and upper (U) triangular matrices, then solving the system using forward and backward substitution.

## Usage

1. **Matrix and Vector Initialization**:
    - The `Matrix` and `Vector` structures define how matrices and vectors are initialized and stored. The matrix is a tridiagonal matrix with specified diagonal and off-diagonal elements.
    - `Matrix` can be initialized using different constructors:
      - Default constructor that initializes matrix with a given index.
      - Custom constructor for arbitrary sizes.
    - `Vector` can be initialized using similar methods, with its values based on a sine function.

2. **Iterative Methods**:
    - **Jacobi Method**: The `jacobiMethod` function iterates until convergence or a maximum number of iterations, solving the system \(Ax = b\).
    - **Gauss-Seidel Method**: The `gaussSeidlerMethod` function similarly solves the system iteratively but uses updated values immediately.
    - **LU Decomposition**: The `LUdecompositionMethod` function decomposes matrix \(A\) into \(L\) and \(U\) matrices and solves for \(x\).

3. **Residual Calculation**:
    - The `calculateResidualNorm` function computes the residual norm, which measures the difference between the left and right sides of the equation \(Ax - b\).
    - This residual is used to determine the convergence of the iterative methods.

4. **CSV Saving**:
    - After running the iterative methods, results are saved into `.csv` files for analysis. The files contain residual values and computation times for each method.

## Example Output

After running the methods, the program will output the time taken for each method and save residual data to CSV files.

