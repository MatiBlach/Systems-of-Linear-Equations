#include <math.h>
#include <String>
#include <chrono>
#include <Vector>
#include <fstream>
#include <iostream>

#define INDEX 193174

using namespace std;


struct Matrix {
    int N; //974
    double a1; // main diagonal
    double a2; //second diagonal
    double a3; //third diagonal
    double** matrix;

    Matrix(int index) {
        N = 900 + (index % 100);
        a1 = 5 + ((index % 1000) / 100);
        a2 = -1;
        a3 = -1;
        createMatrix();

    }
    Matrix(int size, int index) {
        N = size;
        a1 = 5 + ((index % 1000) / 100);
        a2 = -1;
        a3 = -1;
        createMatrix();
    }

    Matrix(int N, bool b) {
        this->N = N;
        this->a1 = b;
        this->a2 = b;
        this->a3 = b;
        createMatrix();
    }

    Matrix(int N, double a1, double a2, double a3) {
        this->N = N;
        this->a1 = a1;
        this->a2 = a2;
        this->a3 = a3;
        createMatrix();
    }

    void createMatrix() {
        matrix = new double* [N];
        for (int i = 0; i < N; i++) {
            matrix[i] = new double[N];
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) {
                    matrix[i][j] = a1;
                }
                else if (i > 0 && (i - 1) == j || j > 0 && (j - 1) == i) {
                    matrix[i][j] = a2;
                }
                else if (i > 1 && (i - 2) == j || j > 0 && (j - 2) == i) {
                    matrix[i][j] = a3;
                }
                else {
                    matrix[i][j] = 0;
                }
            }
        }
    }

    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            for (int i = 0; i < N; i++) {
                delete[] matrix[i];
            }
            delete[] matrix;
            N = other.N;
            a1 = other.a1;
            a2 = other.a2;
            a3 = other.a3;
            matrix = new double* [N];
            for (int i = 0; i < N; i++) {
                matrix[i] = new double[N];
            }
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    matrix[i][j] = other.matrix[i][j];
                }
            }
        }
        return *this;
    }

    ~Matrix() {
        for (int i = 0; i < N; i++) {
            delete[] matrix[i];
        }
        delete[] matrix;
    }

};

struct Vector {
    int N;
    double a;
    double* vector;

    Vector(int index) {
        N = 900 + (index % 100);
        a = 1 + ((index % 10000) / 1000);
        createVector();
    }

    Vector(int size, int index) {
        N = size;
        a = 1 + ((index % 10000) / 1000);
        createVector();
    }

    Vector(int N, double a) {
        this->N = N;
        this->a = a;
        createVector();
    }

    void createVector() {
        vector = new double[N];
        for (int i = 0; i < N; i++) {
            vector[i] = sin(i * a);
        }
    }

    ~Vector() {
        if (vector != nullptr) {
            delete[] vector;
            vector = nullptr;
        }
    }

    Vector& operator=(const Vector& other) {
        if (this != &other) {
            delete[] vector;
            N = other.N;
            a = other.a;
            vector = new double[N];
            for (int i = 0; i < N; i++) {
                vector[i] = other.vector[i];
            }
        }
        return *this;
    }
};

struct SolutionVector : Vector {
    int iterations = 0;
    double time = 0;
    std::vector<double> residual;

    SolutionVector(int index) : Vector(index) {}

    SolutionVector(int index, double a) : Vector(index, a) {}


    SolutionVector& operator=(const SolutionVector& other) {
        if (this != &other) {
            Vector::operator=(other);
            iterations = other.iterations;
            time = other.time;
        }
        return *this;
    }

    ~SolutionVector() {
        if (vector != nullptr) {
            delete[] vector;
            vector = nullptr;
        }

    }
};
//   res = Ax - b
//   |1  4  1|   |1|   |4|   |1*1 + 4*5 + 1*5|   |4|   |26|   |4|   |22|
//   |3  7  8| x |5| - |5| = |3*1 + 7*5 + 8*5| - |5| = |78| - |5| = |73|
//   |8  8  8|   |5|   |6|   |8*1 + 8*5 + 8*5|   |6|   |88|   |6|   |82|
double calculateResidualNorm(const Matrix& A, const Vector& b, const Vector& x) {
    int size = x.N;
    double res = 0;
    std::vector<double> residualVector;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            res += A.matrix[i][j] * x.vector[j];
        }
        res -= b.vector[i];
        residualVector.push_back(res);
    }
    double norm = 0;
    for (int i = 0; i < size; i++) {
        norm += residualVector[i] * residualVector[i];
    }
    return sqrt(norm);



}

//Iteration method calculating the solution based on previous one
SolutionVector jacobiMethod(const Matrix& A, const Vector& b, int max_iterations, double targetError) {
    auto start = std::chrono::high_resolution_clock::now();

    int x_size = b.N;
    SolutionVector x(x_size, 0);
    SolutionVector new_x(x_size, 0);
    std::vector<double> residual;

    for (int k = 0; k < max_iterations; k++) {

        for (int i = 0; i < x_size; i++) {
            double bi = b.vector[i];
            double sum_i = 0;
            double ai = A.matrix[i][i];

            for (int j = 0; j < x_size; j++) {
                if (i != j) //sum without diagonal element
                    sum_i += A.matrix[i][j] * x.vector[j];


            }
            if (ai != 0)
                new_x.vector[i] = (bi - sum_i) / ai;
            else
                printf("Error. Diagonal element cannot equal zero.");


        }

        x = new_x;
        residual.push_back(calculateResidualNorm(A, b, x));
        if (residual[k] < targetError || k + 1 == max_iterations) {
            x.iterations = k + 1;
            x.residual = residual;
            break;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    x.time = duration.count();
    return x;
}

//Iteration method calculating the solution based on the newest data
SolutionVector gaussSeidlerMethod(const Matrix& A, const Vector& b, int max_iterations, double targetError) {
    auto start = std::chrono::high_resolution_clock::now();

    int x_size = b.N;
    SolutionVector x(x_size, 0);
    SolutionVector new_x(x_size, 0);
    std::vector<double> residual;

    for (int k = 0; k < max_iterations; k++) {

        for (int i = 0; i < x_size; i++) {
            double bi = b.vector[i];
            double sum_i = 0;
            double ai = A.matrix[i][i];

            for (int j = 0; j < x_size; j++) {
                if (i > j) {
                    sum_i += A.matrix[i][j] * new_x.vector[j];
                }
                else if (i < j) {
                    sum_i += A.matrix[i][j] * x.vector[j];
                }

            }
            if (ai != 0)
                new_x.vector[i] = (bi - sum_i) / ai;
            else
                printf("Error. Diagonal element cannot equal zero.");


        }

        x = new_x;
        residual.push_back(calculateResidualNorm(A, b, x));
        if (residual[k] < targetError || k + 1 == max_iterations) {
            x.iterations = k + 1;
            x.residual = residual;
            break;
        }
    }


    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    x.time = duration.count();
    return x;
}


//L and U matrixes should be filled with zeros before
void createLU(Matrix& L, Matrix& U, const Matrix& A) {
    int size = A.N;
    U = A;
    for (int i = 0; i < size; i++) {
        L.matrix[i][i] = 1;
    }

    for (int i = 2; i <= size; i++) {
        for (int j = 1; j < i; j++) {
            L.matrix[i - 1][j - 1] = U.matrix[i - 1][j - 1] / U.matrix[j - 1][j - 1];
            for (int k = 1; k <= size; k++) {
                U.matrix[i - 1][k - 1] -= L.matrix[i - 1][j - 1] * U.matrix[j - 1][k - 1];
            }
        }
    }

}

Vector forwardSubstitution(const Matrix& L, const Vector& b) {
    int size = L.N;
    Vector y(size, 0);

    for (int i = 0; i < size; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) {
            sum += L.matrix[i][j] * y.vector[j];
        }
        y.vector[i] = (b.vector[i] - sum) / L.matrix[i][i];
    }

    return y;
}

SolutionVector backwardSubstitution(const Matrix& U, const Vector& y) {
    int size = U.N;
    SolutionVector x(size, 0);

    for (int i = size - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < size; ++j) {
            sum += U.matrix[i][j] * x.vector[j];
        }
        x.vector[i] = (y.vector[i] - sum) / U.matrix[i][i];
    }

    return x;
}


SolutionVector LUdecompositionMethod(const Matrix& A, const Vector& b) {
    auto start = std::chrono::high_resolution_clock::now();

    int size = A.N;
    Matrix L(size, 0);
    Matrix U(size, 0);
    createLU(L, U, A);
    Vector y = forwardSubstitution(L, b);
    SolutionVector x = backwardSubstitution(U, y);
    x.residual.push_back(calculateResidualNorm(A, b, x));

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    x.time = duration.count();
    return x;

}



void saveToCSV(const std::vector<std::vector<double>>& data, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& row : data) {
            for (size_t i = 0; i < row.size(); ++i) {
                file << row[i];
                if (i != row.size() - 1) {
                    file << ",";
                }
            }
            file << "\n";
        }
        file.close();
        std::cout << "Zapisano do " << filename << std::endl;
    }
    else {
        std::cerr << "Unable to open file " << filename << " for writing!" << std::endl;
    }
}

void saveTimesToCSV(const std::vector<int>& n_sizes, const std::vector<double>& jm_time,
    const std::vector<double>& gsm_time, const std::vector<double>& lu_time, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << "n_size,Jacobi_Method_Time,Gauss_Seidel_Method_Time,LU_Factorization_Time\n";

        for (size_t i = 0; i < n_sizes.size(); ++i) {
            file << n_sizes[i] << "," << jm_time[i] << "," << gsm_time[i] << "," << lu_time[i] << "\n";
        }

        file.close();
        std::cout << "Zapisano do " << filename << std::endl;
    }
    else {
        std::cerr << "Unable to open file " << filename << " for writing!" << std::endl;
    }
}

int main() {
    int n_size = 900 + (INDEX % 100);
    double targetError = 10e-9;
    int max_iterations = 1000;

    //Zadanie B
    Matrix A(INDEX);
    Vector b(INDEX);
    SolutionVector xA_JacobiMethod = jacobiMethod(A, b, max_iterations, targetError);
    SolutionVector xA_GaussSeidlerMethod = gaussSeidlerMethod(A, b, max_iterations, targetError);
    printf("Zadanie B\n");
    printf("Jacobi Method: %lfs \n", xA_JacobiMethod.time);
    printf("Gauss-Seidler Method: %lfs\n", xA_GaussSeidlerMethod.time);
    vector<vector<double>> dataB;
    dataB.push_back(xA_JacobiMethod.residual);
    dataB.push_back(xA_GaussSeidlerMethod.residual);
    saveToCSV(dataB, "dataB.csv");

    //Zadanie C
    Matrix C(n_size, 3, -1, -1);
    SolutionVector xC_JacobiMethod = jacobiMethod(C, b, max_iterations, targetError);
    SolutionVector xC_GaussSeidlerMethod = gaussSeidlerMethod(C, b, max_iterations, targetError);
    printf("Zadanie C\n");
    printf("Jacobi Method: %lfs \n", xC_JacobiMethod.time);
    printf("Gauss-Seidler Method: %lfs\n", xC_GaussSeidlerMethod.time);
    vector<vector<double>> dataC;
    dataC.push_back(xC_JacobiMethod.residual);
    dataC.push_back(xC_GaussSeidlerMethod.residual);
    saveToCSV(dataC, "dataC.csv");


    //Zadanie D
    Matrix D(n_size, 3, -1, -1);
    Vector e(INDEX);
    SolutionVector x = LUdecompositionMethod(D, e);
    printf("Zadanie D\n");
    printf("%.30lf %lfs\n", x.residual[0], x.time);


    //Zadanie E

    std::vector<int> n_sizes = { 100,500,1000,2000,3000 };
    std::vector<double> jm_time;
    std::vector<double> gsm_time;
    std::vector<double> lu_time;
    printf("Zadanie E\n");
    printf("N \t Jacobi \t Gauss-Seidler \t LU decomposition \n");
    for (auto n : n_sizes) {
        int max_iterations = INT_MAX;
        Matrix E(n, INDEX);
        Vector f(n, INDEX);
        SolutionVector x_JM = jacobiMethod(E, f, max_iterations, targetError);
        SolutionVector x_GSM = gaussSeidlerMethod(E, f, max_iterations, targetError);
        SolutionVector x_LU = LUdecompositionMethod(E, f);
        printf("%d \t %lf s \t %lf s \t %lf s \n", n, x_JM.time, x_GSM.time, x_LU.time);
        jm_time.push_back(x_JM.time);
        gsm_time.push_back(x_GSM.time);
        lu_time.push_back(x_LU.time);
    }
    saveTimesToCSV(n_sizes, jm_time, gsm_time, lu_time,"dataE.csv");
    



    return 0;
}