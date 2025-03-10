#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <string>
#include <fstream>


class Matrix {
private:
    int columns;
    int rows;
    std::vector<std::vector<double>> values;


public:
    Matrix(int rows, int columns) {
        this->columns = columns;
        this->rows = rows;
        this->values = std::vector<std::vector<double>>(this->rows, std::vector<double>(this->columns, 0));
    }

    Matrix(const Matrix* matrixToCopy) {
        this->columns = matrixToCopy->getColumns();
        this->rows = matrixToCopy->getRows();
        this->values = std::vector<std::vector<double>>(this->rows, std::vector<double>(this->columns, 0));
        for (int i = 0; i < this->rows; i++) {
            for (int j = 0; j < this->columns; j++) {
                this->values[i][j] = matrixToCopy->getValueInPosition(i,j);
            }
        }
    }

    ~Matrix() {

    }

    double getValueInPosition(int row, int column) const{
        if (row < this->rows && column < this->columns) {
            return this->values[row][column];
        }
        std::cout << "Indexes are out of range" << std::endl;
        exit(EXIT_FAILURE);
    }

    void setValueInPosition(int row, int column, double newValue) {
        if (row < this->rows && column < this->columns) {
            this->values[row][column] = newValue;
        }
        else {
            std::cout << "Indexes are out of range" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    void setColumns(int columns) {
        this->columns = columns;
    }

    void setRows(int rows) {
        this->rows = rows;
    }

    int getRows() const{
        return this->rows;
    }

    int getColumns() const{
        return this->columns;
    }

    void printMatrix() {
        for (int i = 0; i < this->rows; i++) {
            std::cout << "|";
            for (int j = 0; j < this->columns; j++) {
                std::cout << this->values[i][j];
                if (j != this->columns - 1) std::cout << " ";
            }
            std::cout << "|" << std::endl;
        }
        std::cout << std::endl;
    }


    void setVariablesForMatrixA(double a1, double a2, double a3) {
        int i = 0, j = 0;
        while (i < this->rows && j < this->columns) {
            this->setValueInPosition(i, j, a1);
            if (j + 1 < this->columns) this->setValueInPosition(i, j+1, a2);
            if (j + 2 < this->columns) this->setValueInPosition(i, j + 2, a3);
            if (j - 1 >= 0) this->setValueInPosition(i, j - 1, a2);
            if (j - 2 >= 0) this->setValueInPosition(i, j - 2, a3);
            i++;
            j++;
        }
    }

    void setVariablesForVectorB(double f) {
        int i = 0;
        while (i < this->rows) {
            this->setValueInPosition(i, 0, sin((i + 1) * (f + 1)));
            i++;
        }
    }

    static Matrix* matrixMultiplication(const Matrix* m1, const Matrix* m2) {
        if (m1->getColumns() != m2->getRows()) {
            std::cout << "Matrix dimensions are incorrect" << std::endl;
            exit(EXIT_FAILURE);
        }
        else {
            Matrix* result = new Matrix(m1->getRows(), m2->getColumns());
            for (int i = 0; i < m1->getRows(); i++) {
                for (int j = 0; j < m2->getColumns(); j++) {
                    double sum = 0;
                    for (int h = 0; h < m1->getColumns(); h++) {
                        sum += m1->getValueInPosition(i, h) * m2->getValueInPosition(h, j);
                    }
                    result->setValueInPosition(i, j, sum);
                }
            }
            return result;
        }
    }

    
    Matrix& operator=(const Matrix& other) {
        if (this != &other) {
            this->rows = other.rows;
            this->columns = other.columns;
            this->values = std::vector<std::vector<double>>(this->rows, std::vector<double>(this->columns, 0));
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0; j < this->columns; j++) {
                    this->values[i][j] = other.values[i][j];
                }
            }
        }
        return *this;
    }
    

    static Matrix* calculateResiduum(const Matrix* A, const Matrix* b, const Matrix* X) {
        Matrix* result = Matrix::matrixMultiplication(A, X);
        for (int i = 0; i < result->getRows(); i++) {
            result->setValueInPosition(i, 0, result->getValueInPosition(i, 0) - b->getValueInPosition(i, 0));
        }
        return result;
    }

    double calculateNorm() {
        double result = 0;
        for (int i = 0; i < this->getRows(); i++) {
            result += this->getValueInPosition(i, 0) * this->getValueInPosition(i, 0);
        }
        return sqrt(result);
    }


    static Matrix* solveGauss(const Matrix* A, const Matrix* b, std::vector<double>* norms, std::vector<double>* times, std::vector<int>* iters) {
        Matrix* X = new Matrix(b->getRows(), 1);
        for (int i = 0; i < X->getRows(); i++) {
            X->setValueInPosition(i, 0, 1);
        }
        int iterations = 0;
        Matrix* residuum = Matrix::calculateResiduum(A, b, X);
        double norm = residuum->calculateNorm();
        auto start = std::chrono::high_resolution_clock::now();
        while (norm > 10e-9 && iterations<100) {
            for (int i = 0; i < X->getRows(); i++) {
                double sum1 = 0, sum2 = 0;
                for (int j = 0; j < i; j++) {
                    sum1 += A->getValueInPosition(i, j) * X->getValueInPosition(j, 0);
                }
                for (int j = i + 1; j < X->getRows(); j++) {
                    sum2 += A->getValueInPosition(i, j) * X->getValueInPosition(j, 0);
                }
                double value = (b->getValueInPosition(i, 0) - sum1 - sum2) / A->getValueInPosition(i, i);
                X->setValueInPosition(i, 0, value);
            }
            iterations++;
            delete residuum;
            residuum = Matrix::calculateResiduum(A, b, X);
            norm = residuum->calculateNorm();
            //std::cout << "Iteration " << iterations << " " << norm << std::endl;
            (*norms).push_back(norm);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration<double>(end - start).count();
        std::cout << "Iterations " << iterations << std::endl;
        std::cout << "Time " << time << std::endl;
        (*times).push_back(time);
        (*iters).push_back(iterations);
        delete residuum;
        return X;
    }

    static Matrix* solveJacobi(const Matrix* A, const Matrix* b, std::vector<double>* norms, std::vector<double>* times, std::vector<int>* iters) {
        Matrix* X = new Matrix(b->getRows(), 1);
        for (int i = 0; i < X->getRows(); i++) {
            X->setValueInPosition(i, 0, 1);
        }
        int iterations = 0;
        Matrix* residuum = Matrix::calculateResiduum(A, b, X);
        double norm = residuum->calculateNorm();
        auto start = std::chrono::high_resolution_clock::now();
        while (norm > 10e-9 && iterations < 100) {
            Matrix* newX = new Matrix(X);
            for (int i = 0; i < X->getRows(); i++) {
                double sum1 = 0, sum2 = 0;
                for (int j = 0; j < i; j++) {
                    sum1 += A->getValueInPosition(i, j) * X->getValueInPosition(j, 0);
                }
                for (int j = i + 1; j < X->getRows(); j++) {
                    sum2 += A->getValueInPosition(i, j) * X->getValueInPosition(j, 0);
                }
                double value = (b->getValueInPosition(i, 0) - sum1 - sum2) / A->getValueInPosition(i, i);
                newX->setValueInPosition(i, 0, value);
            }
            iterations++;
            delete X;
            X = newX;
            delete residuum;
            residuum = Matrix::calculateResiduum(A, b, X);
            norm = residuum->calculateNorm();
            //std::cout << "Iteration " << iterations << " " << norm << std::endl;
            (*norms).push_back(norm);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration<double>(end - start).count();
        std::cout << "Iterations " << iterations << std::endl;
        std::cout << "Time " << time << std::endl;
        (*times).push_back(time);
        (*iters).push_back(iterations);
        delete residuum;
        return X;
    }

    static void createMatrixLU(const Matrix* A, Matrix* U, Matrix* L) {
        for (int j = 0; j < A->getColumns(); j++) {
            for (int i = 0; i < A->getRows(); i++) {
                if (i == j) L->setValueInPosition(i, j, 1);
                if (i <= j) {
                    double value = A->getValueInPosition(i, j);
                    for (int k = 0; k < i; k++) {
                        value -= L->getValueInPosition(i, k) * U->getValueInPosition(k, j);
                    }
                    U->setValueInPosition(i, j, value);
                }
                else {
                    double value = A->getValueInPosition(i, j);
                    for (int k = 0; k < i; k++) {
                        value -= L->getValueInPosition(i, k) * U->getValueInPosition(k, j);
                    }
                    L->setValueInPosition(i, j, value / U->getValueInPosition(j, j));
                }
            }
        }
    }

    static void calculateMatrixY(const Matrix* L, const Matrix* b, Matrix* y) {
        for (int i = 0; i < y->getRows(); i++) {
            double value = b->getValueInPosition(i, 0);
            for (int j = 0; j < i; j++) {
                value -= L->getValueInPosition(i, j) * y->getValueInPosition(j, 0);
            }
            y->setValueInPosition(i, 0, value);
        }
    }

    static void calculateMatrixX(const Matrix* U, const Matrix* y, Matrix* X) {
        for (int i = X->getRows()-1; i >= 0; i--) {
            double value = y->getValueInPosition(i, 0);
            for (int j = U->getColumns()-1; j > i; j--) {
                value -= X->getValueInPosition(j, 0) * U->getValueInPosition(i, j);
            }
            X->setValueInPosition(i, 0, value / U->getValueInPosition(i, i));
        }
    }

    static Matrix* MethodLU(const Matrix* A, const Matrix* b, std::vector<double>* times, std::vector<int>*iters) {
        Matrix* X = new Matrix(b->getRows(), 1);
        for (int i = 0; i < X->getRows(); i++) {
            X->setValueInPosition(i, 0, 1);
        }
        auto start = std::chrono::high_resolution_clock::now();
        Matrix* U = new Matrix(A->getRows(), A -> getColumns());
        Matrix* L = new Matrix(A->getRows(), A->getColumns());
        Matrix::createMatrixLU(A, U, L);
        Matrix* y = new Matrix(b->getRows(), 1);
        calculateMatrixY(L, b, y);
        calculateMatrixX(U, y, X);
        auto end = std::chrono::high_resolution_clock::now();
        delete U;
        delete L;
        delete y;
        auto time = std::chrono::duration<double>(end - start).count();
        std::cout << "Time " << time << std::endl;
        Matrix* residuum = Matrix::calculateResiduum(A, b, X);
        std::cout << "Norm " << residuum->calculateNorm() << std::endl;
        (*iters).push_back(1);
        (*times).push_back(time);
        delete residuum;
        return X;
    }

};


void saveNormsToFile(const std::string nameFile, const std::vector<double> v) {
    std::ofstream file(nameFile);

    if (!file.is_open()) {
        std::cout << "Error with open file" << std::endl;
        exit(EXIT_FAILURE);
    }

    file << "Iterations,Norms" << std::endl;

    for (int i = 0; i < v.size(); i++) {
        file << i + 1 << "," << v[i] << std::endl;
    }

    file.close();
}


void saveTimesToFile(const std::string nameFile, const std::vector<double> times, const std::vector<int> iterations, int size[]) {
    std::ofstream file(nameFile);

    if (!file.is_open()) {
        std::cout << "Error with open file" << std::endl;
        exit(EXIT_FAILURE);
    }

    file << "Sizes,Times,Iterations" << std::endl;

    for (int i = 0; i < times.size(); i++) {
        file << size[i] << "," << times[i] << "," << iterations[i] << std::endl;
    }

    file.close();
}


struct indexData {
    int size = 996;
    double a1 = 9;
    double a2 = -1;
    double a3 = -1;
    double f = 3;
} myIndexData;


void zadA(Matrix* A, Matrix* b) {
    A->setVariablesForMatrixA(myIndexData.a1, myIndexData.a2, myIndexData.a3);
    b->setVariablesForVectorB(myIndexData.f);
}


void zadB(Matrix* A, Matrix* b, std::vector<double> norms, std::vector<double> times, std::vector<int> iterations) {
    std::cout << "Gauss-Seidel: " << std::endl;
    Matrix* X = Matrix::solveGauss(A, b, &norms, &times, &iterations);
    saveNormsToFile("gaussNorms.csv", norms);
    norms.clear();
    times.clear();
    iterations.clear();
    delete X;


    std::cout << "Jacobi: " << std::endl;
    X = Matrix::solveJacobi(A, b, &norms, &times, &iterations);
    saveNormsToFile("jacobiNorms.csv", norms);
    norms.clear();
    times.clear();
    iterations.clear();
    delete X;


    std::cout << "Faktoryzacja LU: " << std::endl;
    X = Matrix::MethodLU(A, b, &times, &iterations);
    times.clear();
    iterations.clear();
    delete X;
}

void changeValuesInA(Matrix* A, int newA1) {
    myIndexData.a1 = newA1;
    A->setVariablesForMatrixA(myIndexData.a1, myIndexData.a2, myIndexData.a3);
}


void zadC(Matrix* A, Matrix* b, std::vector<double> norms, std::vector<double> times, std::vector<int> iterations) {
    std::cout << "Gauss-Seidel with a1 = 3: " << std::endl;
    Matrix* X = Matrix::solveGauss(A, b, &norms, &times, &iterations);
    saveNormsToFile("gaussNormsC.csv", norms);
    norms.clear();
    times.clear();
    iterations.clear();
    delete X;


    std::cout << "Jacobi with a1 = 3: " << std::endl;
    X = Matrix::solveJacobi(A, b, &norms, &times, &iterations);
    saveNormsToFile("jacobiNormsC.csv", norms);
    norms.clear();
    times.clear();
    iterations.clear();
    delete X;
}

void zadD(Matrix* A, Matrix* b, std::vector<double> times, std::vector<int> iterations) {
    std::cout << "Faktoryzacja LU with a1 = 3: " << std::endl;
    Matrix* X = Matrix::MethodLU(A, b, &times, &iterations);
    delete X;
    times.clear();
    iterations.clear();
}

void countGauss(std::vector<double> norms, std::vector<double> times, std::vector<int> iterations) {
    int size[] = { 100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000 };
    Matrix* A2;
    Matrix* b2;
    Matrix* X;

    for (int i = 0; i < sizeof(size)/sizeof(int); i++) {
        A2 = new Matrix(size[i], size[i]);
        b2 = new Matrix(size[i], 1);
        std::cout << "Gauss method for N = " << size[i] << std::endl;
        A2->setVariablesForMatrixA(myIndexData.a1, myIndexData.a2, myIndexData.a3);
        b2->setVariablesForVectorB(myIndexData.f);
        X = Matrix::solveGauss(A2, b2, &norms, &times, &iterations);
        delete A2;
        delete b2;
        delete X;
    }

    saveTimesToFile("gaussTimes.csv", times, iterations, size);
}

void countJacobi(std::vector<double> norms, std::vector<double> times, std::vector<int> iterations) {
    int size[] = { 100, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000 };
    Matrix* A2;
    Matrix* b2;
    Matrix* X;

    for (int i = 0; i < sizeof(size) / sizeof(int); i++) {
        A2 = new Matrix(size[i], size[i]);
        b2 = new Matrix(size[i], 1);
        std::cout << "Jacobi method for N = " << size[i] << std::endl;
        A2->setVariablesForMatrixA(myIndexData.a1, myIndexData.a2, myIndexData.a3);
        b2->setVariablesForVectorB(myIndexData.f);
        X = Matrix::solveJacobi(A2, b2, &norms, &times, &iterations);
        delete A2;
        delete b2;
        delete X;
    }

    saveTimesToFile("jacobiTimes.csv", times, iterations, size);
}

void countLU(std::vector<double> times, std::vector<int> iterations) {
    int size[] = { 100, 500, 1000, 1500, 2000, 2500, 3000 };
    Matrix* A2;
    Matrix* b2;
    Matrix* X;

    for (int i = 0; i < sizeof(size) / sizeof(int); i++) {
        A2 = new Matrix(size[i], size[i]);
        b2 = new Matrix(size[i], 1);
        std::cout << "LU method for N = " << size[i] << std::endl;
        A2->setVariablesForMatrixA(myIndexData.a1, myIndexData.a2, myIndexData.a3);
        b2->setVariablesForVectorB(myIndexData.f);
        X = Matrix::MethodLU(A2, b2, &times, &iterations);
        delete A2;
        delete b2;
        delete X;
    }

    saveTimesToFile("LUTimes.csv", times, iterations, size);
}

void zadE(std::vector<double> norms, std::vector<double> times, std::vector<int> iterations) {
    countGauss(norms, times, iterations);
    countJacobi(norms, times, iterations);
    countLU(times, iterations);
}

int main()
{
    
    Matrix* A = new Matrix(myIndexData.size, myIndexData.size);
    Matrix* b = new Matrix(myIndexData.size, 1);
    std::vector<double> norms;
    std::vector<double> times;
    std::vector<int> iterations;

    /*Zadanie A*/
    zadA(A, b);
    
    /*Zadanie B*/
    zadB(A, b, norms, times, iterations);
    
    /*Zadanie C*/
    changeValuesInA(A, 3);
    zadC(A, b, norms, times, iterations);
    
    /*Zadanie D*/
    zadD(A, b, times, iterations);
    
    /*Zadanie E*/
    changeValuesInA(A, 9);
    zadE(norms, times, iterations);
    
    delete b;
    delete A;
    return 0;
}
