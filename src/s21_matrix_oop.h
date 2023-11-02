#ifndef S21_MATRIX_OOP_H_
#define S21_MATRIX_OOP_H_




#include <new>
#include <iostream>
#include <cmath>
#include <cstring>

static int count = 0;

class S21Matrix {
    friend S21Matrix operator*(const double num, const S21Matrix& matrix);
    private:
        double* matrix_;
        int rows_, cols_;

        void AllocateMatrix();
        void RemoveMatrix();
        void SwapRows(S21Matrix& tmp, int k, int m);

        bool CheckSize(const S21Matrix& other) const;

        double Addition(int i, int j) const;

    public:
        S21Matrix();
        S21Matrix(int rows, int cols);
        S21Matrix(const S21Matrix& other);
        S21Matrix(S21Matrix&& other);
        ~S21Matrix();

        int GetRows() const;
        int GetColumns() const;
        void SetRows(const int n);
        void SetColumns(const int m);

        bool EqMatrix(const S21Matrix& other) const;

        void MulNumber(const double num);
        void SumMatrix(const S21Matrix& other);
        void SubMatrix(const S21Matrix& other);
        void MulMatrix(const S21Matrix& other);

        S21Matrix Transpose() const;
        S21Matrix InverseMatrix() const;
        S21Matrix CalcComplements() const;
        
        double Determinant() const; 

        double& operator()(int i, int ) const;
        S21Matrix& operator+=(const S21Matrix& other);
        S21Matrix& operator-=(const S21Matrix& other);
        S21Matrix& operator*=(const S21Matrix& other);
        S21Matrix& operator*=(const double num);
        S21Matrix& operator=(const S21Matrix& other);
        S21Matrix operator+(const S21Matrix& other);
        S21Matrix operator-(const S21Matrix& other);
        S21Matrix operator*(const S21Matrix& other);
        S21Matrix operator*(const double num);

        bool operator==(const S21Matrix& other);
        void Print();   
};

#endif // S21_MATRIX_OOP_H_