#include "s21_matrix_oop.h"
#include <iostream>

void outputMatrix(const S21Matrix& matrix) {
    
}

int main() {
    S21Matrix matrix(3, 3);
    matrix(0, 0) = 9;
    matrix(0, 1) = 2;
    matrix(0, 2) = 3;
    matrix(1, 0) = 4;
    matrix(1, 1) = 0;
    matrix(1, 2) = 6;
    matrix(2, 0) = 7;
    matrix(2, 1) = 8;
    matrix(2, 2) = 0;
    matrix.Print();
    S21Matrix cc(matrix.InverseMatrix());
    cc.Print();
    return 0;
}