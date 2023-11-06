#include <iostream>

#include "s21_matrix_oop.h"

void outputMatrix(const S21Matrix& matrix) {}

int main() {
  S21Matrix my_matrix(6, 6);
  my_matrix(0, 0) = 1;
  my_matrix(0, 1) = 1;
  my_matrix(0, 2) = 1;
  my_matrix(0, 3) = 1;
  my_matrix(0, 4) = 1;
  my_matrix(0, 5) = 3;

  my_matrix(1, 0) = 2;
  my_matrix(1, 1) = 2;
  my_matrix(1, 2) = 2;
  my_matrix(1, 3) = 2;
  my_matrix(1, 4) = 2;
  my_matrix(1, 5) = 3;

  my_matrix(2, 0) = 4;
  my_matrix(2, 1) = 4;
  my_matrix(2, 2) = 4;
  my_matrix(2, 3) = 4;
  my_matrix(2, 4) = 4;
  my_matrix(2, 5) = 5;

  my_matrix(3, 0) = 5;
  my_matrix(3, 1) = 5;
  my_matrix(3, 2) = 5;
  my_matrix(3, 3) = 5;
  my_matrix(3, 4) = 5;
  my_matrix(3, 5) = 6;

  my_matrix(4, 0) = 6;
  my_matrix(4, 1) = 6;
  my_matrix(4, 2) = 6;
  my_matrix(4, 3) = 6;
  my_matrix(4, 4) = 6;
  my_matrix(4, 5) = 7;

  my_matrix(5, 0) = 7;
  my_matrix(5, 1) = 7;
  my_matrix(5, 2) = 7;
  my_matrix(5, 3) = 7;
  my_matrix(5, 4) = 7;
  my_matrix(5, 5) = 8;

  double det = my_matrix.Determinant();
  return 0;
}