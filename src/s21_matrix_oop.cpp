#include "s21_matrix_oop.h"

/*
Contstrustors and destructor
*/
//===============================================================================================
// default constructor
S21Matrix::S21Matrix(): rows_(3), cols_(3), matrix_(nullptr) {
    AllocateMatrix();
}

// constructor with parameters
S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols), matrix_(nullptr) {
    AllocateMatrix();
}

// copy constructor
S21Matrix::S21Matrix(const S21Matrix& other) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = new double[rows_ * cols_];
    std::memcpy(matrix_, other.matrix_, sizeof(double) * rows_ * cols_);   
}

// move constructor
S21Matrix::S21Matrix(S21Matrix&& other) {
    double *tmp_ = other.matrix_;
    other.matrix_ = nullptr;
    int rows = other.rows_;
    int cols = other.cols_;
    // delete[] matrix_;
    matrix_ = tmp_;
    rows_ = rows;
    cols_ = cols;
}

//destructor
S21Matrix::~S21Matrix() {
    RemoveMatrix();
}

//===============================================================================================
/*
Accessors and mutators

*/
int S21Matrix::GetRows() const {
    return rows_;
}

int S21Matrix::GetColumns() const {
    return cols_;
}

void S21Matrix::SetRows(const int n) {
    if (n < 1) { std::out_of_range("In SetRows() incorrect rows"); }
    double* tmp_ = new double[n * cols_]();
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            if (i < n) tmp_[i * cols_ + j] = matrix_[i * cols_ + j];
        }
    }
    rows_ = n;
    delete[] matrix_;
    matrix_ = tmp_;
}

void S21Matrix::SetColumns(const int m) {
    if (m < 1) {std::out_of_range("In SetColumns() incorrect cols"); }
    double *tmp_ = new double[rows_ * m]();
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            if (j < m) tmp_[i * cols_ + j] = matrix_[i * cols_ + j];
        }
    }
    cols_ = m;
    delete[] matrix_;
    matrix_ = tmp_;
}

//================================================================================================
/*
Task functions

*/
bool S21Matrix::EqMatrix(const S21Matrix& other) const {
    if (!CheckSize(other) || !other.matrix_) { return false; }
    bool answer = true;
    for (int i = 0; i < rows_ && answer; i++) {
        for (int j = 0; j < cols_ && answer; j++) {
            if (fabs(matrix_[i * cols_ + j] - other.matrix_[i * cols_ + j]) > 1e-7) answer = false;
        }
    }
    return answer;
}

void S21Matrix::MulNumber(const double num) {
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j <cols_; j++) {
            matrix_[i * cols_ + j] *= num;
        }
    }
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
    if (!CheckSize(other) && other.matrix_) { std::runtime_error("Incorrect matrix in SumMatrix"); }
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            matrix_[i * cols_ + j] += other.matrix_[i * cols_ + j];
        }
    }
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
    if (!CheckSize(other) && other.matrix_) { std::runtime_error("Incorrect matrix in SubMatrix"); }
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            matrix_[i * cols_ + j] -= other.matrix_[i * cols_ + j];
        }
    }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
    if (cols_ != other.rows_) {
        throw std::runtime_error("Incorrect Matrix in MulMatrix");
    }
    S21Matrix tmp(rows_, other.cols_);
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < other.cols_; j++) {
            tmp.matrix_[i * other.cols_ + j] = 0;
            for (int k = 0; k < cols_; k++) {
                tmp.matrix_[i * cols_ + j] += matrix_[i * cols_ + k] * other.matrix_[k * other.cols_ + j];
            }
        }
    }
    std::swap(tmp.matrix_, matrix_);
    std::swap(tmp.rows_, rows_);
    std::swap(tmp.cols_, cols_);
}

S21Matrix S21Matrix::Transpose() const {
    S21Matrix tmp(cols_, rows_);
    for (int i = 0; i < cols_; i++) {
        for (int j = 0; j < rows_; j++) {
            tmp.matrix_[i * rows_ + j] = matrix_[j * cols_ + i]; 
        }
    }
    return tmp;
}

double S21Matrix::Determinant() const {
    if (rows_ != cols_) {throw std::runtime_error("Non quadratic matrix in Determinant"); }
    S21Matrix tmp(*this);
    double determinant = 1;
    for (int i = 0; i < rows_; i++) {
        int maxRow = i;
        for (int j = i + 1; j < rows_; j++) {
            if (std::abs(tmp.matrix_[j * cols_ + i]) > std::abs(tmp.matrix_[maxRow * cols_ + i])) {
                maxRow = j;
            }
        }
        if (i != maxRow) {
            SwapRows(tmp, i, maxRow);
            determinant *= -1.0; 
        }
        if (tmp.matrix_[i * cols_ + i] == 0.0) {
            return 0.0; 
        }
        determinant *= tmp.matrix_[i * cols_ + i];
        for (int j = i + 1; j < rows_; j++) {
            double factor = tmp.matrix_[j * cols_ + i] / tmp.matrix_[i * cols_ + i];
            for (int k = i; k < cols_; k++) {
                tmp.matrix_[j * cols_ + k] -= factor * tmp.matrix_[i * cols_ + k];
            }
        }
    }
    return determinant;
}

S21Matrix S21Matrix::CalcComplements() const {
    if (rows_ != cols_) { std::runtime_error("Non quadratic matrix in CalcComplements"); }
    S21Matrix tmp(rows_, cols_);
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            tmp.matrix_[i * cols_ + j] = Addition(i, j);     
        }
    }
    return tmp;
}

S21Matrix S21Matrix::InverseMatrix() const {
    if (rows_ != cols_) {throw std::runtime_error("Non quadratic matrix in InverseMatrix"); }
    double det = this->Determinant();
    std::cout << det << std::endl;
    if (fabs(det) < 1e-7) {throw std::runtime_error("Zero determinant in InverseMatrix"); }
    return S21Matrix(this->CalcComplements().Transpose() * (1 / det));
}

//=================================================================================================
/*
Task operators
*/
double& S21Matrix::operator()(int i, int j) const {
    if (0 <= i && i < rows_ && 0 <= j && j < cols_) {
        return matrix_[i * cols_ + j];
    } else {
        throw std::out_of_range("incorrect indexes in operator()");
    }
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
    this->SumMatrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
    this->SubMatrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
    this->MulMatrix(other);
    return *this;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
    S21Matrix tmp(other);
    std::swap(matrix_, tmp.matrix_);
    std::swap(cols_, tmp.cols_);
    std::swap(rows_, tmp.rows_);
    return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
    S21Matrix tmp(*this);
    tmp.SumMatrix(other);
    return tmp;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
    S21Matrix tmp(*this);
    tmp.SumMatrix(other);
    return tmp;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
    S21Matrix tmp(*this);
    tmp.MulMatrix(other);
    return tmp;
}

S21Matrix S21Matrix::operator*(const double num) {
    S21Matrix tmp(*this);
    tmp.MulNumber(num);
    return tmp;
}

bool S21Matrix::operator==(const S21Matrix& other) {
    return this->EqMatrix(other);
}

//=================================================================================================
/*
Helps functions

*/
void S21Matrix::AllocateMatrix() {
    try {
        matrix_ = new double[rows_ * cols_]();
    } catch (std::bad_alloc& ex) {
        std::cout << "Memory allocation error" << std::endl;
        throw;
    }
}

void S21Matrix::RemoveMatrix() {
    if (matrix_) {
        delete[] matrix_;
    }
}

void S21Matrix::SwapRows(S21Matrix& tmp, int k, int m) {
    for (int j = 0; j < cols_; j++) {
        std::swap(tmp.matrix_[m * cols_ + j], tmp.matrix_[k * cols_ + j]);
    }
}

bool S21Matrix::CheckSize(const S21Matrix& other) const {
    return cols_ == other.cols_ && rows_ == other.rows_;
}

void S21Matrix::Print() {
    std::cout << matrix_ << std::endl;
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            std::cout << matrix_[i * cols_ + j] << " ";
        }
        std::cout << std::endl;
    }
}

double S21Matrix::Addition(int i, int j) const {
    S21Matrix tmp(rows_ - 1, cols_ - 1);
    for (int k = 0; k < rows_ - 1; k++) {
        for (int m = 0; m < cols_ - 1; m++) {
            tmp.matrix_[k * (cols_ - 1) + m] = matrix_[((k < i) ? k : k + 1) * cols_ + ((m < j) ? m : m + 1)];
        }
    }
    return ((i + j) % 2 ? -1 : 1) * tmp.Determinant();
}
//=================================================================================================

S21Matrix operator*(const double num, const S21Matrix& matrix) {
    S21Matrix tmp(matrix);
    tmp.MulNumber(num);
    return tmp;
}