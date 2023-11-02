#include <gtest/gtest.h>
#include "../s21_matrix_oop.h"

class S21MatrixTest : public testing::Test {
protected:
    S21MatrixTest() {}

    virtual ~S21MatrixTest() {}

};

TEST_F(S21MatrixTest, ConstructorTest) {
    S21Matrix matrix; 
    ASSERT_EQ(matrix.GetRows(), 3);
    ASSERT_EQ(matrix.GetColumns(), 3);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}