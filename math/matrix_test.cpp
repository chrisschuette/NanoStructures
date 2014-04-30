#include "matrix.h"

#include <gtest/gtest.h>

namespace math {
TEST(MatrixTest, DefaultConstructor) {
    const Matrix a;
    EXPECT_EQ(NULL, a.m_data);
    EXPECT_EQ(0, a.m_columns);
    EXPECT_EQ(0, a.m_rows);

}

TEST(MatrixTest, CopyConstructor) {
    const Matrix a((rand() % 50) + 1, (rand() % 50) + 1);

    int rows = a.getRows();
    int cols = a.getColumns();

    for(int i = 0; i < rows*cols; i++) {
        a.m_data[i] = rand();
    }

    const Matrix b(a);

    for(int i = 0; i < rows*cols; i++)
        EXPECT_EQ(b.m_data[i], a.m_data[i]);

    EXPECT_EQ(b.m_columns, a.m_columns);
    EXPECT_EQ(b.m_rows, a.m_rows);
}



}
