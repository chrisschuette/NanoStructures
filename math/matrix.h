#ifndef MATRIX_H
#define MATRIX_H

#include <string.h>

namespace math {
class Matrix
{
public:
    Matrix();
    Matrix(int rows, int cols);
    Matrix(const Matrix& orig);
    Matrix& operator=(const Matrix& orig);

    inline int getColumns() const { return m_columns; }
    inline int getRows() const { return m_rows; }
    inline double* getData() { return m_data; }

    double get(int row, int column) const;
    void set(int row, int column, double value);

    void resize(int rows, int columns);

    //convenience
    void zero() { memset(m_data, 0, sizeof(double) * m_rows * m_columns); }

    void copy(Matrix& source, int offset_rows, int offset_cols, int size_rows, int size_columns, double prefactor);
    Matrix* crop(int offsetm_rows, int offsetm_cols, int sizem_rows, int sizem_cols);


    ~Matrix();
protected:
    double * m_data;
    int m_rows;
    int m_columns;

    //testing
 //   FRIEND_TEST(MatrixTest, DefaultConstructor);
 //   FRIEND_TEST(MatrixTest, CopyConstructor);

};
}

#endif // MATRIX_H
