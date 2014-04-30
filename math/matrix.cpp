#include <assert.h>
#include "matrix.h"
#include "../err/exception.h"

namespace math {
Matrix::Matrix()
    : m_data(0)
    , m_rows(0)
    , m_columns(0)
{
}

Matrix::Matrix(int rows, int cols)
    : m_rows(rows)
    , m_columns(cols)
{
    assert(rows > 0);
    assert(cols > 0);
    m_data = 0;
    m_data = new double[rows * cols];
    if(m_data == 0)
        ERROR("matrix allocation failed.");
}

Matrix::Matrix(const Matrix& orig)
    : m_rows(orig.m_rows)
    , m_columns(orig.m_columns)
{
    int size = m_columns * m_rows;
    if(size > 0) {
        m_data = new double[size];
        if(m_data == 0)
            ERROR("matrix allocation failed.");
        assert(orig.m_data != 0);
        memcpy(m_data, orig.m_data, sizeof(double) * size);
    } else
        m_data = 0;
}

Matrix& Matrix::operator=(const Matrix& orig) {
    if(m_data != 0) {
        delete [] m_data;
        m_data = 0;
        m_columns = 0;
        m_rows = 0;
    }

    m_rows = orig.m_rows;
    m_columns = orig.m_columns;

    int size = m_columns * m_rows;
    if(size > 0) {
        m_data = new double[size];
        if(m_data == 0)
            ERROR("matrix allocation failed.");
        assert(orig.m_data != 0);
        memcpy(m_data, orig.m_data, sizeof(double) * size);
    } else
        m_data = 0;
    return *this;
}

Matrix::~Matrix() {
    if(m_data != 0)
        delete [] m_data;
}

void Matrix::set(int row, int col, double value) {
    m_data[row + col * m_rows] = value;
}

double Matrix::get(int row, int col) const {
    return m_data[row + col * m_rows];
}

void Matrix::resize(int rows, int columns) {
    int size = rows * columns;
    if(m_data == 0) {
        m_data = new double[size];
        m_rows = rows;
        m_columns = columns;
    } else if((m_rows != rows) || (m_columns != columns)) {
        delete [] m_data;
        m_rows = rows;
        m_columns = columns;
        if(size > 0)
            m_data = new double[size];
        else
            m_data = 0;
    }
}

void Matrix::copy(Matrix& source, int offset_rows, int offset_cols, int size_rows, int size_columns, double prefactor) {
    assert(source.m_rows <= size_rows);
    assert(source.m_columns <= size_columns);
    assert(m_rows >= offset_rows + size_rows);
    assert(m_columns >= offset_cols + size_columns);

    for (int i = 0; i < size_rows; i++)
        for (int j = 0; j < size_columns; j++)
            m_data[offset_rows + i + (offset_cols + j) * m_rows] = prefactor * source.m_data[i + j * source.m_rows];
}

Matrix* Matrix::crop(int offsetm_rows, int offsetm_cols, int sizem_rows, int sizem_cols) {
    Matrix* out = new Matrix(sizem_rows, sizem_cols);
    for (int i = 0; i < sizem_rows; i++)
        for (int j = 0; j < sizem_cols; j++)
            out->set(i, j, get(offsetm_rows + i, offsetm_cols + j));
    return out;
}


}
