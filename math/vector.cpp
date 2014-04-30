#include "vector.h"
#include <iostream>
#include <string.h>

namespace math {
Vector::Vector()
    : m_data(0)
    , m_length(0)
{
}

Vector::Vector(int length)
    : m_data(new double[length])
    , m_length(length)
{

}

Vector::Vector(const Vector& orig)
    : m_length(orig.m_length)
{
    if(orig.m_data == 0) {
        m_length = 0;
        m_data = 0;
    } else {
        m_data = new double[m_length];
        memcpy(m_data, orig.m_data, sizeof(double) * m_length);
    }
}

Vector& Vector::operator =(const Vector& orig) {
    if(m_data != 0) {
        delete [] m_data;
        m_data = 0;
        m_length = 0;
    }

    m_length = orig.m_length;

    if(orig.m_data == 0) {
        m_length = 0;
        m_data = 0;
    } else {
        m_data = new double[m_length];
        memcpy(m_data, orig.m_data, sizeof(double) * m_length);
    }
    return *this;
}

Vector::~Vector() {
    if(m_data)
        delete [] m_data;
}

void Vector::resize(int length) {
    if(m_data == 0) {
        m_data = new double[length];
        m_length = length;
    } else if(m_length != length) {
        delete [] m_data;
        m_length = length;
        if(length > 0)
            m_data = new double[length];
        else
            m_data = 0;
    }
}

}
