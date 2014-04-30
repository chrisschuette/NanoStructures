#ifndef VECTOR_H
#define VECTOR_H

namespace math {
class Vector
{
public:
    Vector();
    Vector(int length);
    Vector(const Vector& orig);
    Vector& operator =(const Vector& orig);
    ~Vector();

    inline double* getData() { return m_data; }

    inline double get(int index) const { return m_data[index]; }
    void set(int index, double value) { m_data[index] = value; }

    inline int getLength() { return m_length; }

    void resize(int length);
protected:
    double * m_data;
    int m_length;
};
}
#endif // VECTOR_H
