#ifndef PFUNCTION_H
#define PFUNCTION_H

#include "../utils/mpreal.h"

namespace math {
class CFunction;
class PFunction
{
public:
    PFunction();
    PFunction(int length);
    PFunction(const PFunction& orig);
    PFunction& operator=(const PFunction& orig);

    PFunction(const CFunction& orig);
    PFunction& operator=(const CFunction& orig);

    virtual ~PFunction();

    void resize(int length);
    int getSize() const { return m_length; }

    void set(int index, mpfr::mpreal argument, mpfr::mpreal value);
    void setArgument(int index, mpfr::mpreal argument);
    void setValue(int index, mpfr::mpreal value);
    mpfr::mpreal getArgument(int index) const;
    mpfr::mpreal getValue(int index) const;

    //convenience
    mpfr::mpreal integrate(mpfr::mpreal x_1, mpfr::mpreal x_2) const;
    mpfr::mpreal integrate_arg(mpfr::mpreal x_1, mpfr::mpreal x_2) const;

protected:
    mpfr::mpreal* m_arguments;
    mpfr::mpreal* m_values;
    int m_length;

    //befriend output operator
    friend std::ostream& operator<<(std::ostream& stream, const math::PFunction& p);

};
std::ostream& operator<<(std::ostream& stream, const math::PFunction& p);
}

#endif // PFUNCTION_H
