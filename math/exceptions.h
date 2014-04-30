#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include "../err/exception.h"

namespace error {
namespace math {
class MathException : error::Exception
{
public:
    MathException(std::string msg) : error::Exception(msg) {}
    virtual ~MathException() throw() {}
};
class OutOfBounds : MathException
{
public:
    OutOfBounds() : MathException("Out of bounds.") {}
    OutOfBounds(std::string msg) : MathException(msg) {}
    virtual ~OutOfBounds() throw() {}
};
class Nan : MathException
{
public:
    Nan() : MathException("Out of bounds.") {}
    Nan(std::string msg) : MathException(msg) {}
    virtual ~Nan() throw() {}
};
}
}

#endif // EXCEPTIONS_H
