#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>
#include <exception>

#define ERROR(message) throw error::Exception(message)

namespace error {
class Exception : public std::exception
{
public:
    Exception(std::string message);
     virtual const char* what() const throw() { return m_message.c_str(); }
    virtual ~Exception() throw() {}
protected:
    std::string m_message;
};
}
#endif // EXCEPTION_H
