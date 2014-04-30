#ifndef UTILS_H
#define UTILS_H

#include <sstream>

namespace utils {
std::string timestamp();

template <class T>
std::string toString (T a) {
    std::stringstream ss;//create a stringstream
    ss << a;//add number to the stream
    return ss.str();//return a string with the contents of the stream
}
std::string toUpperCase(std::string strToConvert);
template <class T>
void swap (T& a, T& b) {
    T c = a;
    a = b;
    b = c;
}
}
#endif // UTILS_H
