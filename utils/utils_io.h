#ifndef UTILS_IO_H
#define UTILS_IO_H

#include <istream>
#include <ostream>
#include <sys/stat.h>
#include <sys/types.h>

namespace utils {
namespace io {
template<typename T>
std::ostream& binary_write(std::ostream* outStream, const T& value){
    return outStream->write(reinterpret_cast<const char*>(&value), sizeof(T));
}
template<typename T>
std::istream & binary_read(std::istream* stream, T& value){
    return stream->read(reinterpret_cast<char*>(&value), sizeof(T));
}

void createDirectoryStructure(std::string filePath);
int mkpath(const char *path, mode_t mode);
static int do_mkdir(const char *path, mode_t mode);
void appendToFile(std::string filename, int num, ...);


}
}
#endif // UTILS_IO_H
