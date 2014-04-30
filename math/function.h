#ifndef FUNC_H
#define FUNC_H
#include <string>

namespace mpi {
    class OpenMPI;
}
namespace math {
class Function
{
public:
    Function();
    Function(int length);
    Function(const Function& orig);
    Function& operator=(const Function& orig);
    virtual ~Function();

    void resize(int length);
    int getSize() const { return m_length; }

    void set(int index, double argument, double value);
    void setArgument(int index, double argument);
    void setValue(int index, double value);
    double getArgument(int index) const;
    double getValue(int index) const;

    //interpolation
    double interpolate(double x);

    //IO
    void read(std::string filename);
    void write(std::string filename);
    friend class mpi::OpenMPI;

protected:
    double* m_arguments;
    double* m_values;
    int m_length;
};
}
#endif // FUNCTION_H
