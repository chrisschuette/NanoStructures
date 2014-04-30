#ifndef FUNCTION_H
#define FUNCTION_H

#include <complex>
#include <vector>

namespace mpi {
    class OpenMPI;
}
namespace math {
class PFunction;
class CFunction
{
public:
    typedef std::complex<double> (*tIntegrandFunction)(double, std::complex<double>, void*);
    inline static std::complex<double> one(double omega, std::complex<double> fomega, void* parameters) {
        return 1.0;
    }
    inline static std::complex<double> arg(double omega, std::complex<double> fomega, void* parameters) {
        return omega;
    }

public:
    CFunction();
    CFunction(int length);
    CFunction(const CFunction& orig);
    CFunction& operator=(const CFunction& orig);
    virtual ~CFunction();

    void resize(int length);
    int getSize() const { return m_length; }

    void set(int index, double argument, std::complex<double> value);
    void setArgument(int index, double argument);
    void setValue(int index, std::complex<double> value);
    double getArgument(int index) const;
    std::complex<double> getValue(int index) const;
    double getValueReal(int index) const;
    double getValueImag(int index) const;

    //IO
    void write(std::string filename) const;
    void read(std::string filename);
    void writeBinary(std::string filename) const;
    void readBinary(std::string filename);
    void dump() const;

    // operators
    void rescale(double factor);
    bool operator==(const CFunction& rhs) const;

    //convenience
    void createLogGrid(double maxFrequency, double lambda, int n);
    math::CFunction transformPH() const;

    // interpolation
    std::complex<double> interpolate(double x);
    math::CFunction interpolate(double x1, double x2) const;
    math::CFunction interpolate(double x1, double x2,  std::complex<double> C1, std::complex<double> C2 ,std::complex<double> D1, std::complex<double> D2) const;

    //integration
    std::complex<double> integrate(double x_1, double x_2) const;
    std::complex<double> integrate(double x_1, double x_2, tIntegrandFunction g, void * parameters) const;

protected:
    double* m_arguments;
    std::complex<double>* m_values;
    int m_length;

    //befriend PFunction
    friend class PFunction;

    //befriend OpenMPI
    friend class mpi::OpenMPI;

    //befriend the operators
    friend math::CFunction operator*(std::complex<double> factor, const math::CFunction& f);
    friend math::CFunction operator/(const math::CFunction& f, const math::CFunction& g);
    friend math::CFunction operator+(const math::CFunction& f, const math::CFunction& g);
    friend math::CFunction operator/(const math::CFunction& f, double denominator);
};

// operators

math::CFunction operator*(std::complex<double> factor, const math::CFunction& f);
math::CFunction operator/(const math::CFunction& f, const math::CFunction& g);
math::CFunction operator+(const math::CFunction& f, const math::CFunction& g);
math::CFunction operator/(const math::CFunction& f, double denominator);
double similarity(const math::CFunction& a, const math::CFunction& b);

}
#endif // FUNCTION_H
