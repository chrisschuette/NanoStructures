#include "cfunction.h"
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "../utils/utils_io.h"
#include "exceptions.h"

namespace math {
CFunction::CFunction()
    : m_arguments(0)
    , m_values(0)
    , m_length(0)
{
}

CFunction::CFunction(int length)
    : m_arguments(new double[length])
    , m_values(new std::complex<double>[length])
    , m_length(length)
{
}

CFunction::CFunction(const CFunction& orig)
    : m_length(orig.m_length)
{
    if(orig.m_arguments != 0) {
        m_arguments = new double[m_length];
        memcpy(m_arguments, orig.m_arguments, sizeof(double) * m_length);
    }
    if(orig.m_values != 0) {
        m_values = new std::complex<double>[m_length];
        memcpy(m_values, orig.m_values, sizeof(std::complex<double>) * m_length);
    }
}

CFunction& CFunction::operator=(const CFunction& orig) {
    m_length = orig.m_length;
    if(m_arguments != 0)
        delete [] m_arguments;
    if(m_values != 0)
        delete [] m_values;

    if(orig.m_arguments != 0) {
        m_arguments = new double[m_length];
        memcpy(m_arguments, orig.m_arguments, sizeof(double) * m_length);
    }
    if(orig.m_values != 0) {
        m_values = new std::complex<double>[m_length];
        memcpy(m_values, orig.m_values, sizeof(std::complex<double>) * m_length);
    }
    return *this;
}

CFunction::~CFunction() {
    if(m_arguments != 0)
        delete [] m_arguments;
    if(m_values != 0)
        delete [] m_values;
}
void CFunction::set(int index, double argument, std::complex<double> value) {
    assert(index >= 0);
    assert(index < m_length);
    m_arguments[index] = argument;
    m_values[index] = value;
}

void CFunction::setArgument(int index, double argument) {
    assert(index >= 0);
    assert(index < m_length);
    m_arguments[index] = argument;
}

void CFunction::setValue(int index, std::complex<double> value) {
    assert(index >= 0);
    assert(index < m_length);
    m_values[index] = value;
}

double CFunction::getArgument(int index) const {
    assert(index >= 0);
    assert(index < m_length);
    return m_arguments[index];
}

std::complex<double> CFunction::getValue(int index) const {
    assert(index >= 0);
    assert(index < m_length);
    return m_values[index];
}

double CFunction::getValueReal(int index) const {
    assert(index >= 0);
    assert(index < m_length);
    return m_values[index].real();
}

double CFunction::getValueImag(int index) const {
    assert(index >= 0);
    assert(index < m_length);
    return m_values[index].imag();
}

void CFunction::read(std::string filename) {
    std::ifstream inFile(filename.c_str());
    if(!inFile.is_open()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        throw std::exception();
    }
    std::vector<double> data[3];
    double f, r, i;
    while (inFile >> f >> i >> r)
    {
        data[0].push_back(f);
        data[1].push_back(r);
        data[2].push_back(i);
    }
    //    std::cout << "Read " << data[0].size() << " entries." << std::endl;
    resize(data[0].size());
    if(data[0].size() == 0)
        std::cout << "WARNING: Function contains no data." << std::endl;
    for(int i = 0; i < data[0].size(); i++) {
        m_arguments[i] = data[0][i];
        m_values[i] = std::complex<double>(data[1][i], data[2][i]);
    }
}

void CFunction::write(std::string filename) const {
    utils::io::createDirectoryStructure(filename);
    std::ofstream outFile(filename.c_str());
    if(!outFile.is_open()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        throw std::exception();
    }
    if(m_length == 0)
        std::cout << "WARNING: Function contains no data." << std::endl;
    for(int i = 0; i < m_length; i++)
        outFile << m_arguments[i] << " " << m_values[i].imag() << " " << m_values[i].real() << std::endl;
    outFile.close();
}

void CFunction::writeBinary(std::string filename) const {
    utils::io::createDirectoryStructure(filename);
    using namespace utils::io;
    std::ofstream outFile(filename.c_str(), std::ios::binary);
    if(!outFile.is_open()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        throw std::exception();
    }
    if(m_length == 0)
        std::cout << "WARNING: Function contains no data." << std::endl;

    binary_write(&outFile, m_length);
    for(int i = 0; i < m_length; i++)
    {
        binary_write(&outFile, m_arguments[i]);
        binary_write(&outFile, m_values[i]);
    }
    outFile.close();
}

void CFunction::readBinary(std::string filename) {
    using namespace utils::io;
    std::ifstream inFile(filename.c_str(), std::ios::binary);
    if(!inFile.is_open()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        throw std::exception();
    }
    int length;
    binary_read(&inFile, length);
    resize(length);
    if(length == 0)
        std::cout << "WARNING: Function contains no data." << std::endl;

    for(int i = 0; i < length; i++) {
        binary_read(&inFile, m_arguments[i]);
        binary_read(&inFile, m_values[i]);
    }
    inFile.close();
}

void CFunction::dump() const {
    for(int i = 0; i < getSize(); i++)
        std::cout << getArgument(i) << " " << getValue(i) << std::endl;
}

void CFunction::rescale(double factor) {
    for(int n = 0; n < m_length; n++) {
        m_values[n] *= factor;
        m_arguments[n] *= factor;
    }
}


bool CFunction::operator==(const CFunction& rhs) const {
    if(m_length != rhs.m_length)
        return false;
    for(int i = 0; i < m_length; i++)
        if((m_arguments[i] != rhs.m_arguments[i]) || (m_values[i] != rhs.m_values[i]))
            return false;
    return true;
}



void CFunction::resize(int length) {
    if(m_length == length)
        return;
    if(m_arguments != 0)
        delete [] m_arguments;
    if(m_values != 0)
        delete [] m_values;

    if(length == 0)
    {
        m_length = length;
        m_arguments = 0;
        m_values = 0;
        return;
    } else {
        m_length = length;
        m_arguments = new double[m_length];
        m_values = new std::complex<double>[m_length];
    }
}

math::CFunction operator*(std::complex<double> factor, const math::CFunction& f) {
    math::CFunction result(f.getSize());
    for(int i = 0; i < f.getSize(); i++)
        result.set(i, f.getArgument(i), f.getValue(i) * factor);
    return result;
}

math::CFunction operator/(const math::CFunction& f, const math::CFunction& g) {
    assert(g.getSize() == f.getSize());
    math::CFunction result(f.getSize());
    for(int i = 0; i < f.getSize(); i++) {
        assert(f.getArgument(i) == g.getArgument(i));
        result.set(i, f.getArgument(i), f.getValue(i) / g.getValue(i));
    }
    return result;
}

math::CFunction operator+(const math::CFunction& f, const math::CFunction& g) {
    assert(f.getSize() == g.getSize());
    math::CFunction result(f.getSize());
    for(int i = 0; i < f.getSize(); i++) {
        assert(f.getArgument(i) == g.getArgument(i));
        result.set(i, f.getArgument(i), f.getValue(i) + g.getValue(i));
    }
    return result;
}

math::CFunction operator/(const math::CFunction& f, double denominator) {
    math::CFunction result(f.getSize());
    for(int i = 0; i < f.getSize(); i++) {
        result.set(i, f.getArgument(i), f.getValue(i) / denominator);
    }
    return result;
}

std::complex<double> CFunction::interpolate(double x) {
    if(x < m_arguments[0])
    {
        std::string message = "Left boundary outside value interval.";
        std::cerr << x << std::endl;
        throw error::math::OutOfBounds();
    }

    if(x > m_arguments[m_length-1])
    {
        std::string message = "Right boundary outside value interval.";
        std::cerr << x << std::endl;

        throw error::math::OutOfBounds();
    }

    int i = 0;
    while( (i < m_length) && (m_arguments[i] < x) )
    {
        i++;
    }
    if(i == m_length)
    {
        std::string message = "Left boundary outside value interval.";
        std::cerr << x << std::endl;

        throw error::math::OutOfBounds();
    }

    //interpolate between i and i-1

    double arg_left = m_arguments[i-1];
    double arg_right = m_arguments[i];

    std::complex<double> left = m_values[i-1];
    std::complex<double> right = m_values[i];

    //		std::cerr << "Interpolation(right) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

    std::complex<double> value = (right - left) / ( arg_right - arg_left ) * (x - arg_left) + left;
    //		std::cerr << "interpolated right value: (" << x_2 << "," << value_x_2 << ")" << std::endl;
    return value;
}


math::CFunction CFunction::interpolate(double x1, double x2) const {
    int i1 = 0;
    while(getArgument(i1) < x1)
        i1++;
    int i2 = 0;
    while(getArgument(i2) < x2)
        i2++;
    std::complex<double> C1 = (getValue(i1) - getValue(i1-1)) / (getArgument(i1) - getArgument(i1-1)) * (x1 - getArgument(i1-1)) + getValue(i1-1);
    std::complex<double> C2 = (getValue(i2) - getValue(i2-1)) / (getArgument(i2) - getArgument(i2-1)) * (x2 - getArgument(i2-1)) + getValue(i2-1);
    std::complex<double> D1 = (getValue(i1) - getValue(i1-1)) / (getArgument(i1) - getArgument(i1-1));
    std::complex<double> D2 = (getValue(i2) - getValue(i2-1)) / (getArgument(i2) - getArgument(i2-1));
    return interpolate(x1, x2, C1, C2, D1, D2);
}

math::CFunction CFunction::interpolate(double x1, double x2,  std::complex<double> C1, std::complex<double> C2 ,std::complex<double> D1, std::complex<double> D2) const {
    assert(x1 < x2);

    // find coefficients first
    std::complex<double> a = -((2.0*C1 - 2.0*C2 - D1*x1 - D2*x1 + D1*x2 + D2*x2)/std::pow(x1 - x2,3));
    std::complex<double> b = -((-3.0*C1*x1 + 3.0*C2*x1 + D1*std::pow(x1,2) + 2.0*D2*std::pow(x1,2) - 3.0*C1*x2 + 3.0*C2*x2 +
                                D1*x1*x2 - D2*x1*x2 - 2.0*D1*std::pow(x2,2) - D2*std::pow(x2,2))/std::pow(x1 - x2,3));
    std::complex<double> c = -((-(D2*x1) + D1*x2)/(x1 - x2)) - (3.0*x1*x2*
                                                                (2.0*C1 - 2.0*C2 - D1*x1 - D2*x1 + D1*x2 + D2*x2))/std::pow(x1 - x2,3);
    std::complex<double> d = -((-(C2*std::pow(x1,3)) + 3.0*C2*std::pow(x1,2)*x2 + D2*std::pow(x1,3)*x2 - 3.0*C1*x1*std::pow(x2,2) +
                                D1*std::pow(x1,2)*std::pow(x2,2) - D2*std::pow(x1,2)*std::pow(x2,2) + C1*std::pow(x2,3) -
                                D1*x1*std::pow(x2,3))/std::pow(x1 - x2,3));

    math::CFunction interpolated(getSize());
    for(int n = 0; n < getSize(); n++) {
        double omega = getArgument(n);
        if( (omega > x1) && (omega < x2) ) {
            //interpolate
            interpolated.set(n, omega, a * std::pow(omega, 3) + b * std::pow(omega,2) + c * omega + d);
        } else {
            // copy
            interpolated.set(n, omega, getValue(n));
        }
    }
    return interpolated;
}

math::CFunction CFunction::transformPH() const {
    int N = getSize();
    if(N==0)
        return math::CFunction();
    math::CFunction result(N);
    for(int i = 0; i < N; i++)
        result.set(i, getArgument(i), std::complex<double>(-getValueReal(N-1-i),getValueImag(N-1-i)));
    return result;
}


void CFunction::createLogGrid(double maxFrequency, double lambda, int n) {
    assert(lambda > 0);
    assert(maxFrequency > 0);
    assert(n > 0);
    resize(n);
    for(int i = 0; i < n/2; i++)
        set(i, -maxFrequency * std::pow(lambda, -i),0);
    if(n % 2 == 1)
        set(n/2,0,0);
    for(int i = 0; i < n/2; i++)
        set(n-1-i, maxFrequency * std::pow(lambda, -i),0);
}

std::complex<double> CFunction::integrate(double x_1, double x_2) const {
    std::complex<double> integral = 0.0;

    int nmax = m_length;


    if(x_1 < m_arguments[0])
    {
        std::string message = "Left boundary outside value interval.";
        std::cerr << x_1 << std::endl;
        throw std::exception();
    }

    if(x_2 > m_arguments[nmax-1])
    {
        std::string message = "Right boundary outside value interval.";
        std::cerr << x_2 << std::endl;

        throw std::exception();
    }

    int i_1 = 0;
    while( (i_1 < nmax) && (m_arguments[i_1] < x_1) )
    {
        i_1++;
    }
    if(i_1 == nmax)
    {
        std::string message = "Left boundary outside value interval.";
        std::cerr << x_1 << std::endl;

        throw std::exception();
    }

     //   std::cerr << i_1 << " " << m_arguments[i_1] << " " << m_values[i_1]  << std::endl;

    int i_2 = nmax-1;
    while( (i_2 >= 0) && (m_arguments[i_2] > x_2) )
    {
        i_2--;
    }
    if(i_2 == -1)
    {
        std::string message = "Right boundary outside value interval.";
        std::cerr << x_2 << std::endl;

        throw std::exception();
    }

    //    std::cerr << i_2 << " " << m_arguments[i_2] << " " << m_values[i_2]  << std::endl;

    if(i_1 < i_2)
    {
        //integrate
        std::complex<double> left;
        std::complex<double> right;
        double darg;
        for( int i = i_1; i < i_2; i++)
        {
            //		std::cerr << i << std::endl;
            left = m_values[i];
            right = m_values[i + 1];
            darg = m_arguments[i+1] - m_arguments[i];
            integral += darg * (right + left)/2.0;
        }
    }
    else if(i_1 > i_2)
    {
        //	std::cerr << "No discretization point in interval" << std::endl;

        double arg_left = m_arguments[i_1-1];
        double arg_right = m_arguments[i_1];


        std::complex<double> left = m_values[i_1-1];
        std::complex<double> right = m_values[i_1];

        //		std::cerr << "Interpolation(left) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

        std::complex<double> value_x_1 = (right - left) / ( arg_right - arg_left ) * (x_1 - arg_left) + left;
        //		std::cerr << "interpolated left value: (" << x_1 << "," << value_x_1 << ")" << std::endl;

        arg_left = m_arguments[i_2];
        arg_right = m_arguments[i_2+1];


        left = m_values[i_2];
        right = m_values[i_2+1];

        //		std::cerr << "Interpolation(right) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

        std::complex<double> value_x_2 = (right - left) / ( arg_right - arg_left ) * (x_2 - arg_left) + left;
        //		std::cerr << "interpolated right value: (" << x_2 << "," << value_x_2 << ")" << std::endl;

        integral = ((value_x_2 + value_x_1) * (x_2 - x_1))/2.0;
        //		std::cerr << "integral = " << integral << std::endl;

        return integral;

    }

    //	std::cerr << "integral_middle = " << integral << std::endl;

    if(x_1 < m_arguments[i_1])
    {
        //		std::cerr << "add left" << std::endl;
        double arg_left = m_arguments[i_1-1];
        double arg_right = m_arguments[i_1];


        std::complex<double> left = m_values[i_1-1];
        std::complex<double> right = m_values[i_1];

        //		std::cerr << "Interpolation(left) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

        std::complex<double> value_x_1 = (right - left) / ( arg_right - arg_left ) * (x_1 - arg_left) + left;
        //		std::cerr << "interpolated left value: " << value_x_1 << std::endl;

        std::complex<double> integralL = ((value_x_1  + right) * (arg_right - x_1))/2.0;
        //	std::cerr << "Left contribution = " << integralL << std::endl;
        integral += integralL;




    }

    if(x_2 > m_arguments[i_2])
    {
        //		std::cerr << "add right" << std::endl;

        double arg_left = m_arguments[i_2];
        double arg_right = m_arguments[i_2+1];


        std::complex<double> left = m_values[i_2];
        std::complex<double> right = m_values[i_2+1];

        //		std::cerr << "Interpolation(right) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

        std::complex<double> value_x_2 = (right - left) / ( arg_right - arg_left ) * (x_2 - arg_left) + left;
        //		std::cerr << "interpolated left value: " << value_x_2 << std::endl;

        std::complex<double> integralR = ((value_x_2 + left) * (x_2 - arg_left) )/2.0;
        //		std::cerr << "Right contribution = " << integralR << std::endl;
        integral += integralR;

    }
    //	std::cerr << "final integral = " << integral << std::endl;

    return integral;
}


std::complex<double> CFunction::integrate(double x_1, double x_2, tIntegrandFunction g, void * parameters) const {
    std::complex<double> integral = 0.0;

    int nmax = m_length;


    if(x_1 < m_arguments[0])
    {
        std::string message = "Left boundary outside value interval.";
        std::cerr << x_1 << std::endl;
        throw std::exception();
    }

    if(x_2 > m_arguments[nmax-1])
    {
        std::string message = "Right boundary outside value interval.";
        std::cerr << x_2 << std::endl;

        throw std::exception();
    }

    int i_1 = 0;
    while( (i_1 < nmax) && (m_arguments[i_1] < x_1) )
    {
        i_1++;
    }
    if(i_1 == nmax)
    {
        std::string message = "Left boundary outside value interval.";
        std::cerr << x_1 << std::endl;

        throw std::exception();
    }

    //	std::cerr << i_1 << " " << m_arguments[i_1] << " " << m_values[i_1]  << std::endl;

    int i_2 = nmax-1;
    while( (i_2 >= 0) && (m_arguments[i_2] > x_2) )
    {
        i_2--;
    }
    if(i_2 == -1)
    {
        std::string message = "Right boundary outside value interval.";
        std::cerr << x_2 << std::endl;

        throw std::exception();
    }

    //	std::cerr << i_2 << " " << m_arguments[i_2] << " " << m_values[i_2]  << std::endl;

    if(i_1 < i_2)
    {
        //integrate
        std::complex<double> left;
        std::complex<double> right;
        double darg;
        for( int i = i_1; i < i_2; i++)
        {
            //		std::cerr << i << std::endl;
            left = m_values[i] * (*g)(m_arguments[i], m_values[i], parameters);
            right = m_values[i + 1] * (*g)(m_arguments[i+1], m_values[i+1], parameters);
            darg = m_arguments[i+1] - m_arguments[i];
            integral += darg * (right + left)/2.0;
        }
    }
    else if(i_1 > i_2)
    {
        //	std::cerr << "No discretization point in interval" << std::endl;

        double arg_left = m_arguments[i_1-1];
        double arg_right = m_arguments[i_1];


        std::complex<double> left = m_values[i_1-1] * (*g)(m_arguments[i_1-1], m_values[i_1-1], parameters);
        std::complex<double> right = m_values[i_1] * (*g)(m_arguments[i_1], m_values[i_1], parameters);

        //		std::cerr << "Interpolation(left) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

        std::complex<double> value_x_1 = (right - left) / ( arg_right - arg_left ) * (x_1 - arg_left) + left;
        //		std::cerr << "interpolated left value: (" << x_1 << "," << value_x_1 << ")" << std::endl;

        arg_left = m_arguments[i_2];
        arg_right = m_arguments[i_2+1];


        left = m_values[i_2] * (*g)(m_arguments[i_2], m_values[i_2], parameters);
        right = m_values[i_2+1] * (*g)(m_arguments[i_2+1], m_values[i_2+1], parameters);

        //		std::cerr << "Interpolation(right) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

        std::complex<double> value_x_2 = (right - left) / ( arg_right - arg_left ) * (x_2 - arg_left) + left;
        //		std::cerr << "interpolated right value: (" << x_2 << "," << value_x_2 << ")" << std::endl;

        integral = (( value_x_2 + value_x_1) * (x_2 - x_1))/2.0;
        //		std::cerr << "integral = " << integral << std::endl;

        return integral;

    }

    //	std::cerr << "integral_middle = " << integral << std::endl;

    if(x_1 < m_arguments[i_1])
    {
        //		std::cerr << "add left" << std::endl;
        double arg_left = m_arguments[i_1-1];
        double arg_right = m_arguments[i_1];


        std::complex<double> left = m_values[i_1-1] * (*g)(m_arguments[i_1-1], m_values[i_1-1], parameters);
        std::complex<double> right = m_values[i_1] * (*g)(m_arguments[i_1], m_values[i_1], parameters);

        //		std::cerr << "Interpolation(left) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

        std::complex<double> value_x_1 = (right - left) / ( arg_right - arg_left ) * (x_1 - arg_left) + left;
        //		std::cerr << "interpolated left value: " << value_x_1 << std::endl;

        std::complex<double> integralL = ((value_x_1  + right ) * (arg_right - x_1))/2.0;
        //	std::cerr << "Left contribution = " << integralL << std::endl;
        integral += integralL;




    }

    if(x_2 > m_arguments[i_2])
    {
        //		std::cerr << "add right" << std::endl;

        double arg_left = m_arguments[i_2];
        double arg_right = m_arguments[i_2+1];


        std::complex<double> left = m_values[i_2] * (*g)(m_arguments[i_2], m_values[i_2], parameters);
        std::complex<double> right = m_values[i_2+1] * (*g)(m_arguments[i_2+1], m_values[i_2+1], parameters);

        //		std::cerr << "Interpolation(right) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

        std::complex<double> value_x_2 = (right - left) / ( arg_right - arg_left ) * (x_2 - arg_left) + left;
        //		std::cerr << "interpolated left value: " << value_x_2 << std::endl;

        std::complex<double> integralR = ((value_x_2  + left ) * (x_2 - arg_left) )/2.0;
        //		std::cerr << "Right contribution = " << integralR << std::endl;
        integral += integralR;

    }
    //	std::cerr << "final integral = " << integral << std::endl;

    return integral;
}


double similarity(const math::CFunction& a, const math::CFunction& b) {
  // assert(a.getSize() == b.getSize());
  if(a.getSize() != b.getSize())
    return 0.0;
  double intDiff = 0;
    double aInt = 0;
    for(int n = 0; n < a.getSize()-1; n++) {
        aInt += a.getValue(n).imag() * (a.getArgument(n+1) - a.getArgument(n));
        intDiff += std::abs(a.getValue(n) - b.getValue(n)) * (a.getArgument(n+1) - a.getArgument(n));
    }
    if(std::abs(aInt) < intDiff)
        return 0;
    return 1.0 - intDiff / std::abs(aInt);
}


}
