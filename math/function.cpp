#include "function.h"
#include <string.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>

namespace math {
Function::Function()
    :m_arguments(0)
    ,m_values(0)
    ,m_length(0)
{
}

Function::Function(int length)
    : m_arguments(new double[length])
    , m_values(new double[length])
    , m_length(length)
{

}

Function::Function(const Function& orig)
    : m_length(orig.m_length)
{
    if(orig.m_arguments != 0) {
        m_arguments = new double[m_length];
        memcpy(m_arguments, orig.m_arguments, sizeof(double) * m_length);
    }
    if(orig.m_values != 0) {
        m_values = new double[m_length];
        memcpy(m_values, orig.m_values, sizeof(double) * m_length);
    }
}

Function& Function::operator=(const Function& orig) {
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
        m_values = new double[m_length];
        memcpy(m_values, orig.m_values, sizeof(double) * m_length);
    }
    return *this;
}

Function::~Function() {
    if(m_arguments != 0)
        delete [] m_arguments;
    if(m_values != 0)
        delete [] m_values;
}

void Function::resize(int length) {
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
        m_values = new double[m_length];
    }
}

void Function::set(int index, double argument, double value) {
    assert(index >= 0);
    assert(index < m_length);
    m_arguments[index] = argument;
    m_values[index] = value;
}

void Function::setArgument(int index, double argument) {
    assert(index >= 0);
    assert(index < m_length);
    m_arguments[index] = argument;
}

void Function::setValue(int index, double value) {
    assert(index >= 0);
    assert(index < m_length);
    m_values[index] = value;
}

double Function::getArgument(int index) const {
    assert(index >= 0);
    assert(index < m_length);
    return m_arguments[index];
}
double Function::getValue(int index) const {
    assert(index >= 0);
    assert(index < m_length);
    return m_values[index];
}

double Function::interpolate(double x) {
    if(x < m_arguments[0])
    {
        std::string message = "Left boundary outside value interval.";
        std::cerr << x << std::endl;
        throw std::exception();
    }

    if(x > m_arguments[m_length-1])
    {
        std::string message = "Right boundary outside value interval.";
        std::cerr << x << std::endl;

        throw std::exception();
    }

    int i = 0;
    while( (i < m_length) && (m_arguments[i] <= x) )
    {
        i++;
    }
    if(i == m_length)
    {
        std::string message = "Left boundary outside value interval.";
        std::cerr << x << std::endl;

        throw std::exception();
    }

    //interpolate between i and i-1

    double arg_left = m_arguments[i-1];
    double arg_right = m_arguments[i];

    double left = m_values[i-1];
    double right = m_values[i];

    //		std::cerr << "Interpolation(right) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

    double value = (right - left) / ( arg_right - arg_left ) * (x - arg_left) + left;
    //		std::cerr << "interpolated right value: (" << x_2 << "," << value_x_2 << ")" << std::endl;
    return value;

}


//IO
void Function::read(std::string filename) {
    std::ifstream inFile(filename.c_str());
    if(!inFile.is_open()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        throw std::exception();
    }
    std::vector<double> data[2];
    double f, v;
    while (inFile >> f >> v)
    {
        data[0].push_back(f);
        data[1].push_back(v);
    }
    //    std::cout << "Read " << data[0].size() << " entries." << std::endl;
    resize(data[0].size());
    if(data[0].size() == 0)
        std::cout << "WARNING: Function contains no data." << std::endl;
    for(int i = 0; i < data[0].size(); i++) {
        m_arguments[i] = data[0][i];
        m_values[i] = data[1][i];
    }
}


void Function::write(std::string filename) {
    std::ofstream outFile(filename.c_str());
    if(!outFile.is_open()) {
        std::cerr << "Cannot open file " << filename << std::endl;
        throw std::exception();
    }
    for(int i = 0; i < m_length; i++)
        outFile << m_arguments[i] << " " << m_values[i] << std::endl;
    outFile.close();
}
}
