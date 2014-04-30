#include "pfunction.h"
#include "cfunction.h"
#include "exceptions.h"
#include <assert.h>

namespace math {
PFunction::PFunction()
    : m_arguments(0)
    , m_values(0)
    , m_length(0)
{
}

PFunction::PFunction(int length)
    : m_arguments(new mpfr::mpreal[length])
    , m_values(new mpfr::mpreal[length])
    , m_length(length)
{
}

PFunction::PFunction(const PFunction& orig)
    : m_length(orig.m_length)
{
    if(m_length > 0) {
        m_arguments = new mpfr::mpreal[m_length];
        m_values = new mpfr::mpreal[m_length];
        for(int i = 0; i < m_length; i++) {
            m_arguments[i] = orig.m_arguments[i];
            m_values[i] = orig.m_values[i];
        }
    } else {
        m_arguments = 0;
        m_values = 0;
    }
}

PFunction& PFunction::operator=(const PFunction& orig) {
    m_length = orig.m_length;
    if(m_arguments != 0)
        delete [] m_arguments;
    if(m_values != 0)
        delete [] m_values;

    if(m_length > 0) {
        m_arguments = new mpfr::mpreal[m_length];
        m_values = new mpfr::mpreal[m_length];
        for(int i = 0; i < m_length; i++) {
            m_arguments[i] = orig.m_arguments[i];
            m_values[i] = orig.m_values[i];
        }
    } else {
        m_arguments = 0;
        m_values = 0;
    }
    return *this;
}

PFunction::PFunction(const CFunction& orig)
    : m_length(orig.m_length)
{
    if(m_length > 0) {
        m_arguments = new mpfr::mpreal[m_length];
        m_values = new mpfr::mpreal[m_length];
        for(int i = 0; i < m_length; i++) {
            m_arguments[i] = orig.m_arguments[i];
            m_values[i] = orig.m_values[i].imag();
        }
    } else {
        m_arguments = 0;
        m_values = 0;
    }
}

PFunction& PFunction::operator=(const CFunction& orig) {
    m_length = orig.m_length;
    if(m_arguments != 0)
        delete [] m_arguments;
    if(m_values != 0)
        delete [] m_values;

    if(m_length > 0) {
        m_arguments = new mpfr::mpreal[m_length];
        m_values = new mpfr::mpreal[m_length];
        for(int i = 0; i < m_length; i++) {
            m_arguments[i] = orig.m_arguments[i];
            m_values[i] = orig.m_values[i].imag();
        }
    } else {
        m_arguments = 0;
        m_values = 0;
    }

    return *this;
}

PFunction::~PFunction() {
    if(m_arguments != 0)
        delete [] m_arguments;
    if(m_values != 0)
        delete [] m_values;
}

void PFunction::resize(int length) {
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
        m_arguments = new mpfr::mpreal[m_length];
        m_values = new mpfr::mpreal[m_length];
    }
}

void PFunction::set(int index, mpfr::mpreal argument, mpfr::mpreal value) {
    assert(index >= 0);
    assert(index < m_length);
    m_arguments[index] = argument;
    m_values[index] = value;
}

void PFunction::setArgument(int index, mpfr::mpreal argument) {
    assert(index >= 0);
    assert(index < m_length);
    m_arguments[index] = argument;
}

void PFunction::setValue(int index, mpfr::mpreal value) {
    assert(index >= 0);
    assert(index < m_length);
    m_values[index] = value;
}

mpfr::mpreal PFunction::getArgument(int index) const {
    assert(index >= 0);
    assert(index < m_length);
    return m_arguments[index];
}

mpfr::mpreal PFunction::getValue(int index) const {
    assert(index >= 0);
    assert(index < m_length);
    return m_values[index];
}


mpfr::mpreal PFunction::integrate_arg(mpfr::mpreal x_1, mpfr::mpreal x_2) const
{
    mpfr::mpreal integral = 0.0;

    int nmax = m_length;


    if(x_1 < m_arguments[0])
    {
        std::string message = "Left boundary outside value interval.";
        std::cerr << x_1 << std::endl;
        throw error::math::OutOfBounds(message.c_str());
    }

    if(x_2 > m_arguments[nmax-1])
    {
        std::string message = "Right boundary outside value interval.";
        std::cerr << x_2 << std::endl;

        throw error::math::OutOfBounds(message.c_str());
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

        throw error::math::OutOfBounds(message.c_str());
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

        throw error::math::OutOfBounds(message.c_str());
    }

    //	std::cerr << i_2 << " " << m_arguments[i_2] << " " << m_values[i_2]  << std::endl;

    if(i_1 < i_2)
    {
        //integrate
        mpfr::mpreal left;
        mpfr::mpreal right;
        mpfr::mpreal darg;
        for( int i = i_1; i < i_2; i++)
        {
            //		std::cerr << i << std::endl;
            left = m_values[i];
            right = m_values[i + 1];
            darg = m_arguments[i+1] - m_arguments[i];
            integral += darg * (m_arguments[i+1] * right + m_arguments[i] * left)/2.0;
        }
    }
    else if(i_1 > i_2)
    {
        //	std::cerr << "No discretization point in interval" << std::endl;

        mpfr::mpreal arg_left = m_arguments[i_1-1];
        mpfr::mpreal arg_right = m_arguments[i_1];


        mpfr::mpreal left = m_values[i_1-1];
        mpfr::mpreal right = m_values[i_1];

        //		std::cerr << "Interpolation(left) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

        mpfr::mpreal value_x_1 = (right - left) / ( arg_right - arg_left ) * (x_1 - arg_left) + left;
        //		std::cerr << "interpolated left value: (" << x_1 << "," << value_x_1 << ")" << std::endl;

        arg_left = m_arguments[i_2];
        arg_right = m_arguments[i_2+1];


        left = m_values[i_2];
        right = m_values[i_2+1];

        //		std::cerr << "Interpolation(right) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

        mpfr::mpreal value_x_2 = (right - left) / ( arg_right - arg_left ) * (x_2 - arg_left) + left;
        //		std::cerr << "interpolated right value: (" << x_2 << "," << value_x_2 << ")" << std::endl;

        integral = ((x_2 * value_x_2 + x_1 * value_x_1) * (x_2 - x_1))/2.0;
        //		std::cerr << "integral = " << integral << std::endl;

        return integral;

    }

    //	std::cerr << "integral_middle = " << integral << std::endl;

    if(x_1 < m_arguments[i_1])
    {
        //		std::cerr << "add left" << std::endl;
        mpfr::mpreal arg_left = m_arguments[i_1-1];
        mpfr::mpreal arg_right = m_arguments[i_1];


        mpfr::mpreal left = m_values[i_1-1];
        mpfr::mpreal right = m_values[i_1];

        //		std::cerr << "Interpolation(left) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

        mpfr::mpreal value_x_1 = (right - left) / ( arg_right - arg_left ) * (x_1 - arg_left) + left;
        //		std::cerr << "interpolated left value: " << value_x_1 << std::endl;

        mpfr::mpreal integralL = ((value_x_1 * x_1 + right * arg_right) * (arg_right - x_1))/2.0;
        //	std::cerr << "Left contribution = " << integralL << std::endl;
        integral += integralL;




    }

    if(x_2 > m_arguments[i_2])
    {
        //		std::cerr << "add right" << std::endl;

        mpfr::mpreal arg_left = m_arguments[i_2];
        mpfr::mpreal arg_right = m_arguments[i_2+1];


        mpfr::mpreal left = m_values[i_2];
        mpfr::mpreal right = m_values[i_2+1];

        //		std::cerr << "Interpolation(right) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << std::endl;

        mpfr::mpreal value_x_2 = (right - left) / ( arg_right - arg_left ) * (x_2 - arg_left) + left;
        //		std::cerr << "interpolated left value: " << value_x_2 << std::endl;

        mpfr::mpreal integralR = ((value_x_2 * x_2 + left * arg_left) * (x_2 - arg_left) )/2.0;
        //		std::cerr << "Right contribution = " << integralR << std::endl;
        integral += integralR;

    }
    //	std::cerr << "final integral = " << integral << std::endl;

    return integral;
}

mpfr::mpreal PFunction::integrate(mpfr::mpreal x_1, mpfr::mpreal x_2) const
{
    mpfr::mpreal integral = 0.0;

    int nmax = m_length;


    if(x_1 < m_arguments[0])
    {
        std::string message = "Left boundary outside value interval.";
        std::cerr << x_1 << std::endl;
        std::cerr << m_arguments[0] << std::endl;
        throw error::math::OutOfBounds(message);
    }

    if(x_2 > m_arguments[nmax-1])
    {
        std::string message = "Left boundary outside value interval.";
        std::cerr << x_2 << std::endl;
        throw error::math::OutOfBounds(message);
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
        throw error::math::OutOfBounds(message);
    }


    int i_2 = nmax-1;
    while( (i_2 >= 0) && (m_arguments[i_2] > x_2) )
    {
        i_2--;
    }
    if(i_2 == -1)
    {
        std::string message = "Right boundary outside value interval.";
        std::cerr << x_2 << std::endl;

        throw error::math::OutOfBounds(message);
    }

    if(i_1 < i_2)
    {
        //integrate
        mpfr::mpreal left;
        mpfr::mpreal right;
        mpfr::mpreal darg;
        for( int i = i_1; i < i_2; i++)
        {
            //		cout << i << endl;
            left = m_values[i];
            right = m_values[i + 1];
            darg = m_arguments[i+1] - m_arguments[i];
            integral += darg * (right + left)/2.0;
        }
    }
    else if(i_1 > i_2)
    {
        //	cout << "No discretization point in interval" << endl;

        mpfr::mpreal arg_left = m_arguments[i_1-1];
        mpfr::mpreal arg_right = m_arguments[i_1];


        mpfr::mpreal left = m_values[i_1-1];
        mpfr::mpreal right = m_values[i_1];

        //		cout << "Interpolation(left) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << endl;

        mpfr::mpreal value_x_1 = (right - left) / ( arg_right - arg_left ) * (x_1 - arg_left) + left;
        //	cout << "interpolated left value: (" << x_1 << "," << value_x_1 << ")" << endl;

        arg_left = m_arguments[i_2];
        arg_right = m_arguments[i_2+1];


        left = m_values[i_2];
        right = m_values[i_2+1];

        //		cout << "Interpolation(right) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << endl;

        mpfr::mpreal value_x_2 = (right - left) / ( arg_right - arg_left ) * (x_2 - arg_left) + left;
        //		cout << "interpolated right value: (" << x_2 << "," << value_x_2 << ")" << endl;

        integral = (value_x_2 + value_x_1)/2.0 * (x_2 - x_1);
        //		cout << "integral = " << integral << endl;

        return integral;

    }

    //	cout << "integral = " << integral << endl;

    if(x_1 < m_arguments[i_1])
    {
        //		cout << "add left" << endl;
        mpfr::mpreal arg_left = m_arguments[i_1-1];
        mpfr::mpreal arg_right = m_arguments[i_1];


        mpfr::mpreal left = m_values[i_1-1];
        mpfr::mpreal right = m_values[i_1];

        //		cout << "Interpolation(left) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << endl;

        mpfr::mpreal value_x_1 = (right - left) / ( arg_right - arg_left ) * (x_1 - arg_left) + left;
        //	cout << "interpolated left value: " << value_x_1 << endl;

        mpfr::mpreal integralL = (value_x_1 + right)/2.0 * (arg_right - x_1);
        //	cout << "Left contribution = " << integralL << endl;
        integral += integralL;




    }

    if(x_2 > m_arguments[i_2])
    {
        //	cout << "add right" << endl;

        mpfr::mpreal arg_left = m_arguments[i_2];
        mpfr::mpreal arg_right = m_arguments[i_2+1];


        mpfr::mpreal left = m_values[i_2];
        mpfr::mpreal right = m_values[i_2+1];

        //	cout << "Interpolation(right) between (" << arg_left << "," << left << ") and (" << arg_right << "," << right << ")" << endl;

        mpfr::mpreal value_x_2 = (right - left) / ( arg_right - arg_left ) * (x_2 - arg_left) + left;
        //	cout << "interpolated left value: " << value_x_2 << endl;

        mpfr::mpreal integralR = (value_x_2 + left)/2.0 * (x_2 - arg_left);
        //	cout << "Right contribution = " << integralR << endl;
        integral += integralR;

    }
    //	cout << "final integral = " << integral << endl;

    return integral;
}

std::ostream& operator<<(std::ostream& stream, const math::PFunction& p) {
    for(int i = 0; i < p.m_length; i++)
        stream << p.m_arguments[i] << " " << p.m_values[i] << std::endl;
    return stream;
}


}
