#ifndef FUNCTION_H
#define FUNCTION_H

/**
 * @file
 *
 * @ingroup math
 *
 * @brief Discretized complex valued function for real arguments
 */

/*
 * Copyright (c) 2014 Christoph Schuette.
 *
 * The license and distribution terms for this file may be
 * found in the file LICENSE in this distribution
 */

#include <complex>
#include <vector>

namespace mpi {
    class OpenMPI;
}
namespace math {
    class PFunction;

    /**
     * @brief Discretized complex-valued function for real arguments
     * 
     * The arguments should be added in ascending order. It is the programmer's
     * job to ensure this.
     */
    class CFunction {
    public:
        typedef std::complex<double> (*tIntegrandFunction)(double, std::complex<double>, void*);

        inline static std::complex<double> one(double omega, std::complex<double> fomega, void* parameters) {
            return 1.0;
        }

        inline static std::complex<double> arg(double omega, std::complex<double> fomega, void* parameters) {
            return omega;
        }

    public:
        /**
         * @brief empty, default constructor
         */
        CFunction();
        /**
         * @brief allocates an instance for <i>length</i> data points.
         * @param length # of data points
         */
        CFunction(int length);
        /**
         * @brief copy constructor
         * @param[in] orig reference to original instance
         */
        CFunction(const CFunction& orig);
        
        /**
         * @brief operator=
         * @param[in] orig reference to original instance
         * @return reference to self
         */
        CFunction& operator=(const CFunction& orig);
        
        /**
         * @brief destructor frees all previously allocated memory.
         */
        virtual ~CFunction();

        /**
         * @brief changes the size of the storage container to hold <i>length</i>
         *        argument-value pairs.
         * 
         * All previous data is discarded.
         * 
         * @param[in] length new size of the storage container.
         */
        void resize(int length);

        /**
         * @brief returns current size of the storage container
         * @return size of the storage container
         */
        int getSize() const {
            return m_length;
        }

        /**
         * @brief sets the datapoint with indicated index.
         * 
         * The function sets argument of the argument-value pair with index <i>index</i> 
         * to <i>argument</i> and value <i>value</i>.
         * 
         * @param[in] index index of name-value pair
         * @param[in] argument new argument
         * @param value new value
         */
        void set(int index, double argument, std::complex<double> value);
        
        /**
         * @brief sets argument of the argument-value pair with index <i>index</i> 
         * to <i>argument</i>
         * @param[in] index index of name-value pair
         * @param[in] argument new argument
         */
        void setArgument(int index, double argument);
        
        /**
         * @brief sets value of the argument-value pair with index <i>index</i> 
         * to <i>value</i>
         * @param[in] index index of argument-value pair
         * @param[in] value new value
         */
        void setValue(int index, std::complex<double> value);
        
        /**
         * @brief returns the argument of the pair with index <i>index</i>.
         * @param[in] index index of name-value pair
         * @return argument of the argument-value pair with index <i>index</i>.
         */
        double getArgument(int index) const;
        
        /**
         * @brief returns the value of the pair with index <i>index</i>.
         * @param[in] index index of name-value pair
         * @return value of the argument-value pair with index <i>index</i>.
         */
        std::complex<double> getValue(int index) const;
        
        /**
         * @brief returns the value's real part of the pair with index <i>index</i>.
         * @param[in] index index of name-value pair
         * @return value's real part of the argument-value pair with index <i>index</i>.
         */        
        double getValueReal(int index) const;
        
        /**
         * @brief returns the value's imaginary part of the pair with index <i>index</i>.
         * @param[in] index index of name-value pair
         * @return value's imaginary part of the argument-value pair with index <i>index</i>.
         */
        double getValueImag(int index) const;

        //IO
        /**
         * @brief writes the data of the instance as a text file named <i>filename</i>.
         * 
         * Every line of the file contains an argument-value pair with a space as
         * a separator.
         * 
         * @param[in] filename filename for the written file
         */
        void write(std::string filename) const;
        
        /**
         * @brief reads data from the text file named <i>filename</i>.
         * 
         * Every line of the file should contain an argument-value pair with a space as
         * a separator. The size of the storage container is adjusted. All previous
         * data is discarded.
         * 
         * @param[in] filename filename for the written file
         */
        void read(std::string filename);
        
        /**
         * @brief writes the data of the instance as a binary file named <i>filename</i>.
         * 
         * The first integer in the file denotes the # of name-value pairs.
         * 
         * @param filename filename for the written file
         */
        void writeBinary(std::string filename) const;
        
        /**
         * @brief reads data to the instance as a binary file named <i>filename</i>.
         * 
         * The first integer in the file should denote the # of name-value pairs. The
         * size of the storage container is adjusted. All previous
         * data is discarded.
         * 
         * @param filename filename for the written file
         */
        void readBinary(std::string filename);
        
        /**
         * @brief dumps the data of the instance to std::cout.
         */
        void dump() const;

        // operators
        /**
         * @brief multiplies all arguments and values by the factor <i>factor</i>.
         * @param factor scaling factor
         */
        void rescale(double factor);
        
        /**
         * @brief compares the current instance to the instance <i>rhs</i>
         * 
         * This function returns true if the # of argument-value pairs and all
         * arguments and values of the two instances coincide.
         * 
         * @param[in] rhs instance to compare to
         * @return true if the two are same, false otherwise.
         */
        bool operator==(const CFunction& rhs) const;

        //convenience
        /**
         * @brief creates a logarithmic grid in argument space about 0.
         * 
         * This function creates a function with \f$n\f$ argument-value pairs
         * with arguments
         * \f$-A \lambda^0, -A \lambda^1, \dots, -A \lambda^{-n/2},A \lambda^{-n/2},\dots,A \lambda^1,A \lambda^0\f$.
         * If \f$n%2==1\f$ an additional datapoint with argument 0 is created. The
         * values are initialized with 0's.
         * 
         * @param[in] maxFrequency maximum frequency
         * @param[in] lambda discretization parameter
         * @param[in] n # of datapoints
         */
        void createLogGrid(double maxFrequency, double lambda, int n);
 
        math::CFunction transformPH() const;

        // interpolation
        /**
         * @brief finds the linear interpolation of the discretized function at <i>x</i>.
         * 
         * If <i>x</i> lies outside the range of arguments, the function throws an
         * error::math::OutOfBounds exception.
         * 
         * @param[in] x interpolation point
         * @return interpolated value
         */
        std::complex<double> interpolate(double x);
        math::CFunction interpolate(double x1, double x2) const;
        math::CFunction interpolate(double x1, double x2, std::complex<double> C1, std::complex<double> C2, std::complex<double> D1, std::complex<double> D2) const;

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
