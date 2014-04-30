#include "integral.h"
#include <string.h>
#include <iostream>
#include "../mpi/openmpi.h"
#include <gsl/gsl_errno.h>
#include <stdarg.h>     /* va_list, va_start, va_arg, va_end */

//#define DEBUG

#ifdef DEBUG
#include "../debug/timer.h"
#endif

namespace math {

void Integral::integrationErrorHandler(const char * x, const char * y, int a, int b) {
    std::cout << x << " " << y << " " << a << " " << b << std::endl;
}


Integral::Integral()
    : m_bounds(0)
    , m_nBounds(0)
    , m_kernelReal(0)
    , m_kernelImag(0)
    , m_integrationWorkspace(gsl_integration_workspace_alloc(1000))
    , m_n(0)
    , m_integrationParameters(0)
    , m_real(true)
    , m_imag(true)
{
    m_gslIntegrationFunctionReal.function = &Integral::integrationHelper_real;
    m_gslIntegrationFunctionReal.params = this;
    m_gslIntegrationFunctionImag.function = &Integral::integrationHelper_imag;
    m_gslIntegrationFunctionImag.params = this;
}

/*
    double * m_bounds;
    int m_nBounds;
    tIntegralKernel m_kernelReal;
    tIntegralKernel m_kernelImag;
    CFunction* m_parameterFunction;
    gsl_function m_gslIntegrationFunctionReal;
    gsl_function m_gslIntegrationFunctionImag;
    Accuracy m_accuracy;
    gsl_integration_workspace* m_integrationWorkspace;
    int m_n;
    void * m_integrationParameters;
  */

Integral::Integral(const Integral& orig) {
    m_nBounds = orig.m_nBounds;
    m_bounds = new double[m_nBounds];
    memcpy(m_bounds, orig.m_bounds, sizeof(double) * m_nBounds);
    m_kernelReal = orig.m_kernelReal;
    m_kernelImag = orig.m_kernelImag;
    m_parameterFunction = orig.m_parameterFunction;
    m_gslIntegrationFunctionReal.function = &Integral::integrationHelper_real;
    m_gslIntegrationFunctionReal.params = this;
    m_gslIntegrationFunctionImag.function = &Integral::integrationHelper_imag;
    m_gslIntegrationFunctionImag.params = this;
    m_accuracy = orig.m_accuracy;
    m_integrationWorkspace = gsl_integration_workspace_alloc(1000);
    m_n = orig.m_n;
    m_integrationParameters = orig.m_integrationParameters;
    m_real = orig.m_real;
    m_imag = orig.m_imag;
}

Integral& Integral::operator=(const Integral& orig) {
    if(m_bounds != 0)
        delete [] m_bounds;
    m_nBounds = orig.m_nBounds;
    m_bounds = new double[m_nBounds];
    memcpy(m_bounds, orig.m_bounds, sizeof(double) * m_nBounds);
    m_kernelReal = orig.m_kernelReal;
    m_kernelImag = orig.m_kernelImag;
    m_parameterFunction = orig.m_parameterFunction;
    m_gslIntegrationFunctionReal.function = &Integral::integrationHelper_real;
    m_gslIntegrationFunctionReal.params = this;
    m_gslIntegrationFunctionImag.function = &Integral::integrationHelper_imag;
    m_gslIntegrationFunctionImag.params = this;
    m_accuracy = orig.m_accuracy;
    m_integrationWorkspace = gsl_integration_workspace_alloc(1000);
    m_n = orig.m_n;
    m_integrationParameters = orig.m_integrationParameters;
    m_real = orig.m_real;
    m_imag = orig.m_imag;
    return *this;
}

Integral::~Integral() {
    gsl_integration_workspace_free(m_integrationWorkspace);
    if(m_bounds != 0)
        delete [] m_bounds;
}

void Integral::setBounds(int n, ...) {
    m_nBounds = n;
    m_bounds = new double [n];
    va_list argptr;
    va_start( argptr, n );
    for(int i = 0; i < n; i++)
        m_bounds[i] = va_arg( argptr, double );
    va_end( argptr );
}

double Integral::integrationHelper_real(double eps, void * integral) {
    Integral* i = static_cast<Integral*>(integral);
    return (*i->m_kernelReal)(eps,i->m_n,i->m_parameterFunction->getArgument(i->m_n),i->m_parameterFunction->getValue(i->m_n),i->m_integrationParameters);
}

double Integral::integrationHelper_imag(double eps, void * integral) {
    Integral* i = static_cast<Integral*>(integral);
    return (*i->m_kernelImag)(eps,i->m_n,i->m_parameterFunction->getArgument(i->m_n),i->m_parameterFunction->getValue(i->m_n),i->m_integrationParameters);
}

void Integral::operator()(CFunction& result) {
    int N = m_parameterFunction->getSize();
    result.resize(N);
    double real;
    double imag;
    double error;
#ifdef DEBUG
    debug::Timer timer;
#endif
    for(int n = 0; n < N; n++) {
        m_n = n;
        gsl_set_error_handler(&Integral::integrationErrorHandler);
        if(m_real)
            gsl_integration_qagp(&m_gslIntegrationFunctionReal, m_bounds, m_nBounds, m_accuracy.absolute, m_accuracy.relative, 1000, m_integrationWorkspace, &real, &error);
        else
            real = 0;
        if(m_imag)
            gsl_integration_qagp(&m_gslIntegrationFunctionImag, m_bounds, m_nBounds, m_accuracy.absolute, m_accuracy.relative, 1000, m_integrationWorkspace, &imag, &error);
        else
            imag = 0;
        result.set(n, m_parameterFunction->getArgument(n), std::complex<double>(real, imag));
#ifdef DEBUG
        std::cout << m_parameterFunction->getArgument(m_n) << " - " << m_parameterFunction->getValue(m_n) << " error: " << error << std::endl;
#endif
    }
#ifdef DEBUG
    std::cout << " * took " << timer.elapsed() << " seconds." << std::endl;
#endif
}

void Integral::operator()(CFunction& result, mpi::OpenMPI& com) {
    int N = m_parameterFunction->getSize()/com.getSize();
    if((m_parameterFunction->getSize() % com.getSize()) > com.getRank())
        N++;
    result.resize(N);
    double real;
    double imag;
    double error;
#ifdef DEBUG
    debug::Timer timer;
#endif
    int n = com.getRank();
    for(int i = 0; i < N; i++) {
        m_n = n;
        gsl_set_error_handler(&Integral::integrationErrorHandler);
        if(m_real)
            gsl_integration_qagp(&m_gslIntegrationFunctionReal, m_bounds, m_nBounds, m_accuracy.absolute, m_accuracy.relative, 1000, m_integrationWorkspace, &real, &error);
        else
            real = 0;
        if(m_imag)
            gsl_integration_qagp(&m_gslIntegrationFunctionImag, m_bounds, m_nBounds, m_accuracy.absolute, m_accuracy.relative, 1000, m_integrationWorkspace, &imag, &error);
        else imag = 0;
        result.set(i, m_parameterFunction->getArgument(n), std::complex<double>(real, imag));
        n += com.getSize();
#ifdef DEBUG
        std::cout << "<" << com.getRank() << "> " << " * " << m_parameterFunction->getArgument(m_n) << " - " << m_parameterFunction->getValue(m_n) << " error: " << error << std::endl;
#endif
    }
#ifdef DEBUG
    std::cout << "<" << com.getRank() << "> " << " * did " << N << " datapoints." << std::endl;
    std::cout << "<" << com.getRank() << "> " << " * took " << timer.elapsed() << " seconds." << std::endl;
#endif

    com.combine(result);

}

}
