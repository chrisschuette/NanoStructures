#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "cfunction.h"

#include <gsl/gsl_integration.h>

namespace mpi {
    class OpenMPI;
}

namespace math {
class Integral
{
    typedef void (*tIntegrationErrorHandler)(const char*, const char*, int, int);
public:
    static void integrationErrorHandler(const char * x, const char * y, int a, int b);
    typedef struct {
        double absolute;
        double relative;
    } Accuracy;
public:
    typedef double (*tIntegralKernel)(double eps, int nOmega, double omega, std::complex<double> fomega, void * parameters);
    Integral();
    Integral(const Integral& orig);
    Integral& operator=(const Integral& orig);
    ~Integral();

    void setBounds(int n, ...);
    void setIntegralKernels(tIntegralKernel real, tIntegralKernel imag) { m_kernelReal = real; m_kernelImag = imag; }
    void setAccuracy(double abs, double rel) { m_accuracy.absolute = abs; m_accuracy.relative = rel; }
    void setParameterFunction(const CFunction& pf) { m_parameterFunction = &pf; }
    void setIntegrationParameters(void * parameters) { m_integrationParameters = parameters; }
    void setIntegrals(bool real, bool imag) { m_real = real; m_imag = imag; }

    void operator()(CFunction& result);
    void operator()(CFunction& result, mpi::OpenMPI& com);
protected:
    static double integrationHelper_real(double eps, void * integral);
    static double integrationHelper_imag(double eps, void * integral);

    bool m_real;
    bool m_imag;

    double * m_bounds;
    int m_nBounds;
    tIntegralKernel m_kernelReal;
    tIntegralKernel m_kernelImag;
    const CFunction* m_parameterFunction;
    gsl_function m_gslIntegrationFunctionReal;
    gsl_function m_gslIntegrationFunctionImag;
    Accuracy m_accuracy;
    gsl_integration_workspace* m_integrationWorkspace;
    int m_n;
    void * m_integrationParameters;
};
}
#endif // INTEGRAL_H
