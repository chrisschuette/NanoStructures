#ifndef HILBERTTRANSFORMER_H
#define HILBERTTRANSFORMER_H

#include "mkl_types.h"
#include "function.h"

#include <vector>

namespace math {
class HilbertTransformer
{
public:
    HilbertTransformer();
    void init();
    void setupKernel();
    void setupWeights();
    void doConvolution();
    void broaden();
    void setPoles(math::Function positive, math::Function negative) { posPoles = positive; negPoles = negative; }
    std::vector<double> posHT;
    std::vector<double> negHT;
protected:
    int polesSize;
    int polesPosSize;
    int polesNegSize;
    static int halfKernelSize;
    static int kernelSize;
    int numPointsFT;

    math::Function posPoles;
    math::Function negPoles;

    std::vector<MKL_Complex16> weightsPos;
    std::vector<MKL_Complex16> weightsNeg;
    std::vector<MKL_Complex16> kernelF;
    std::vector<MKL_Complex16> kernelB;
    std::vector<MKL_Complex16> broadenedF;
    std::vector<MKL_Complex16> broadenedB;

};
}
#endif // HILBERTTRANSFORMER_H
