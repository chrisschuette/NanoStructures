#include "hilberttransformer.h"
#include "fourier.h"

namespace math {

static MKL_Complex16 complexZero = {0.0, 0.0};
int HilbertTransformer::halfKernelSize = 32000;
int HilbertTransformer::kernelSize = 2 * HilbertTransformer::halfKernelSize + 1;

HilbertTransformer::HilbertTransformer()
{
}
void HilbertTransformer::init() {
    polesSize = posPoles.getSize() + negPoles.getSize();
    numPointsFT = (int) std::pow(2.0,
          int(log((polesSize+kernelSize-1.0))/log(2.0)+0.5));
    setupKernel();
}

void HilbertTransformer::setupKernel() {
    double logr = log(posPoles.getArgument(1)/posPoles.getArgument(0));/// fix to be general
   // OUT(3) << "posPoles.grid[1] = " << posPoles.getArgument(1) << EOM;
  //  OUT(3) << "posPoles.grid[0] = " << posPoles.getArgument(0) << EOM;

   // OUT(3) << "logr = " << logr << EOM;
    kernelB.clear();
    kernelB.resize(numPointsFT, complexZero);
    kernelF.clear();
    kernelF.resize(numPointsFT, complexZero);

    for(int i = 0; i < kernelSize; ++i)
       kernelF[i].real = logr/(exp((i-halfKernelSize)*logr)+1.0);

    for(int i = 0; i < halfKernelSize; ++i)
       kernelB[i].real = logr/(exp((i-halfKernelSize)*logr)-1.0);
    for(int i = halfKernelSize + 1; i < kernelSize; ++i)
       kernelB[i].real = logr/(exp((i-halfKernelSize)*logr)-1.0);
    kernelF[halfKernelSize].real = 0.0;

    fourier::fourierTransform(kernelB);
    fourier::fourierTransform(kernelF);
}

void HilbertTransformer::setupWeights() {
    weightsPos.clear();
    weightsPos.resize(numPointsFT, complexZero);
    weightsNeg.clear();
    weightsNeg.resize(numPointsFT, complexZero);

    for(int i = 0; i < posPoles.getSize(); ++i)
       weightsPos[i].real = posPoles.getValue(i);

    for(int i = 0; i < negPoles.getSize(); ++i)
       weightsNeg[i].real = negPoles.getValue(i);

    fourier::fourierTransform(weightsPos);
    fourier::fourierTransform(weightsNeg);
}

void HilbertTransformer::doConvolution() {
    posHT.clear();
    posHT.resize(posPoles.getSize());
    negHT.clear();
    negHT.resize(negPoles.getSize());

    broadenedB.clear();
    broadenedB.resize(numPointsFT);
    broadenedF.clear();
    broadenedF.resize(numPointsFT);

    for(int i = 0; i < numPointsFT; ++i)
    {
       broadenedB[i].real = weightsPos[i].real * kernelB[i].real
          - weightsPos[i].imag * kernelB[i].imag;
       broadenedB[i].imag = weightsPos[i].real * kernelB[i].imag
          + weightsPos[i].imag * kernelB[i].real;
       broadenedF[i].real = weightsNeg[i].real * kernelF[i].real
          - weightsNeg[i].imag * kernelF[i].imag;
       broadenedF[i].imag = weightsNeg[i].real * kernelF[i].imag
          + weightsNeg[i].imag * kernelF[i].real;
    }

    fourier::inverseFourierTransform(broadenedB);
    fourier::inverseFourierTransform(broadenedF);

    for(int i = 0; i < posPoles.getSize(); ++i)
    {
       posHT[i] = broadenedB[i+halfKernelSize].real
          + broadenedF[i+halfKernelSize].real;
    }

    broadenedB.clear();
    broadenedB.resize(numPointsFT);
    broadenedF.clear();
    broadenedF.resize(numPointsFT);

    for(int i = 0; i < numPointsFT; ++i)
    {
       broadenedB[i].real = weightsNeg[i].real * kernelB[i].real
          - weightsNeg[i].imag * kernelB[i].imag;
       broadenedB[i].imag = weightsNeg[i].real * kernelB[i].imag
          + weightsNeg[i].imag * kernelB[i].real;
       broadenedF[i].real = weightsPos[i].real * kernelF[i].real
          - weightsPos[i].imag * kernelF[i].imag;
       broadenedF[i].imag = weightsPos[i].real * kernelF[i].imag
          + weightsPos[i].imag * kernelF[i].real;
    }

    fourier::inverseFourierTransform(broadenedB);
    fourier::inverseFourierTransform(broadenedF);

    for(int i = 0; i < negPoles.getSize(); ++i)
    {
       negHT[i] = -broadenedB[i+halfKernelSize].real
          - broadenedF[i+halfKernelSize].real;
    }
}

void HilbertTransformer::broaden() {
    setupWeights();
    doConvolution();
     for(unsigned int i = 0; i < posHT.size() ; i++)
         posHT[i] /= M_PI;
     for(unsigned int i = 0; i < negHT.size() ; i++)
         negHT[i] /= M_PI;
}
}
