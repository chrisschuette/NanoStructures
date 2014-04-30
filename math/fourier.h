#ifndef FOURIER_H
#define FOURIER_H

#include <vector>
#include <complex>
#include "mkl_types.h"

using namespace std;
namespace math {
    namespace fourier {
        void fourierTransform(vector<std::complex<double> >& vec);
        void inverseFourierTransform(vector< std::complex<double> >& vec);
        void fourierTransform(vector<MKL_Complex16>& vec);
        void inverseFourierTransform(vector<MKL_Complex16>& vec);
    }
}
#endif // FOURIER_H
