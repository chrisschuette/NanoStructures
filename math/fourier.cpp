#include "fourier.h"

#include <fftw3.h>

namespace math {
namespace fourier {

void fourierTransform(vector<MKL_Complex16>& vec) {
    int N = vec.size();

    fftw_plan plan;
    fftw_complex *In;
    fftw_complex *Out;
    In = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * N);
    Out = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * N);
    plan = fftw_plan_dft_1d(N, In, Out, FFTW_FORWARD, FFTW_ESTIMATE);
    //move data into memory
    for (unsigned int i = 0; i < vec.size(); i++) {
        In[i][0] = vec.at(i).real;
        In[i][1] = vec.at(i).imag;
    }

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (unsigned int i = 0; i < vec.size(); i++) {
        vec[i].real = Out[i][0];
        vec[i].imag = Out[i][1];
    }

    fftw_free(In);
    fftw_free(Out);
}

void inverseFourierTransform(vector<MKL_Complex16>& vec) {
    int N = vec.size();

    fftw_plan plan;
    fftw_complex *In;
    fftw_complex *Out;
    In = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * N);
    Out = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * N);
    plan = fftw_plan_dft_1d(N, In, Out, FFTW_BACKWARD, FFTW_ESTIMATE);
    //move data into memory
    for (unsigned int i = 0; i < vec.size(); i++) {
        In[i][0] = vec.at(i).real;
        In[i][1] = vec.at(i).imag;
    }

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (unsigned int i = 0; i < vec.size(); i++) {
        vec[i].real = Out[i][0] / (double) N;
        vec[i].imag = Out[i][1] / (double) N;
    }

    fftw_free(In);
    fftw_free(Out);
}

void fourierTransform(vector<std::complex<double> >& vec) {
    int N = vec.size();

    fftw_plan plan;
    fftw_complex *In;
    fftw_complex *Out;
    In = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * N);
    Out = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * N);
    plan = fftw_plan_dft_1d(N, In, Out, FFTW_FORWARD, FFTW_ESTIMATE);
    //move data into memory
    for (unsigned int i = 0; i < vec.size(); i++) {
        In[i][0] = vec.at(i).real();
        In[i][1] = vec.at(i).imag();
    }

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (unsigned int i = 0; i < vec.size(); i++) {
        vec[i].real() = Out[i][0];
        vec[i].imag() = Out[i][1];
    }

    fftw_free(In);
    fftw_free(Out);
}

void inverseFourierTransform(vector< std::complex<double> >& vec) {
    int N = vec.size();

    fftw_plan plan;
    fftw_complex *In;
    fftw_complex *Out;
    In = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * N);
    Out = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * N);
    plan = fftw_plan_dft_1d(N, In, Out, FFTW_BACKWARD, FFTW_ESTIMATE);
    //move data into memory
    for (unsigned int i = 0; i < vec.size(); i++) {
        In[i][0] = vec.at(i).real();
        In[i][1] = vec.at(i).imag();
    }

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for (unsigned int i = 0; i < vec.size(); i++) {
        vec[i].real() = Out[i][0] / (double) N;
        vec[i].imag() = Out[i][1] / (double) N;
    }

    fftw_free(In);
    fftw_free(Out);
}
}
}
