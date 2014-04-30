#ifndef LAPACK_H
#define LAPACK_H

namespace lapack {
extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a,
               int *lda, double *w, double *work, int *lwork,
               int *info);
extern "C" void dgemm_(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *a,int *lda,double *b,int *ldb,double *beta,double *c,int *ldc);

}

#endif // LAPACK_H
