#include "lapack.hpp"
#include <complex>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include "tensor.hpp"
using namespace std;

extern "C" {
    void zgetrf_(int* M, int *N, complex<double>* A, int* lda, int* IPIV, int* INFO);
    void zgetri_(int* N, complex<double>* A, int* lda, int* IPIV, complex<double>* WORK, int* lwork, int* INFO);
}

struct lapackInvertHandle {
    int n;
    int * ipiv;
    int lwork;
    complex<double> * work;
};
lapackInvertHandle * lapackInvertInit(int n) {
    lapackInvertHandle * h = new lapackInvertHandle;
    h->n = n;
    h->ipiv = new int[n+1];
    h->lwork = n*n;
    h->work = new complex<double>[n*n];
    return h;
}
void lapackInvertFree(lapackInvertHandle * h) {
    if (h != 0) {
        if (h->ipiv != 0) {
            delete[] h->ipiv;
        }
        if (h->work != 0) {
            delete[] h->work;
        }
        delete h;
    }
}
void lapackInvert(lapackInvertHandle * h, TensorBase2< std::complex<double> > & A) {
    int INFO;
    complex<double> * p = A.getMemory();
    zgetrf_(&h->n, &h->n, p, &h->n, h->ipiv, &INFO);
    zgetri_(&h->n, p, &h->n, h->ipiv, h->work, &h->lwork, &INFO);
}


extern "C" {
    void zheev_(char* jobz, char* uplo, int* n, complex<double>* a, int* lda, double* w, complex<double>* work, int* lwork, double* rwork, int* info);
}

struct lapackHermitianEigensystemHandle {
    int n;
    complex<double> * backupMem;
    double * rwork;
    complex<double> * work;
    int lwork;
    lapackHermitianEigensystemResult res;
};
lapackHermitianEigensystemHandle * lapackHermitianEigensystemInit(int n) {
    if (n < 1) {
        return 0;
    }
    lapackHermitianEigensystemHandle * h = new lapackHermitianEigensystemHandle;
    h->n = n;
    h->backupMem = new complex<double>[n * n];
    h->rwork = new double[3*n-2];
    h->work = 0;
    h->lwork = 0;
    h->res.eigenvalues.resize(n);
    h->res.eigenvectors.resize(n, n);
    return h;
}
void lapackHermitianEigensystemFree(lapackHermitianEigensystemHandle * h) {
    if (h != 0) {
        if (h->backupMem != 0) {
            free(h->backupMem);
        }
        if (h->rwork != 0) {
            free(h->rwork);
        }
        if (h->work != 0) {
            free(h->work);
        }
        delete h;
    }
}
lapackHermitianEigensystemResult * lapackHermitianEigensystem(
    lapackHermitianEigensystemHandle * h,
    TensorBase2< complex<double> > & A) {
    
    int info;
    complex<double> * p = A.getMemory();
    complex<double> * pevec = h->res.eigenvectors.getMemory();
    double * peval = h->res.eigenvalues.getMemory();
    memcpy(h->backupMem, p, sizeof(complex<double>)*h->n*h->n);
    int lworktmp = -1;
    complex<double> opt;
    char job = 'V';
    char uplo = 'L';
    zheev_(&job, &uplo, &h->n, p, &h->n, peval, &opt, &lworktmp, h->rwork, &info);
    lworktmp = (int)opt.real();
    if (h->work == 0) {
        h->work = new complex<double>[lworktmp];
        h->lwork = lworktmp;
    } else if (h->lwork < lworktmp) {
        delete[] h->work;
        h->work = new complex<double>[lworktmp];
        h->lwork = lworktmp;
    }
    zheev_(&job, &uplo, &h->n, p, &h->n, peval, h->work, &lworktmp, h->rwork, &info);
    memcpy(pevec, p, sizeof(complex<double>)*h->n*h->n);
    memcpy(p, h->backupMem, sizeof(complex<double>)*h->n*h->n);
    
    return &h->res;
}
