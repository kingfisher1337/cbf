/* 
 * File:   lapack.hpp
 * Author: kingfisher1337
 *
 * Created on December 29, 2013, 3:51 PM
 */

#ifndef LAPACK_HPP
#define	LAPACK_HPP

#include <complex>
#include "tensor.hpp"

struct lapackInvertHandle;
lapackInvertHandle * lapackInvertInit(int n);
void lapackInvertFree(lapackInvertHandle * h);
void lapackInvert(lapackInvertHandle * h, TensorBase2< std::complex<double> > & A);

struct lapackHermitianEigensystemHandle;
struct lapackHermitianEigensystemResult {
    Tensor1<double> eigenvalues;
    Tensor2< std::complex<double> > eigenvectors;
};
lapackHermitianEigensystemHandle * lapackHermitianEigensystemInit(int n);
void lapackHermitianEigensystemFree(lapackHermitianEigensystemHandle * h);
lapackHermitianEigensystemResult * lapackHermitianEigensystem(
    lapackHermitianEigensystemHandle * h,
    TensorBase2< std::complex<double> > & A);

#endif	/* LAPACK_HPP */

