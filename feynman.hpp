/* 
 * File:   feynman.hpp
 * Author: kingfisher1337
 *
 * Created on December 29, 2013, 12:06 PM
 */

#ifndef FEYNMAN_HPP
#define	FEYNMAN_HPP

#include <complex>
#include "lapack.hpp"
#include "tensor.hpp"

struct FeynmanResult {
    Tensor3<double> omega;
    Tensor4< std::complex<double> > psi;
    Tensor4< std::complex<double> > phi;
    Tensor4< std::complex<double> > zeta;
};

void solveFeynmanProblem(
    Tensor4< std::complex<double> > & S,
    const Tensor1<double> & m,
    const Tensor1<double> & rho,
    double deltak,
    FeynmanResult & result,
    lapackHermitianEigensystemHandle * h);

#endif	/* FEYNMAN_HPP */

