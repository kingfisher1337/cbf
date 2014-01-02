/* 
 * File:   cbf.hpp
 * Author: kingfisher1337
 *
 * Created on December 29, 2013, 12:06 PM
 */

#ifndef CBF_HPP
#define	CBF_HPP

#include <complex>
#include "tensor.hpp"

void solveCBFSelfEnergy(
    const TensorBase1<double> & rho,
    const TensorBase1<double> & masses,
    const TensorBase3<double> & omega,
    const TensorBase4< std::complex<double> > & psi,
    const TensorBase4< std::complex<double> > & phi,
    const TensorBase4< std::complex<double> > & zeta,
    double deltak, double domega, double epsilon,
    int kx, int ky, int m, int n,
    TensorBase1< std::complex<double> > & Sigma);

#endif	/* CBF_HPP */

