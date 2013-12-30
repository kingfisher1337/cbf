/* 
 * File:   io.hpp
 * Author: kingfisher1337
 *
 * Created on December 27, 2013, 7:34 PM
 */

#ifndef IO_HPP
#define	IO_HPP

#include <complex>
#include <string>
#include "tensor.hpp"

void readInfoFileMartin(
    const std::string & filename, 
    int & gridSize,
    int & numLayers,
    Tensor1<double> & densities, 
    Tensor1<double> & masses);

void readGroundstateMartin(
    const std::string & filename,
    int interpolationLevel,
    int gridSize,
    int kmax,
    int numLayers,
    Tensor1<double> & kValues,
    Tensor4< std::complex<double> > & staticStructureFactor,
    void (*kValuesCallback)(TensorBase1<double> &),
    void (*progressCallback)(TensorBase4< std::complex<double> > &, int, int));

#endif	/* IO_HPP */
