/* 
 * File:   main.cpp
 * Author: kingfisher1337
 *
 * Created on December 27, 2013, 4:45 PM
 */

#include <cstdlib>
#include <fstream>
#include "io.hpp"
#include <iostream>
//#include <mpi.h>
using namespace std;

int mainMaster(int argc, char ** argv) {
    return 0;
}

int mainSlave(int rank) {
    return 0;
}

int main(int argc, char ** argv) {
//    int rank;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    
//    int error = rank == 0 ? mainMaster(argc, argv) : mainSlave(rank);
//    
//    MPI_Finalize();
//    
//    return error;
    
    int gridSize;
    int numLayers;
    Tensor1<double> densities;
    Tensor1<double> masses;
    readInfoFileMartin(
        "/media/michael/OS/tmp/bac/dplayer/256density/2layer/ungekippt/0.06distance/daten.txt", 
        gridSize, numLayers, densities, masses);
    
    Tensor4< complex<double> > S;
    Tensor1<double> kValues;
    readGroundstateMartin(
        "/media/michael/OS/tmp/bac/dplayer/256density/2layer/ungekippt/0.06distance/s.dat",
        1, gridSize, 5, numLayers, kValues, S);
    
    int kmax = (S.size(0)-1)/2;
    ofstream o("/media/michael/OS/tmp/bac/dplayer/256density/2layer/ungekippt/0.06distance/Stest.dat");
    for (int kx = -kmax; kx <= kmax; ++kx) {
        for (int ky = -kmax; ky <= kmax; ++ky) {
            o << kValues(kx) << " " << kValues(ky);
            for (int a = 0; a < numLayers; ++a) {
                for (int b = 0; b < numLayers; ++b) {
                    o << " " << S(kx, ky, a, b).real()
                      << " " << S(kx, ky, a, b).imag();
                }
            }
            o << endl;
        }
        o << endl;
    }
    o.close();
    
    return 0;
}

