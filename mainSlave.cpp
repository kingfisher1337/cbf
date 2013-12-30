#include <complex>
#include "feynman.hpp"
#include <iostream>
#include <mpi/mpi.h>
#include "tensor.hpp"
using namespace std;

int mainSlave(
    int rank,
    int interpolationLevel,
    int Nk, int Nkomega) {
    
    int M;
    Tensor1<double> densities;
    Tensor1<double> masses;
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    densities.resize(M);
    masses.resize(M);
    MPI_Bcast(densities.getMemory(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(masses.getMemory(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    int kmax;
    Tensor1<double> kValues;
    MPI_Bcast(&kmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    kmax = (kmax-1)/2;
    kValues.resize(-kmax, +kmax);
    MPI_Bcast(kValues.getMemory(), kValues.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double deltak = (kValues(kValues.baseIndex(0)+kValues.size()+1) - kValues(kValues.baseIndex(0))) / (kValues.size() - 1);
    
    {
        int ctrl[3];
        MPI_Status status;
        Tensor4< complex<double> > Spart(0, interpolationLevel+2, -kmax, +kmax, 0, M-1, 0, M-1);
        FeynmanResult feynmanResult;
        lapackHermitianEigensystemHandle * h = lapackHermitianEigensystemInit(M);
        
        do {
            MPI_Recv(ctrl, 3, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (ctrl[0] != 0) {
                int kx1 = ctrl[1];
                int kx2 = ctrl[2];
                MPI_Recv(Spart.getMemory(), (kx2-kx1+1)*(2*kmax+1)*M*M, MPI_DOUBLE_COMPLEX, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                Spart.setBaseIndex(0, kx1);
                solveFeynmanProblem(Spart, masses, densities, deltak, feynmanResult, h);
                MPI_Send(&ctrl[1], 2, MPI_INT, 0, 1, MPI_COMM_WORLD);
                int n = (kx2-kx1+1) * (2*kmax+1) * M;
                MPI_Send(feynmanResult.omega.getMemory(), n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                n *= M;
                MPI_Send(feynmanResult.psi.getMemory(), n, MPI_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
                MPI_Send(feynmanResult.phi.getMemory(), n, MPI_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
                MPI_Send(feynmanResult.zeta.getMemory(), n, MPI_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
        } while (ctrl[0] != 0);
        
        lapackHermitianEigensystemFree(h);
    }
    
    Tensor3<double> omega(-kmax, +kmax, -kmax, +kmax, 0, M-1);
    Tensor4< complex<double> > psi(-kmax, +kmax, -kmax, +kmax, 0, M-1, 0, M-1);
    Tensor4< complex<double> > phi(-kmax, +kmax, -kmax, +kmax, 0, M-1, 0, M-1);
    Tensor4< complex<double> > zeta(-kmax, +kmax, -kmax, +kmax, 0, M-1, 0, M-1);
    MPI_Bcast(omega.getMemory(), omega.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(psi.getMemory(), psi.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(phi.getMemory(), phi.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(zeta.getMemory(), zeta.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    
    return 0;
}
