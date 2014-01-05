#include "cbf.hpp"
#include <complex>
#include "feynman.hpp"
#include <iostream>
#include <mpi.h>
#include <fstream>
#include "tensor.hpp"
using namespace std;

int mainSlave(
    int rank,
    int interpolationLevel,
    int Nk, int Nomega,
    bool feynmansOnly,
    bool selfenergyOnly) {
    
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
                
//                {
//                    stringstream s;
//                    s << "/media/michael/OS/tmp/bac/excitations/256density/2layer/ungekippt/0.06distance/S/" << kx1 << ".dat";
//                    ofstream o(s.str().c_str());
//                    for (int kx = kx1; kx <= kx2; ++kx) {
//                        for (int ky = -kmax; ky <= kmax; ++ky) {
//                            o << kValues(kx) << " " << kValues(ky);
//                            for (int a = 0; a < M; ++a) {
//                                for (int b = 0; b < M; ++b) {
//                                    o << " " << real(Spart(kx, ky, a, b))
//                                      << " " << imag(Spart(kx, ky, a, b));
//                                }
//                            }
//                            o << endl;
//                        }
//                        o << endl;
//                    }
//                    o.close();
//                }
                
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
    
    if (feynmansOnly) {
        return 0;
    }
    
    double domega;
    double epsilon;
    Tensor3<double> omega(-kmax, +kmax, -kmax, +kmax, 0, M-1);
    Tensor4< complex<double> > psi(-kmax, +kmax, -kmax, +kmax, 0, M-1, 0, M-1);
    Tensor4< complex<double> > phi(-kmax, +kmax, -kmax, +kmax, 0, M-1, 0, M-1);
    Tensor4< complex<double> > zeta(-kmax, +kmax, -kmax, +kmax, 0, M-1, 0, M-1);
    MPI_Bcast(&domega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(omega.getMemory(), omega.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(psi.getMemory(), psi.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(phi.getMemory(), phi.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(zeta.getMemory(), zeta.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    
    // domega
    
    {
        int ctrl[5];
        MPI_Status status;
        Tensor1< complex<double> > Sigma(Nomega);
        do {
            MPI_Recv(ctrl, 5, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            
            if (ctrl[0] != 0) {                
                solveCBFSelfEnergy(
                    densities, masses, omega, psi, phi, zeta, 
                    deltak, domega, epsilon, 
                    ctrl[1], ctrl[2], ctrl[3], ctrl[4], 
                    Sigma);
                
                MPI_Send(&ctrl[1], 4, MPI_INT, 0, 1, MPI_COMM_WORLD);
                MPI_Send(Sigma.getMemory(), Sigma.size(), MPI_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
            }
        } while (ctrl[0] != 0);
    }
    
    return 0;
}
