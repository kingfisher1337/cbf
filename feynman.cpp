#include "feynman.hpp"
#include <iostream>
using namespace std;

void solveFeynmanProblem(
    Tensor4< complex<double> > & S,
    const Tensor1<double> & m,
    const Tensor1<double> & rho,
    double deltak,
    FeynmanResult & result,
    lapackHermitianEigensystemHandle * h) {
    
    const int kx1 = S.baseIndex(0);
    const int kx2 = kx1 + S.size(0) - 1;
    const int kmax = (S.size(1) - 1) / 2;
    const int M = S.size(2);
    
    if (result.omega.size(0) < (kx2 - kx1 + 1)) {
        result.omega.resize(0, kx2 - kx1, -kmax, +kmax, 0, M-1);
        result.psi.resize(0, kx2 - kx1, -kmax, +kmax, 0, M-1, 0, M-1);
        result.phi.resize(0, kx2 - kx1, -kmax, +kmax, 0, M-1, 0, M-1);
        result.zeta.resize(0, kx2 - kx1, -kmax, +kmax, 0, M-1, 0, M-1);
    }
    
    result.omega.setBaseIndex(0, kx1);
    result.psi.setBaseIndex(0, kx1);
    result.phi.setBaseIndex(0, kx1);
    result.zeta.setBaseIndex(0, kx1);
    
    for (int kx = kx1; kx <= kx2; ++kx) {
        for (int ky = -kmax; ky <= +kmax; ++ky) {
            for (int a = 0; a < M; ++a) {
                for (int b = 0; b < M; ++b) {
                    S(kx, ky, a, b) *= sqrt(m(a) * m(b));
                }
            }
            
            TensorBase2< complex<double> > tmp = S(kx, ky);            
            lapackHermitianEigensystemResult * r = lapackHermitianEigensystem(h, tmp);
                        
            for (int n = 0; n < M; ++n) {
                double k2 = (kx*kx + ky*ky) * deltak*deltak;
                result.omega(kx, ky, n) = k2 / 2 / r->eigenvalues(n);
                for (int a = 0; a < M; ++a) {
                    double f = (kx == 0 && ky == 0) ? 
                        10e10 : 
                        sqrt(2.0 * m(a) * result.omega(kx, ky, n) / k2);
                    result.psi(kx, ky, n, a) = r->eigenvectors(n, a) * f;
                    result.phi(kx, ky, n, a) = r->eigenvectors(n, a) / f;
                    result.zeta(kx, ky, n, a) = (result.phi(kx, ky, n, a) - result.psi(kx, ky, n, a)) / sqrt(rho(a));
                }
            }
        }
    }
}
