#include <complex>
#include "tensor.hpp"
using namespace std;

#define pi 3.1415926535897932384626433832795028841971693993751058

#define kv(k) k##x, k##y
#define kdot(k1, k2) ((k1##x * k2##x + k1##y * k2##y) * deltak * deltak)
#define ksquare(k) dot(k, k)


//void solveCBF(
//    const TensorBase3<double> & omega,
//    const TensorBase4< complex<double> > & psi,
//    const TensorBase4< complex<double> > & phi,
//    const TensorBase4< complex<double> > & zeta,
//    int Nk, int Nomega, double domega,
//    int kx, int ky, int a, int b,
//    Tensor3< complex<double> > & chi) {
//    
//    const int kmax = (omega.size(0)-1)/2;
//    const int M = omega.size(1);
//    
//    if (chi.size() == 0) {
//        chi.resize(Nomega, M, M);
//    }
//    
//}

static complex<double> Vstn(
    const TensorBase1<double> & rho,
    const TensorBase1<double> & masses,
    const TensorBase3<double> & omega,
    const TensorBase4< complex<double> > & psi,
    const TensorBase4< complex<double> > & phi,
    const TensorBase4< complex<double> > & zeta,
    int ksx, int ksy, 
    int ktx, int kty, 
    int knx, int kny,
    int s, int t, int n,
    double deltak) {
    
    complex<double> result(0, 0);
    
    const int M = omega.size(1);
    
    for (int a = 0; a < M; ++a) {
        result -= (
            conj(psi(kv(kn), n, a)) * (
                phi(kv(ks), s, a) * zeta(kv(kt), t, a) * kdot(kt, kn) +
                phi(kv(kt), t, a) * zeta(kv(ks), s, a) * kdot(ks, kn)
            ) - sqrt(rho(a)) * omega(kv(kn), n) / (omega(kv(kn), n) + omega(kv(ks), s) + omega(kv(kt), t)) * (
                conj(phi(kv(kn), n, a)) * zeta(kv(ks), s, a) * zeta(kv(kt), t, a) * kdot(ks, kt) -
                phi(kv(ks), s, a) * conj(zeta(kv(kn), n, a)) * zeta(kv(kt), t, a) * kdot(kn, kt) -
                phi(kv(kt), t, a) * conj(zeta(kv(kn), n, a)) * zeta(kv(ks), s, a) * kdot(kn, ks)
            )
        ) / masses(a);
    }
    
    return result / 2.0;
}

void solveCBFSelfEnergy(
    const TensorBase1<double> & rho,
    const TensorBase1<double> & masses,
    const TensorBase3<double> & omega,
    const TensorBase4< complex<double> > & psi,
    const TensorBase4< complex<double> > & phi,
    const TensorBase4< complex<double> > & zeta,
    double deltak, int Nomega, double domega, double epsilon,
    int kx, int ky, int m, int n,
    TensorBase1< complex<double> > & Sigma) {
    
    const int M = omega.size(1);
    const int kmax = (omega.size(0)-1)/2;
    const complex<double> iepsilon(0, epsilon);
    
    const int ksx1 = kx - kmax;
    const int ksx2 = kx + kmax;
    const int ksy1 = ky - kmax;
    const int ksy2 = ky + kmax;
    
    double f = deltak / (2*pi);
    f = f*f / 2;
    //memset(Sigma.getMemory(), 0, NOMEGA * sizeof(complex<double>));
    for (int s = 0; s < M; ++s) {
        for (int t = 0; t < M; ++t) {
            for (int ksx = ksx1; ksx <= ksx2; ++ksx) {
                for (int ksy = ksy1; ksy <= ksy2; ++ksy) {
//            for (int ksx = -kmax; ksx <= kmax; ++ksx) {
//                for (int ksy = -kmax; ksy <= kmax; ++ksy) {
                    int ktx = kx - ksx;
                    int kty = ky - ksy;
//                    if (ktx < -kmax || ktx > kmax || kty < -kmax || kty > kmax) {
//                        continue; // todo: do it more efficient with a break
//                    }
                    
                    complex<double> x = Vstn(rho, masses, omega, phi, psi, zeta, kv(ks), kv(kt), kv(k), s, t, m, deltak) *
                                        conj(Vstn(rho, masses, omega, phi, psi, zeta, kv(ks), kv(kt), kv(k), s, t, n, deltak)) * f;
                    for (int iomega = 0; iomega < Nomega; ++iomega) {
                        Sigma(iomega) += x / (iomega*domega - omega(kv(ks), s) - omega(kv(kt), t) - iepsilon);
                    }
                }
            }
        }
    }
}
