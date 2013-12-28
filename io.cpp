#include "io.hpp"
#include <fstream>
using namespace std;

static void skipDoubles(istream & in, int n) {
    double x;
    for (int i = 0; i < n; ++i) {
        in >> x;
    }
}
static void skipLine(istream & in, int M) {
    skipDoubles(in, 2*M*M+2);
}
static void skipLines(istream & in, int n, int M) {
    for (int i = 0; i < n; ++i) {
        skipLine(in, M);
    }
}

static void interpolateKVals(Tensor1<double> & kvals, int level) {
    int n = kvals.baseIndex(0) + kvals.size();
    for (int i = kvals.baseIndex(0) + level + 1; i < n; i += (level + 1)) {
        for (int j = 0; j < level; ++j) {
            kvals(i-level+j) = kvals(i-level-1) + (kvals(i) - kvals(i-level-1)) * (j+1) / (level+1);
        }
    }
}
static void interpolateSLine(TensorBase3< complex<double> > & S, int level) {
    const int kmax = (S.size(0)-1)/2;
    const int M = S.size(1);
    for (int k = -kmax+level+1; k <= kmax; k += (level+1)) {
        for (int k2 = 0; k2 < level; ++k2) {
            for (int a = 0; a < M; ++a) {
                for (int b = 0; b < M; ++b) {
                    S(k-level+k2, a, b) = S(k-level-1, a, b) + (double)(k2+1) / (level+1) * (S(k, a, b) - S(k-level-1, a, b));
                }
            }
        }
    }
}
static void interpolateSStripe(TensorBase4< complex<double> > & S, int kx, int level) {
    const int kmax = (S.size(0)-1)/2;
    const int M = S.size(2);
    for (int ky = -kmax; ky <= kmax; ++ky) {
        for (int kx2 = 0; kx2 < level; ++kx2) {
            for (int a = 0; a < M; ++a) {
                for (int b = 0; b < M; ++b) {
                    S(kx-level+kx2, ky, a, b) = S(kx-level-1, ky, a, b) + (double)(kx2+1) / (level+1) * (S(kx, ky, a, b) - S(kx-level-1, ky, a, b));
                }
            }
        }
    }
}

static void interpolateKVals(
    const TensorBase1<double> & kvals,
    Tensor1<double> & kvalsInterpolated,
    const size_t level,
    const size_t kMin,
    const size_t kMax) {
    
    const size_t N = kMax - kMin + 1;
    const size_t N2 = N + level * (N - 1);
    const double dk = 1.0 / (level+1);
    
    kvalsInterpolated.resize(N2);
    
    for (size_t k = kMin; k <= kMax; ++k) {
        kvalsInterpolated((k-kMin) * (level+1)) = kvals(k);
        if (k < kMax) {
            for (size_t ki = 1; ki <= level; ++ki) {
                kvalsInterpolated((k-kMin) * (level+1) + ki) = kvals(k) + (kvals(k+1) - kvals(k)) * ki * dk;
            }
        }
    }
}

void readInfoFileMartin(
    const string & filename, 
    int & gridSize,
    int & numLayers,
    Tensor1<double> & densities, 
    Tensor1<double> & masses) {

    string line;
    ifstream i(filename.c_str());
    int N, M;
    int j = 0;
    
    M = 0; // sufficient to loop until M is set
    while (j < 4 + 2 * M) {
        getline(i, line);
        if (line[0] != '#') {
            if (j == 1) {
                N = atoi(line.c_str());
                gridSize = N;
            } else if (j == 3) {
                M = atoi(line.c_str());
                numLayers = M;
                densities.resize(M);
                masses.resize(M);
            } else if (j >= 4 && j < 4 + M) {
                densities(j-4) = atoi(line.c_str());
            } else if (j >= 4 + M && j < 4 + 2 * M) {
                masses(j-4-M) = atoi(line.c_str());
            }
            ++j;
        }
    }
    
    i.close();
    
}

void readGroundstateMartin(
    const string & filename,
    int interpolationLevel,
    int gridSize,
    int kmax,
    int numLayers,
    Tensor1<double> & kValues,
    Tensor4< complex<double> > & staticStructureFactor) {
    // reads static structure factor from file,
    // where it is given at gridSize² grid points in momentum space;
    // but we are only interested in (2kmax+1)² grid points
    // (symmetric around the origin in momentum space);
    
    const int M = numLayers;
    
    ifstream in(filename.c_str());
    
    int k0 = 0;
    double x, y;
    skipDoubles(in, 1);
    in >> y;
    skipDoubles(in, 2*M*M+1);
    in >> x;
    while (fabs(x) < fabs(y)) {
        skipDoubles(in, 2*M*M+1);
        y = x;
        in >> x;
        ++k0;
    }
    
    int N; // output is given on a grid of (2N+1)² points
    if (interpolationLevel < 0) {
        N = kmax / (1 - interpolationLevel);
    } else {
        N = kmax * (interpolationLevel + 1);
    }
    staticStructureFactor.resize(-N, N, -N, N, 0, M-1, 0, M-1);
    kValues.resize(-N, N);
    
    in.seekg(0);
    
    int outputStride = interpolationLevel > 0 ? interpolationLevel+1 : 1;
    int kx = -N;
    int ky = -N;
    for (int i = 0; i < gridSize && kx <= N; ++i) {
        for (int j = 0; j < gridSize && kx <= N; ++j) {
            if (i-k0 < -kmax || i-k0 > kmax || j-k0 < -kmax || j-k0 > kmax) {
                skipLine(in, M);
            } else if (interpolationLevel < 0 && ((i-k0)%(1-interpolationLevel) != 0 || (j-k0)%(1-interpolationLevel) != 0)) {
                skipLine(in, M);
            } else {
                double kxVal, kyVal;
                in >> kxVal >> kyVal;
                if (kx == -N) {
                    kValues(ky) = kyVal;
                }
                for (int a = 0; a < M; ++a) {
                    for (int b = 0; b < M; ++b) {
                        double re, im;
                        in >> re >> im;
                        staticStructureFactor(kx, ky, a, b) = complex<double>(re, im);
                    }
                }
                ky += outputStride;
                if (ky > N) {
                    if (interpolationLevel > 0) {
                        if (kx == -N) {
                            interpolateKVals(kValues, interpolationLevel);
                        }
                        TensorBase3< complex<double> > Sline = staticStructureFactor(kx);
                        interpolateSLine(Sline, interpolationLevel);
                        if (kx > -N) {
                            interpolateSStripe(staticStructureFactor, kx, interpolationLevel);
                        }
                    }
                    ky = -N;
                    kx += outputStride;
                }
            }
        }
    }
    
    in.close();
    
}