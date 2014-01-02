#include <complex>
#include <fstream>
#include "io.hpp"
#include <mpi/mpi.h>
#include "tensor.hpp"
#include <queue>
#include <set>
using namespace std;

static queue<int> workers;
static int M;
static int Nk;
static int Nomega;
static int kmax;
static Tensor3<double> omega;
static Tensor4< complex<double> > psi;
static Tensor4< complex<double> > phi;
static Tensor4< complex<double> > zeta;
static Tensor1< complex<double> > Sigma;
static int jobsStarted;
static int jobsDone;
static ofstream ofsSelfEnergy;

static void waitForFeynmanWorker() {
    int ctrl[2];
    MPI_Status status;
    MPI_Recv(ctrl, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    int sender = status.MPI_SOURCE;
    
    int n = (ctrl[1]-ctrl[0]+1) * (2*kmax+1) * M;
    MPI_Recv(&omega(ctrl[0], -kmax, 0), n, MPI_DOUBLE, sender, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    n *= M;
    MPI_Recv(&psi(ctrl[0], -kmax, 0, 0), n, MPI_DOUBLE_COMPLEX, sender, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&phi(ctrl[0], -kmax, 0, 0), n, MPI_DOUBLE_COMPLEX, sender, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&zeta(ctrl[0], -kmax, 0, 0), n, MPI_DOUBLE_COMPLEX, sender, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    workers.push(sender);
    ++jobsDone;
}
static void waitForCBFSelfEnergyWorker() {
    int ctrl[4];
    MPI_Status status;
    MPI_Recv(ctrl, 4, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    int sender = status.MPI_SOURCE;
    MPI_Recv(Sigma.getMemory(), Nomega, MPI_DOUBLE_COMPLEX, sender, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    for (int iomega = 0; iomega < Nomega; ++iomega) {
        writeMatrixFileEntry(ofsSelfEnergy, Nk, Nomega, M, ctrl[0], iomega, ctrl[2]*M + ctrl[3], Sigma(iomega));
    }
    workers.push(sender);    
    ++jobsDone;
    cout << "[MASTER]: Sigma calculation for k = (" << ctrl[0] << ", " 
         << ctrl[1] << "), a = " << ctrl[2] << ", b = " << ctrl[3] << " done" 
         << endl;
}
static void kValuesCallback(TensorBase1<double> & kValues) {
    cout << "[MASTER]: k-values determined" << endl;
    int n = kValues.size();
    kmax = -kValues.baseIndex(0);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(kValues.getMemory(), kValues.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    cout << "[MASTER]: resizing feynman result storages for k € {" << -kmax << ", " << kmax << "}" << endl;
    omega.resize(-kmax, +kmax, -kmax, +kmax, 0, M-1);
    psi.resize(-kmax, +kmax, -kmax, +kmax, 0, M-1, 0, M-1);
    phi.resize(-kmax, +kmax, -kmax, +kmax, 0, M-1, 0, M-1);
    zeta.resize(-kmax, +kmax, -kmax, +kmax, 0, M-1, 0, M-1);
}
static void staticStructureFactorStripeCallback(TensorBase4< complex<double> > & S, int kx1, int kx2) {
    if (workers.size() == 0) {
        waitForFeynmanWorker();
    }
    int kmaxint = (S.size(0)-1)/2;
    int worker = workers.front();
    workers.pop();
    int ctrl[] = { 1, kx1, kx2 };
    MPI_Send(ctrl, 3, MPI_INT, worker, 1, MPI_COMM_WORLD);
    MPI_Send(&S(kx1, -kmaxint, 0, 0), S.size() / S.size(0) * (kx2-kx1+1), MPI_DOUBLE_COMPLEX, worker, 1, MPI_COMM_WORLD);
    ++jobsStarted;
}

int mainMaster(
    const string & infoFileName, const string & sFileName, 
    const string & opath,
    int interpolationLevel,
    int Nk_, int Nkint, int Nomega_,
    bool feynmansOnly,
    bool selfenergyOnly) {
    
    Nomega = Nomega_;
    Nk = Nk_;
    
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    for (int i = 1; i < worldSize; ++i) {
        workers.push(i);
    }
    
    cout << "[MASTER]: having " << worldSize << " workers" << endl;
    
    int gridSize;
    Tensor1<double> densities;
    Tensor1<double> masses;
    readInfoFileMartin(infoFileName, gridSize, M, densities, masses);
    
    cout << "[MASTER]: info file read; gridSize=" << gridSize << "; M=" << M << endl;
    
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(densities.getMemory(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(masses.getMemory(), M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    cout << "[MASTER]: bcasted basic info to workers" << endl;
    
    Tensor1<double> kValues;
    {
        Tensor4< complex<double> > S;
        jobsStarted = 0;
        jobsDone = 0;
        readGroundstateMartin(
            sFileName,
            interpolationLevel, gridSize, Nkint, M, kValues, S,
            kValuesCallback, staticStructureFactorStripeCallback);
    }
    
    cout << "[MASTER]: read groundstate file" << endl;
        
    while (jobsStarted > jobsDone) {
        waitForFeynmanWorker();
    }
    
    cout << "[MASTER]: all " << jobsDone << " jobs have finished" << endl;
    
    for (int i = 1; i < worldSize; ++i) {
        int ctrl = 0;
        MPI_Send(&ctrl, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
    }
    
    // todo: write feynmans to file
    
    double maxOmega = 0;
    for (int a = 0; a < M; ++a) {
        if (densities(a) > maxOmega) {
            maxOmega = densities(a);
        }
    }
    maxOmega *= 80; // rho_max --> omega_max
    double domega = maxOmega / (Nomega - 1);
    double epsilon = 3 * domega;
    
    MPI_Bcast(&domega, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(omega.getMemory(), omega.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(psi.getMemory(), psi.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(phi.getMemory(), phi.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(zeta.getMemory(), zeta.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    
    Sigma.resize(Nomega);
    ofsSelfEnergy.open((opath + "Sigma.dat").c_str());
    prepareMatrixFile(ofsSelfEnergy, Nk, kValues, Nomega, domega, M);
//    writeMatrixFileEntry(ofsSelfEnergy, Nk, Nomega, M, 0, 0, 0, 42.0);
//    writeMatrixFileEntry(ofsSelfEnergy, Nk, Nomega, M, 0, 1, 2, 43.0);
//    writeMatrixFileEntry(ofsSelfEnergy, Nk, Nomega, M, 1, 1, 1, 43.0);
    cout << "[MASTER]: preparing Sigma matrix file finished" << endl;
    
    while (!workers.empty()) {
        workers.pop();
    }
    for (int i = 1; i < worldSize; ++i) {
        workers.push(i);
    }
    
    if (interpolationLevel > 0) {
        Nk *= (interpolationLevel + 1);
    } else if (interpolationLevel < 0) {
        Nk /= (1 - interpolationLevel);
    }
    
    set<int> kxDone;
    set<int> kxCurrent;
    kxCurrent.insert(1);
    kxCurrent.insert(Nk);
    while (kxDone.size() != kmax) {
        for (set<int>::iterator it = kxCurrent.begin(); it != kxCurrent.end(); ++it) {
            int kx = *it;
            for (int a = 0; a < M; ++a) {
                for (int b = 0; b < M; ++b) {
                    if (workers.size() == 0) {
                        waitForCBFSelfEnergyWorker();
                    }
                    int worker = workers.front();
                    workers.pop();
                    int ctrl[] = { 1, kx, 0, a, b };
                    MPI_Send(ctrl, 5, MPI_INT, worker, 1, MPI_COMM_WORLD);
                    ++jobsStarted;
                }
            }
        }
        
        kxDone.insert(kxCurrent.begin(), kxCurrent.end());
        kxCurrent.clear();
        int prev = -1;
        for (set<int>::iterator it = kxDone.begin(); it != kxDone.end(); ++it) {
            if (prev != -1) {
                if (*it != prev + 1) {
                    kxCurrent.insert((*it + prev) / 2);
                }
            }
            prev = *it;
        }
    }
    
    while (jobsStarted > jobsDone) {
        waitForCBFSelfEnergyWorker();
    }
    
    for (int i = 1; i < worldSize; ++i) {
        int ctrl = 0;
        MPI_Send(&ctrl, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
    }
    
    ofsSelfEnergy.close();
    
    return 0;
}
