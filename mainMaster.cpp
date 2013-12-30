#include <complex>
#include "io.hpp"
#include <mpi/mpi.h>
#include "tensor.hpp"
#include <queue>
using namespace std;

static queue<int> workers;
static int M;
static int kmax;
static Tensor3<double> omega;
static Tensor4< complex<double> > psi;
static Tensor4< complex<double> > phi;
static Tensor4< complex<double> > zeta;
static int jobsStarted;
static int jobsDone;

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
    int Nk, int Nkint, int Nkomega) {
    
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    for (int i = 1; i < worldSize; ++i) {
        workers.push(i);
        cout << "[MASTER]: adding slave#" << i << " to queue" << endl;
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
    
    MPI_Bcast(omega.getMemory(), omega.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(psi.getMemory(), psi.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(phi.getMemory(), phi.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    MPI_Bcast(zeta.getMemory(), zeta.size(), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
    
    return 0;
}
