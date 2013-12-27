/* 
 * File:   main.cpp
 * Author: kingfisher1337
 *
 * Created on December 27, 2013, 4:45 PM
 */

#include <cstdlib>
#include <mpi.h>
using namespace std;

int mainMaster(int argc, char ** argv) {
    return 0;
}

int mainSlave(int rank) {
    return 0;
}

int main(int argc, char ** argv) {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int error = rank == 0 ? mainMaster(argc, argv) : mainSlave(rank);
    
    MPI_Finalize();
    
    return error;
}

