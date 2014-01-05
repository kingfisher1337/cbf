/* 
 * File:   main.cpp
 * Author: kingfisher1337
 *
 * Created on December 27, 2013, 4:45 PM
 */

#include <mpi.h>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>

#include "Logger.hpp"
using namespace std;

int mainMaster(
    const string & infoFileName, const string & sFileName, 
    const string & opath,
    int interpolationLevel,
    int Nk_, int Nkint, int Nomega_,
    bool feynmansOnly,
    bool selfenergyOnly,
    Logger * logger);
int mainSlave(
    int rank,
    int interpolationLevel,
    int Nk, int Nomega,
    bool feynmansOnly,
    bool selfenergyOnly);

static bool startsWith(const char * s, const char * p) {
    int i;
    for (i = 0; s[i] != 0 && p[i] != 0; ++i) {
        if (s[i] != p[i]) {
            return false;
        }
    }
    return p[i] == 0;
}

int main(int argc, char ** argv) {
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    string ipath = ""; // optional
    string infofile = "-";
    string sfile = "-";
    string opath = ""; // optional
    int Nk = -1;
    int Nkint = -1;
    int Nomega = -1;
    int interpolationLevel = 0; // optional
    for (int i = 0; i < argc; ++i) {
        if (startsWith(argv[i], "ipath=")) {
            ipath = string(argv[i]).substr(6);
            if (ipath[ipath.length()-1] != '/') {
                ipath += "/";
            }
        } else if (startsWith(argv[i], "opath=")) {
            opath = string(argv[i]).substr(6);
            if (opath[opath.length()-1] != '/') {
                opath += "/";
            }
        } else if (startsWith(argv[i], "infofile=")) {
            infofile = string(argv[i]).substr(9);
        } else if (startsWith(argv[i], "sfile=")) {
            sfile = string(argv[i]).substr(6);
        } else if (startsWith(argv[i], "Nk=")) {
            Nk = atoi(string(argv[i]).substr(3).c_str());
        } else if (startsWith(argv[i], "Nkint=")) {
            Nkint = atoi(string(argv[i]).substr(6).c_str());
        } else if (startsWith(argv[i], "Nomega=")) {
            Nomega = atoi(string(argv[i]).substr(7).c_str());
        } else if (startsWith(argv[i], "interpolationLevel=")) {
            interpolationLevel = atoi(string(argv[i]).substr(19).c_str());
        }
    }
    
    int error = -1;
    if (infofile == "-") {
        if (rank == 0) {
            cerr << "no infofile was given!" << endl;
        }
    } else if (sfile == "-") {
        if (rank == 0) {
            cerr << "no sfile was given!" << endl;
        }
    } else if (Nk == -1) {
        if (rank == 0) {
            cerr << "Nk was not given!" << endl;
        }
    } else if (Nkint == -1) {
        if (rank == 0) {
            cerr << "Nkint was not given!" << endl;
        }
    } else if (Nomega == -1) {
        if (rank == 0) {
            cerr << "Nomega was not given!" << endl;
        }
    } else {
        
        Logger * logger = 0;
        ofstream ofsLogger;
        if (rank == 0) {
            stringstream filename;
            filename << opath << rank << ".log";
            ofsLogger.open(filename.str().c_str());
            logger = new Logger(ofsLogger, rank);
        }
        
        error = rank == 0 ?
            mainMaster(ipath + infofile, ipath + sfile, opath, interpolationLevel, Nk, Nkint, Nomega, false, true, logger) :
            mainSlave(rank, interpolationLevel, Nk, Nomega, false, true);
        
        ofsLogger.close();
        delete logger;
    }
    
    MPI_Finalize();
    
    return error;
}
