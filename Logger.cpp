#include <ctime>
#include <iomanip>
#include "Logger.hpp"
#include <sstream>
#include <vector>
using namespace std;

static void printLineBeginning(ostream & o, const string & n) {
    static char tbuf[64];
    time_t t = time(0);
    tm * t2 = localtime(&t);
    strftime(tbuf, sizeof(tbuf), "%F %T ", t2); // %F = %Y-%m-%d; %T = %H:%M:%S
    o << tbuf << setw(10) << n << " ";
}

Logger & endl(Logger & logger) {
    logger.o << endl;
    logger.newLine = true;
    return logger;
}

Logger::Logger(ostream & o, int rank) : o (o), newLine(true) {
    if (rank == 0) {
        name = "MASTER";
    } else {
        stringstream s;
        s << "SLAVE" << rank;
        name = s.str();
    }
}
Logger::~Logger() { }
Logger & Logger::operator <<(Logger &(*manipulator)(Logger &)) {
    return manipulator(*this);
}
Logger & Logger::operator <<(const std::string & m) {
    if (newLine) {
        printLineBeginning(o, name);
        newLine = false;
    }
    o << m;
    return *this;
}
Logger & Logger::operator <<(int x) {
    if (newLine) {
        printLineBeginning(o, name);
        newLine = false;
    }
    o << x;
    return *this;
}
Logger & Logger::operator <<(double x) {
    if (newLine) {
        printLineBeginning(o, name);
        newLine = false;
    }
    o << x;
    return *this;
}
