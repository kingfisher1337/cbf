/* 
 * File:   Logger.hpp
 * Author: kingfisher1337
 *
 * Created on January 2, 2014, 9:54 PM
 */

#ifndef LOGGER_HPP
#define	LOGGER_HPP

#include <string>

class Logger;

Logger & endl(Logger & logger);

class Logger {
    friend Logger & endl(Logger & logger);
private:
    std::ostream & o;
    std::string name;
    bool newLine;
public:
    Logger(std::ostream & o, int rank);
    ~Logger();
public:
    Logger & operator <<(Logger &(*manipulator)(Logger &));
    Logger & operator <<(const std::string & m);
    Logger & operator <<(int x);
    Logger & operator <<(double x);
};

#endif	/* LOGGER_HPP */
