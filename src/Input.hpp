/*
 * input.hpp
 *
 *  Created on: Oct 9, 2013
 *      Author: john
 */
#pragma once
#ifndef INPUT_HPP_
#define INPUT_HPP_

//#include <boost/filesystem.hpp>
#include <sys/stat.h>

#include "LEM3D.hpp"
#include "MiscTools.hpp"

class Input
{
public:
    bool initLogFile ();
    void configFile(char**);
    int mkpath(std::string, mode_t);
    bool isRead;

private:
    std::ifstream infile;
    std::string LogPathName;
    std::string ErrLogPathName;
    std::vector<unsigned>			inUnsigned;
    std::vector<int> 				inInt;
    std::vector<double> 			inDouble;
    std::vector<bool>				inBool;

    void logConfig ();
};

class RTCommand
{
public:
	RTCommand() {};
	RTCommand	(char*	fileName);
	bool isStop ();
	bool update ();
private:
	std::ifstream 	infile;
	std::string		fileName;
};

#endif /* INPUT_HPP_ */
