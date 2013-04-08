/*
 * globals.cpp
 *
 *  Created on: Feb 17, 2012
 *      Author: sergeyn
 */



#include "globals.h"
#include <cmath>
#ifdef DEBUG
int Log::logger::GLOBAL_VERBOSITY(Log::VALL);
#else
int Log::logger::GLOBAL_VERBOSITY(Log::VOUTPUT);
#endif

void Log::setGlobalVerbosityForAllLoggers(int v) { Log::logger::setGlobalV(v); }
Log::logger& Log::endl(Log::logger& stream) { //stream flush
	stream.flush();
	stream.setV(0);
	return stream;
}

Log::logger& Log::fatal(Log::logger& stream) {
	stream.flush();
	throw(1);
}
