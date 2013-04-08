/*
 * profiling.cpp
 *
 *  Created on: Jan 12, 2012
 *      Author: sergeyn
 */
#include "profiling.h"
#include "globals.h"
profilerC profilerC::GlobalProfiler;

void profilerC::setCosts() {
	profilerC::getInstance().add("mistrust",-42);
	for(size_t i=0; i<100000000; ++i) {
		profilerC::getInstance().start(-42);
		profilerC::getInstance().end(-42);
	}
	cycleCost = timers[-42].getTotal()/1000.0;
	profilerC::getInstance().reset(-42);
}

void profilerC::printReport() const {
	std::map<int,sTimer>::const_iterator it;
	COUT3 << "==============" << Log::endl;
	for ( it=timers.begin() ; it != timers.end(); it++ ) {
		int id = (*it).first;
		const sTimer& t = (*it).second;
		const std::string& cm = (*comments.find(id)).second;
		if(t.getCount()<1)
			continue;

		unsigned long count = t.getCount();
		double total = 	t.getTotal()/1000.0; //in ms
		COUT3 << "[" << id << "]" << cm <<
		" total: " << total-(cycleCost*(t.getCount()/100000000))
		  << "(" << total << ") count: " << count  <<
		" avg: " <<   total/count	<< Log::endl;
	}
}

