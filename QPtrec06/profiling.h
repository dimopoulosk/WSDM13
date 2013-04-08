#ifndef _PROFILING_H_
#define _PROFILING_H_

#include <sys/time.h>
#include <map>
#include <iostream>
#include <vector>

class sTimer {
private:
	timeval startT;
	double total; 
	unsigned long count;
	
public:
		sTimer() { reset();	}

		double getTotal() const { return total; }
		unsigned long getCount() const { return count;}
				
		inline void start() { gettimeofday(&startT,0); }
					
		void reset() { total = count = 0; }

		inline double end() {
			timeval end;
			gettimeofday(&end,0);
			double t0=(startT.tv_sec*1000000.0)+startT.tv_usec;
			double t1=(end.tv_sec*1000000.0)+end.tv_usec;

			double delta=(t1-t0);
			total+= delta;
			++count;
			return  delta;
		}
};


class profilerC { //single threaded
private:
	std::map<int,std::string> comments;	
	std::map<int,sTimer> timers;
	std::map<int,unsigned long> cyclesAtStart;
	std::vector<unsigned long int> counters;

	double cycleCost;
	unsigned long cycles;
	
	static profilerC GlobalProfiler;
	profilerC():cycles(0) { /*setCosts();*/}

	void setCosts();
	
public:
	static inline profilerC& getInstance() { return GlobalProfiler; }

	size_t initNewCounter() {
		counters.push_back(0);
		return counters.size()-1;
	}

	inline void stepCounter(size_t cnt, const unsigned long int inc=1) { counters[cnt]+=inc; }
	inline  unsigned long int getCounter(size_t cnt) { return counters[cnt]; }

	void add(const std::string& comment, int id) {
		comments[id] = comment;
		timers[id] = sTimer();
	} 
	
	void reset(int id) { timers[id].reset();}

	inline void start(int id) {
		#ifdef TIMING
		//cyclesAtStart[id]=cycles;
		timers[id].start();
		#endif
	}

	inline
#ifdef TIMING
	double
#else
	void
#endif
	end(int id) {
		#ifdef TIMING
		//double deduct = cycleCost*(cyclesAtStart[id]-cycles);
		//++cycles;
		return timers[id].end();
		#endif
	}
	void mistrustPrinter();
	void printReport() const;
};

#endif
