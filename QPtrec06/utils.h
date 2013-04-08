/*
 * utils.h
 *
 *  Created on: Mar 22, 2012
 *      Author: sergeyn
 */

#ifndef UTILS_H_
#define UTILS_H_
#include <limits>
#include "globals.h"
#include <algorithm>

#define PARANOIC_DEL(p) if(p) delete p; p=0;
#define PARANOIC_DEL_ARR(p) if(p) delete[] p; p=0;

	// Knuth trick for comparing floating numbers
	inline bool FloatEquality( float a, float b) { // check if a and b are equal with respect to the defined tolerance epsilon
	  const float epsilon = 0.0001; // or some other small number (for our case it is sufficient)
	  if ((std::fabs(a-b) <= epsilon*std::fabs(a)) && (std::fabs(a-b) <= epsilon*std::fabs(b))) return true;
	  else return false; }

	inline bool FloatComparison( float a, float b) { // compare two floats ONLY for greater or smaller and NOT for equality
	    if (a > b) return true;
	    else return false; }

	// Arguments a,b and returns a) 0 in case of equality, b) 1 in case a > b, and c) -1 in case a < b
	inline int Fcompare( float a, float b) {
	  if (FloatEquality(a, b)) return 0; // equal
	  else { if (FloatComparison(a, b)) return 1;  // a > b
	         else return -1; }      // a < b
	}

/*
	inline void Get_Percentiles(lptrArray& lps, std::vector<float>& percentiles, const float& threshold) {
		const float frequency = lps[i]->getFreq();
		const float score = lps[i]->calcScore(frequency,pages[smallest_did]);
	}
*/

	inline uint32_t intlog2(const uint32_t x) {
	  uint32_t y;
	  asm ( "\tbsr %1, %0\n"
	      : "=r"(y)
	      : "r" (x)
	  );
	  return y;
	}

	inline unsigned int NearestPowerOf2( unsigned int x ){ //test this first...
	    return 1 << (1+intlog2(x));
	}

	template <typename T> //stream vectors
	std::ostream& operator<< (std::ostream& s, const std::vector<T>& t) {
	        for(auto& it : t)
	                s << it << ' ';
	        return s;
	}

	// Find live areas, given lps and threshold output a vector of liveareas in format <live_area1 = [start,end]> , ... <live_areak = [start,end]>
	//inline void find_live_areas_with_hash(lptrArray& lps, const float threshold, std::vector <int>& live_areas) {
	//}

	template<class T>
	class PriorityArray{
		const int size;
	    std::vector<T> minData;

	public:
	    PriorityArray(const size_t sz=CONSTS::TOPK):size(sz){
	    	minData.resize(size+1);
	    	std::fill(minData.begin(),minData.end(),std::numeric_limits<T>::min());
	    }
	    const std::vector<T>& getV() const { return minData;}
	    const T& head() const { return minData[0]; }
	    void sortData() { std::sort(minData.begin(), minData.end()); std::reverse(minData.begin(), minData.end()); }
	    bool push(const T& val){
	    	//are we better than the smallest top-k member?
	    	if(minData[0]>val)
	    		return false;
	    	minData[0] = val;
			for (unsigned int i = 0; (2*i+1) < size; ){
				unsigned int next = 2*i + 1;
				if (minData[next] > minData[next+1])
					++next;
				if (minData[i] > minData[next])
					std::swap(minData[i], minData[next]);
				else
					break;
				i = next;
			}
			return true;
	    }
	};



#endif /* UTILS_H_ */
