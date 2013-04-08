/*
 * Wand.h
 *
 *  Created on: May 20, 2012
 *      Author: constantinos
 */

#ifndef WAND_H_
#define WAND_H_

#include "PostingOriented_BMW.h"

class Wand{
private:
	unsigned int* pages;
	inline bool pickCandidate(lptrArray& lps, const float& threshold, int& pivot, const int& threshold_unknown) MYNOINLINE;
	inline void handlePivot(lptrArray& lps, const int& pivot, QpResult* res, const float& threshold, const int& topK, const int& threshold_unknown, bool& update_threshold) MYNOINLINE;

public:
	Wand(unsigned int* pgs) : pages(pgs) {}
	void operator()(lptrArray& lps, const int& k, QpResult* r, const float& threshold);
};

#endif /* WAND_H_ */
