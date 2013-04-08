#ifndef _BMFHH_H_
#define _BMFHH_H_

#include "PostingOriented_BMW.h"

class DocidOriented_BMW{
private:
	unsigned int* pages;
	inline bool pickCandidate(lptrArray& lps, const float threshold, const int numOfTermsInQ, int& pivot, float& sumOfBLockMaxesUptoPivot) MYNOINLINE;
	inline void handleReal(lptrArray& lps, const int pivot, float sumOfBLockMaxesUptoPivot,PriorityArray<QpResult>&, const float threshold, const int topK) MYNOINLINE;
	inline bool handleFake(lptrArray& lps, const int pivot, const int numOfTermsInQ) MYNOINLINE;

public:
	DocidOriented_BMW(unsigned int* pgs) : pages(pgs) {}
	PriorityArray<QpResult> operator()(lptrArray& lps, const int k);
};

#endif
