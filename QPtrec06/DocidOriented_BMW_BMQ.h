/*
 * DocidOriented_BMW_BMQ.h
 *
 *  Created on: Jul 31, 2012
 *      Author: constantinos
 */

#ifndef DOCIDORIENTED_BMW_BMQ_H_
#define DOCIDORIENTED_BMW_BMQ_H_

#include "PostingOriented_BMW.h"

class DocidOriented_BMW_BMQ{
private:
	unsigned int* pages;
	inline bool pickCandidate(lptrArray& lps, const float threshold, const int numOfTermsInQ, int& pivot, float& sumOfBLockMaxesUptoPivot) MYNOINLINE;
	inline void handleReal(lptrArray& lps, const int pivot, float sumOfBLockMaxesUptoPivot,PriorityArray<QpResult>&, float threshold, const int topK) MYNOINLINE;
	inline bool handleFake(lptrArray& lps, const int pivot, const int numOfTermsInQ) MYNOINLINE;
public:
	DocidOriented_BMW_BMQ(unsigned int* pgs) : pages(pgs) {}
	PriorityArray<QpResult> operator()(lptrArray& lps, const int k);
	PriorityArray<QpResult> operator()(lptrArray& lps, const int k, const float& theta);
};

#endif /* DOCIDORIENTED_BMW_BMQ_H_ */
