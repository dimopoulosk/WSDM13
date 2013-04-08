#ifndef _BMF_H_
#define _BMF_H_

#include "pfor.h"
#include "ListIterator.h"

// Result structure
struct QpResult {
	unsigned int did;
	float score;

	inline QpResult& operator=(const QpResult& rhs) {
		did=rhs.did;
		score=rhs.score;
		return *this;
	}

	inline bool operator>(const QpResult& rhs) const { return score>rhs.score; }
	inline bool operator<(const QpResult& rhs) const { return score<rhs.score; }
	QpResult(unsigned int d=CONSTS::MAXD+1, float s=-1.0) { did=d; score=s;}
	inline void setR(unsigned int d, float s) { did=d; score=s;}
};

class PostingOriented_BMW {
private:
	unsigned int* pages;
	inline bool pickCandidate(lptrArray& lps, const float threshold, const int numOfTermsInQ, int& pivot, int& smallestMaxDid, float& sumOfBLockMaxesUptoPivot) MYNOINLINE;
	inline void handleReal(lptrArray& lps, const int pivot, float sumOfBLockMaxesUptoPivot,QpResult* res, const float threshold, const int topK) MYNOINLINE;
	inline bool handleFake(lptrArray& lps, const int pivot, const int smallestMaxDid, const int numOfTermsInQ) MYNOINLINE;

public:
	PostingOriented_BMW(unsigned int* pgs) : pages(pgs) {}
	void operator()(lptrArray& lps, const int k, QpResult* r);
};

#endif
