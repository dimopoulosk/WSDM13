/*
 * DocidOriented_BMW.cpp
 *
 *  Created on: Jan 14, 2012
 *      Author: constantinos
 */


#include "DocidOriented_BMW.h"
#include "globals.h"
#include "profiling.h"

inline bool DocidOriented_BMW::pickCandidate(lptrArray& lps, const float threshold, const int numOfTermsInQ, int& pivot, float& sumOfBLockMaxesUptoPivot) {
	pivot = -1; 		// find the pivot term
	//select pivot candidate like you did in wand -- according to list-max-scores
	float sumOfListMaxScores = 0.0;
	for( int i = 0; i < numOfTermsInQ; ++i)		{
		sumOfListMaxScores += lps.getListMaxScore(i);
		if( sumOfListMaxScores > threshold ) {
			//small optimization: if ALL dids after pivot-list are equal, it safe to move the pivot down
			while( i + 1 < numOfTermsInQ && lps[i + 1]->did == lps[i]->did)
				++i;

			pivot = i;
			break;
		}
	}

	lptr* lpspivot = lps[pivot];
	if(pivot == -1 || lpspivot->did >= CONSTS::MAXD) // pivot is -1 or MAXD, no item will make into topK
		return false;

	//smallestMaxDid is used in handleFake and it is the smallest maxDidPerBlock[] + 1 of the lists up to pivot
	//smallestMaxDid = lpspivot->gen.skipToNextBlock(lpspivot->did) - 1 > CONSTS::MAXD ? CONSTS::MAXD : lpspivot->gen.skipToNextBlock(lpspivot->did);
	//this is overapprox. of the real sum of scores and used in eval. in hadnleReal
	sumOfBLockMaxesUptoPivot = lps[pivot]->gen.getMaxScore(lpspivot->did); // lpspivot->getMaxScoreOfDeepBlock();

	for(int i = 0; i < pivot; ++i) 	{
		sumOfBLockMaxesUptoPivot +=  lps[i]->gen.getMaxScore(lpspivot->did);
	}
	return true;
}


inline void DocidOriented_BMW::handleReal(lptrArray& lps, const int pivot, float sumOfBLockMaxesUptoPivot,PriorityArray<QpResult>&       resultsHeap, const float threshold, const int topK) {
	// evaluate this did only when did is equal to lps[0]->did
	lptr* lpspivot = lps[pivot];
	const int& did = lpspivot->did;
	float sumOfDidScores = 0.0f;
	if( did == lps[0]->did ) { // this dude needs to be evaluated!!
		for(int i = 0 ; i < pivot + 1; ++i)	{
			//PROFILER(CONSTS::GETFREQ);
			//PROFILER(CONSTS::EVAL);
			const float f = lps[i]->getFreq();
			const float scoret = lps[i]->calcScore(f,pages[did]);
			sumOfDidScores += scoret;
			sumOfBLockMaxesUptoPivot -= ( lps[i]->gen.getMaxScore(did)/* was => lps[i]->getMaxScoreOfDeepBlock() => was score[ lps[i]->block ]*/ - scoret );
			//early termination after substr. deltas of real scores from approximated upper bound
			if( !(sumOfBLockMaxesUptoPivot > threshold) )
				goto mvptrs;
		}

		if (sumOfDidScores > threshold) {  //heapify the new result
			QpResult tmp(did,sumOfDidScores);
			if(resultsHeap.push(tmp))
				;	//PROFILER(CONSTS::HEAPIFYOK);
		} //end heapify
mvptrs:
		for(int i = 0 ;i < pivot + 1; ++i)  //move pointers
			lps[i]->did = lps[i]->nextGEQ( did + 1);

		for(int i = pivot ; i > -1; i--)//reorder list
			lps.popdown(i);
	} else {
		int least_ft = lps[0]->lengthOfList;
		int least_i = 0;
		for(int i = 1; i < pivot ; ++i)	{
			if( lps[i]->lengthOfList < least_ft &&  lps[i]->did != did ) {
				least_ft = lps[i]->lengthOfList;
				least_i = i;
			}
		}
		lps[least_i]->did = lps[least_i]->nextGEQ( did );
		lps.popdown(least_i);
	}
}

inline  bool DocidOriented_BMW::handleFake(lptrArray& lps, const int pivot, const int numOfTermsInQ) {
	//PROFILER(CONSTS::BLOCKFAIL);
	int least_ft = lps[0]->lengthOfList;
	int least_i = 0;
	for(int i = 1; i < pivot + 1; ++i) { //find shortest list
		if( lps[i]->lengthOfList < least_ft) {
			least_ft = lps[i]->lengthOfList;
			least_i = i;
		}
	}

	//smallestMaxDid is the smallest maxDidPerBlock[] + 1 of the lists up to pivot
	unsigned int smallestMaxDid = lps[pivot]->gen.skipToNextBlock(lps[pivot]->did) - 1 > CONSTS::MAXD ? CONSTS::MAXD : lps[pivot]->gen.skipToNextBlock(lps[pivot]->did);
	for(int i = 0; i < pivot; ++i) 	{
		int tempd = lps[i]->gen.skipToNextBlock(lps[pivot]->did)-1; // shallowMove(lps[pivot]->did);
		if( smallestMaxDid >  tempd)
			smallestMaxDid = tempd;
	}

	int did = (( pivot < numOfTermsInQ - 1) && smallestMaxDid > lps[ pivot + 1 ]->did)? lps[ pivot + 1 ]->did : smallestMaxDid;

	if( did <= lps[pivot]->did)
		did = lps[pivot]->did + 1;
	if(did > CONSTS::MAXD)
		return false;
	lps[least_i]->did = lps[least_i]->nextGEQ( did );
	lps.popdown(least_i);
	return true;
}

PriorityArray<QpResult> DocidOriented_BMW::operator()( lptrArray& lps, const int topK) {
	int numOfTermsInQ;
		float  sumOfBLockMaxesUptoPivot, scoret;
		int f,pivot,smallestMaxDid,tempd,least_i;
		unsigned int least_ft;
		lptr* lpspivot;

		numOfTermsInQ= lps.size();
		PriorityArray<QpResult> resultsHeap(topK); // initialize result heap
		lps.sort(); //initial sorting by did
		while(true) {
			const float threshold =resultsHeap.head().score;

			if(! pickCandidate(lps,threshold,numOfTermsInQ,pivot,sumOfBLockMaxesUptoPivot))
				break;

			if( sumOfBLockMaxesUptoPivot > threshold ) {
				handleReal(lps,pivot,sumOfBLockMaxesUptoPivot,resultsHeap,threshold,topK);
			}
			else {
				if(! handleFake(lps, pivot,numOfTermsInQ))
					break;
			}
		} //end while
	return resultsHeap;
}
