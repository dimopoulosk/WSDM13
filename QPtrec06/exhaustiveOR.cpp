/*
 * exhaustiveOR.cpp
 *
 *  Created on: May 17, 2012
 *      Author: constantinos
 */

#include "exhaustiveOR.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"

void ExhaustiveOR::operator() (lptrArray& lps, const int topK, QpResult* res) {
	// initial sorting by did
	lps.sort();

	int smallest_did = lps[0]->did;
	float final_score = 0;

	// initialize results heap
	for ( int i = 0; i < topK; ++i)  {
		res[i].did = -1;
		res[i].score = -1.0;
	}

	while(true) {
		// set round's threshold
		const float threshold = res[topK-1].score;

		// check for termination condition
		if (smallest_did >= CONSTS::MAXD)
			break;

		// initialize final score
		final_score = 0.0f;

		// evaluate all dids with did == smallest_did
		for (int i=0; i<lps.size(); ++i){
			if (lps[i]->did == smallest_did) {
				//PROFILER(CONSTS::EVAL);
				//PROFILER(CONSTS::GETFREQ);
				const float frequency = lps[i]->getFreq();
				const float score = lps[i]->calcScore(frequency,pages[smallest_did]);
				final_score += score;

				// move safely to the next did
				lps[i]->did = lps[i]->nextGEQ( smallest_did + 1 );
				//PROFILER(CONSTS::NEXTGEQ);
			}
		}

		// if calculated score more than threshold, heapify
		if (Fcompare(final_score, threshold)==1) {
			//PROFILER(CONSTS::HEAPIFY);
			int j;
			for (j = topK-2; (j >= 0) && (Fcompare(final_score, res[j].score)==1); j--)
				res[j+1]=res[j];
			res[j+1].setR(smallest_did,final_score);
		}

		// sort list by did
		lps.sort();

		// set new smallest did for the next round
		smallest_did = lps[0]->did;
	} //end while
}
