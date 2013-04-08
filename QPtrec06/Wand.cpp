/*
 * Wand.cpp
 *
 *  Created on: May 20, 2012
 *      Author: constantinos
 */

#include "Wand.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"

inline bool Wand::pickCandidate(lptrArray& lps, const float& threshold, int& pivot, const int& threshold_unknown) {
	//PROFILER(CONSTS::CAND);
	// find pivot term until sum of max lists score >= threshold
	float sum_of_list_max_score = 0.0f;
	pivot = -1;
	for (int i = 0; i < lps.size(); ++i)		{
		sum_of_list_max_score += lps.getListMaxScore(i);
		if (Fcompare(sum_of_list_max_score, threshold) >= threshold_unknown) {
			pivot = i;
			break;
		}
	}

	// termination condition: if pivot is -1 or MAXD, no item will make into topK so we can break
	if (pivot == -1 || lps[pivot]->did >= CONSTS::MAXD)
		return false;
	else
		return true;
}

inline void Wand::handlePivot(lptrArray& lps, const int& pivot, QpResult* res, const float& threshold, const int& topK, const int& threshold_unknown, bool& update_threshold) {
	// check alignment + handle real
	if (lps[pivot]->did == lps[0]->did) {
		//PROFILER(CONSTS::ALIGN);

		float final_score = 0.0f;
		const int did = lps[0]->did;

		// we loop through all lists to ensure we pick lists with the same did that appear after the pivot list (in the current sorting)
		for (int i = 0; i < lps.size(); ++i)	{
			if (did == lps[i]->did) {
				//PROFILER(CONSTS::GETFREQ);
				//PROFILER(CONSTS::EVAL);
				const float frequency = lps[i]->getFreq();
				const float score = lps[i]->calcScore(frequency, pages[did]);
				final_score += score;
				// move safely to did + 1, since there is no chance this did to be evaluated
				lps[i]->did = lps[i]->nextGEQ( did + 1);
				//PROFILER(CONSTS::NEXTGEQ);
			}
		}

		// check if this did can make it to the topk and if so heapify
		if (Fcompare(final_score, threshold) >= threshold_unknown) {
			//PROFILER(CONSTS::HEAPIFY);
			update_threshold = true;
			int j;
			for (j = topK-2; (j >= 0) && (Fcompare(final_score, res[j].score) == 1); j--)
				res[j+1] = res[j];
			res[j+1].setR(did, final_score);
		}

		// reorder list
		//PROFILER(CONSTS::SORT);
		// TOFIX
		//for (int i = pivot; i>=0; i--)
      	//	lps.popdown(i);
		lps.sort();
	}
	else {
		// we have misalignment so we find the term-list with the smallest list length (according to wand paper)
		// if misalignment -> handle fake
		//PROFILER(CONSTS::MISALIGNED);
		int smallest_list_length_list = CONSTS::MAXD;
		int smallest_list_length_list_index = -1;
		for (int i = 0; i < pivot; ++i) {
			if (lps[i]->lengthOfList < smallest_list_length_list && lps[i]->did != lps[pivot]->did ) {
				smallest_list_length_list = lps[i]->lengthOfList;
				smallest_list_length_list_index = i;
			}
		}
		// move safely the selected list up to pivot's did
		lps[smallest_list_length_list_index]->did = lps[smallest_list_length_list_index]->nextGEQ( lps[pivot]->did );
		//PROFILER(CONSTS::NEXTGEQ);

		// re-sort the list
		//PROFILER(CONSTS::SORT);
		lps.popdown(smallest_list_length_list_index);
	}
}

void Wand::operator() (lptrArray& lps, const int& topK, QpResult* res, const float& threshold){
	int num_query_terms = lps.size();
	int pivot;
	float sum_of_block_maxes_up_to_pivot = 0.0f;
	float current_threshold = threshold;
	const int threshold_unknown = (Fcompare(threshold, 0.0f) == 0) ? 1 : 0; // treat as bool variable
	bool update_threshold = false;

	// initialize results heap
	for ( int i = 0; i < topK; ++i)  {
		res[i].did = CONSTS::MAXD + 1;
		res[i].score = -1.0;
	}

	// initial sorting by did
	//PROFILER(CONSTS::SORT);
	lps.sort();

	while(true) {
		// if threshold is not known update the current_threshold
		if ((update_threshold) && (threshold_unknown == 1)) {
			current_threshold = res[topK-1].score;
			update_threshold = false;
		}

		// pick candidate and break if pivot did out of docid space
		if(! pickCandidate(lps, current_threshold, pivot, threshold_unknown))
			break;
		else
			// handle real and fake
			handlePivot(lps, pivot, res, current_threshold, topK, threshold_unknown, update_threshold);
	}
}
