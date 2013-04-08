/*
 * Maxscore.cpp
 *
 *  Created on: May 14, 2012
 *      Author: constantinos
 */

#include "Maxscore.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"
#include <vector>

inline bool Maxscore::find_smallest_did_in_essential_set(lptrArray& lps, const int& num_of_essential_lists, int& smallest_did) {
	smallest_did = lps[num_of_essential_lists]->did;
	for (int i = num_of_essential_lists; i < lps.size(); ++i) {
		if (lps[i]->did < smallest_did)
			smallest_did = lps[i]->did;
	}

	// check if out of range
	if (smallest_did >= CONSTS::MAXD)
		return false;
	else
		return true;
}

inline void Maxscore::evaluate(lptrArray& lps, QpResult* res, const float& threshold, const int topK, const int& num_of_essential_lists, int& smallest_did, std::vector <float>& lists_maxscore, const float& current_prefix_sum_max_score, bool& check_for_new_essential_list, const int& threshold_unknown, const int& comparison, int& next_smallest_did) {
	bool failure = false;
	float final_score = 0.0f;
    next_smallest_did = CONSTS::MAXD + 1;

	// evaluate dids == to smallest_did in the essential lists
	//for (int i=num_of_essential_lists; i<lps.size(); i++) {         // In order version
	for (int i = lps.size()-1; i>=num_of_essential_lists; --i) {      // Reverse order version
		if (smallest_did == lps[i]->did) {
			//PROFILER(CONSTS::EVAL);
			//PROFILER(CONSTS::GETFREQ);
			//PROFILER(CONSTS::ESSENTIAL);
			const float frequency = lps[i]->getFreq();
			const float score = lps[i]->calcScore(frequency,pages[smallest_did]);
			final_score += score;
			lps[i]->did = lps[i]->nextGEQ( smallest_did + 1 );
			//PROFILER(CONSTS::NEXTGEQ);
		}

		// pick next smallest did
		if (lps[i]->did < next_smallest_did)
			next_smallest_did = lps[i]->did;
	}

	// early termination = prefix + final score of essential lists
	float early_termination = current_prefix_sum_max_score + final_score;

	// if early termination check is passed, evaluate smallest did in the non-essential set
	if (! (Fcompare(early_termination, threshold) <= comparison)) {
		// evaluate dids == to smallest_did in the non essential lists
		//for (int i = 0; i < num_of_essential_lists; ++i)	{       // In order version
		for (int i=num_of_essential_lists-1; i>=0; --i)	{           // Reverse order version
			// move pointers if needed
			if (lps[i]->did < smallest_did) {
				//PROFILER(CONSTS::NEXTGEQ1);
				lps[i]->did = lps[i]->nextGEQ(smallest_did);
				//PROFILER(CONSTS::NEXTGEQ);
			}

			// check if evaluation is needed
			if (smallest_did == lps[i]->did) {
				//PROFILER(CONSTS::EVAL);
				//PROFILER(CONSTS::GETFREQ);
				//PROFILER(CONSTS::NONESSENTIAL);
				const float frequency = lps[i]->getFreq();
				const float score = lps[i]->calcScore(frequency,pages[smallest_did]);
				final_score += score;
				early_termination -= (lists_maxscore.at(i) - score);
			} else
				early_termination -= lists_maxscore.at(i);

			// early termination
			if (Fcompare(early_termination, threshold) <= comparison) {
				failure = true;
				//PROFILER(CONSTS::EARLYTERMINATION2);
				break;
			}
		}

		// if not failure, heapify new result
		if ((!failure) && (Fcompare(final_score, threshold) >= threshold_unknown)) {
			//PROFILER(CONSTS::HEAPIFY);
			check_for_new_essential_list = true;
			int j;
			for (j = topK-2; (j >= 0) && (Fcompare(final_score, res[j].score)==1); j--)
				res[j+1]=res[j];
			res[j+1].setR(smallest_did,final_score);
		}
	} else {  // togo
		//PROFILER(CONSTS::EARLYTERMINATION1);
	}
}

void Maxscore::find_new_essential_lists(unsigned int& num_essential_lists, unsigned int& num_of_non_essential_lists, std::vector <float>& prefix_sum_max_score, lptrArray& lps, const float& threshold, bool& essential_list_added, float& parameter) {
		while ((num_essential_lists<lps.size()-1) && (Fcompare(prefix_sum_max_score.at(num_essential_lists), threshold*parameter) == -1)) {   // loop until more than threshold
			++num_essential_lists;        // increase number of essential lists
			--num_of_non_essential_lists; // decrease number of non-essential lists
			essential_list_added = true;
		}
}

void Maxscore::operator() (lptrArray& lps, const int topK, QpResult* res, const float& threshold){
	unsigned int num_essential_lists = 0;
	unsigned int num_of_non_essential_lists = lps.size();
	float total_maxscore_sum = 0.0f;
	float current_threshold = threshold;
	float current_prefix_sum_max_score = 0.0f;
	float parameter = 1.0f;
	const int threshold_unknown = (Fcompare(threshold, 0.0f) == 0) ? 1 : 0; // treat as bool variable
	const int comparison = -1;
	bool check_for_new_essential_list = false;
	bool essential_list_added = false;
	std::vector <float> lists_maxscore;
	std::vector <float> prefix_sum_max_score;
	lists_maxscore.reserve(lps.size());
	prefix_sum_max_score.reserve(lps.size()-1);

	// initialize results heap
	for (int i = 0; i < topK; ++i)  {
		res[i].did = CONSTS::MAXD + 1;
		res[i].score = -1.0;
	}

	// initial sorting of lists by max score
	lps.sortbyscore();
	int smallest_did = lps[0]->did;
	int next_smallest_did = CONSTS::MAXD + 1;

	// initialize maxscore array
	for (int i=0; i<lps.size(); i++) {
		lists_maxscore.push_back(lps.getListMaxScore(i));
		total_maxscore_sum += lps.getListMaxScore(i); // useless to check!
	}

	// initialize prefix sum max score array that is used in find_essential_lists
	prefix_sum_max_score.push_back(lps.getListMaxScore(0));
	for (int i=1; i<lps.size()-1; i++)
		prefix_sum_max_score.push_back(prefix_sum_max_score.at(i-1) + lps.getListMaxScore(i));

	// if threshold is known find essential lists
	if (threshold_unknown==0) {
		find_new_essential_lists(num_essential_lists, num_of_non_essential_lists, prefix_sum_max_score, lps, current_threshold, essential_list_added, parameter);
		essential_list_added = false;
	}

	// find first smallest did
	find_smallest_did_in_essential_set(lps, num_essential_lists, smallest_did);

	// set it for the first time only
	if (threshold_unknown==1)
		next_smallest_did = smallest_did;

	while(true) {
		// if threshold is not known and heapify occured, update the current_threshold, find if a new list can be added to the essential set
		if ((threshold_unknown == 1) && (check_for_new_essential_list)) {
			// update threshold
			current_threshold = res[topK-1].score;
			// find essential lists
			find_new_essential_lists(num_essential_lists, num_of_non_essential_lists, prefix_sum_max_score, lps, current_threshold, essential_list_added, parameter);
			check_for_new_essential_list = false;

			// if essential list added, check for new smallest did if needed
			if (!essential_list_added)
				smallest_did = next_smallest_did; // we have already found the smallest did from the previous round
			else {
				essential_list_added = false;
				// if new smallest did out of range break, otherwise return smallest_did in this round
				if (! find_smallest_did_in_essential_set(lps, num_essential_lists, smallest_did))
					break;
			}
		} else
			smallest_did = next_smallest_did; // we have already found the smallest did from the previous round

		// set correct current prefix sum max score
		current_prefix_sum_max_score = num_essential_lists > 0 ? prefix_sum_max_score.at(num_essential_lists-1) : 0.0f;

		// evaluate
		evaluate(lps, res, current_threshold, topK, num_essential_lists, smallest_did, lists_maxscore, current_prefix_sum_max_score, check_for_new_essential_list, threshold_unknown, comparison, next_smallest_did);

		// check if next smallest did is out of docid range
		if (next_smallest_did >= CONSTS::MAXD)
			break;
	}
}
