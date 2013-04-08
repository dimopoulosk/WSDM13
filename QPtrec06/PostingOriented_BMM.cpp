/*
 * PostingOriented_BMM.cpp
 *
 *  Created on: May 22, 2012
 *      Author: constantinos
 */

/* Algorithm description:
 * PostingOriented_BMM: (evaluate method)
 * 1. evaluate essential set
 * 2. check whether (sum(essential set score) + prefix > threshold)
 * 3. if so, get block maxes of non-essential set
 * 4. then check whether (sum(essential set score) + sum(block maxes of non-essential set) > threshold)
 * 5. if so refine non-essential set by evaluating and having early termination
 *
 * Expected avg qp time for query log: 32-33ms for GOV2
 */

#include "PostingOriented_BMM.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"
#include <vector>

inline bool PostingOriented_BMM::find_smallest_did_in_essential_set(lptrArray& lps, const int& num_of_essential_lists, int& smallest_did) {
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

inline void PostingOriented_BMM::evaluate(lptrArray& lps, QpResult* res, const float& threshold, const int topK, const int& num_of_essential_lists, int& smallest_did, std::vector <float>& lists_block_max_score, const float& current_prefix_sum_max_score, bool& check_for_new_essential_list, const int& threshold_unknown, const int& comparison, int& next_smallest_did) {
	float final_score = 0.0f;
    float block_score_sum = 0.0f;
    next_smallest_did = CONSTS::MAXD + 1;

    // candidate population
   //PROFILER(CONSTS::CAND);

    // Evaluate essential set
    for (int i = lps.size()-1; i>=num_of_essential_lists; --i) {
		if (smallest_did == lps[i]->did) {
			//PROFILER(CONSTS::EVAL);
			//PROFILER(CONSTS::GETFREQ);
			//PROFILER(CONSTS::ESSENTIAL);
			const float frequency = lps[i]->getFreq();
			const float score = lps[i]->calcScore(frequency,pages[smallest_did]);
			final_score += score;
			//PROFILER(CONSTS::NEXTGEQ);
			lps[i]->did = lps[i]->nextGEQ( smallest_did + 1 );
		}

		// pick next smallest did
		if (lps[i]->did < next_smallest_did)
			next_smallest_did = lps[i]->did;
	}

	// Early termination 2: if actual score of essential set + sum of list maxscore of non-essential, less than threshold
	if (! (Fcompare(final_score + current_prefix_sum_max_score, threshold) <= comparison)) {

		//PROFILER(CONSTS::EARLYTERMINATION1);
		// Get Block Maxes in the non-essential set
		for (int i=num_of_essential_lists-1; i>=0; --i)	{
			// move shallow pointers to the current smallest did
			//PROFILER(CONSTS::SHALLOWPOINTER1);
			lps[i]->DEPshallowMove(smallest_did);

			// get and store block max score
			//PROFILER(CONSTS::GETBLOCKSCORE1);
			lists_block_max_score.at(i) = lps[i]->getMaxScoreOfBlock();

			// add to sum of non-ess. block maxes
			block_score_sum += lists_block_max_score.at(i);
		}

		// total upper bound of score: sum of non-essential block maxes + essential set actual score
		block_score_sum += final_score;
		//bool fail = false;

		// Early termination 3: if essential actual score + sum of block scores of non-essential set, less than threshold
		if (! (Fcompare(block_score_sum, threshold) <= comparison)) {
			//PROFILER(CONSTS::EARLYTERMINATION3);

			// Evaluate non-essential set
			for (int i=num_of_essential_lists-1; i>=0; --i)	{
				// move pointers if needed
				if (lps[i]->did < smallest_did) {
					//PROFILER(CONSTS::NEXTGEQ1);
					lps[i]->did = lps[i]->nextGEQ( smallest_did );
				}

				// check if evaluation is needed
				if (smallest_did == lps[i]->did) {
					//PROFILER(CONSTS::EVAL);
					//PROFILER(CONSTS::GETFREQ);
					//PROFILER(CONSTS::NONESSENTIAL);
					const float frequency = lps[i]->getFreq();
					const float score = lps[i]->calcScore(frequency,pages[smallest_did]);
					final_score += score;
					block_score_sum -= (lists_block_max_score.at(i) - score);
				} else
					block_score_sum -= lists_block_max_score.at(i);

				// Early terminate during non-essential set max block scores refinement
				if (Fcompare(block_score_sum, threshold) <= comparison) {
					//PROFILER(CONSTS::EARLYTERMINATION5);
					//fail = true;
					break;
				}
			}

			// if smallest did passed all filters, heapify new result
			if (Fcompare(final_score, threshold) >= threshold_unknown) {
				//PROFILER(CONSTS::HEAPIFY);
				check_for_new_essential_list = true;
				int j;
				for (j = topK-2; (j >= 0) && (Fcompare(final_score, res[j].score)==1); j--)
					res[j+1]=res[j];
				res[j+1].setR(smallest_did,final_score);
			}/* else {
				if (!fail) {
					PROFILER(CONSTS::GETBLOCKSCORE);
					//std::cout << final_score << " , theta: " << threshold << ", blocksum: " << block_score_sum << std::endl;
				}
			}
*/		} //else
		//	PROFILER(CONSTS::EARLYTERMINATION4);
	} //else
		//PROFILER(CONSTS::EARLYTERMINATION2);
}

void PostingOriented_BMM::find_new_essential_lists(unsigned int& num_essential_lists, unsigned int& num_of_non_essential_lists, std::vector <float>& prefix_sum_max_score, lptrArray& lps, const float& threshold, bool& essential_list_added) {
		while ((num_essential_lists<lps.size()-1) && (Fcompare(prefix_sum_max_score.at(num_essential_lists), threshold) == -1)) {   // loop until more than threshold
			++num_essential_lists;        // increase number of essential lists
			--num_of_non_essential_lists; // decrease number of non-essential lists
			essential_list_added = true;
		}
}

void PostingOriented_BMM::operator() (lptrArray& lps, const int topK, QpResult* res, const float& threshold){
	unsigned int num_essential_lists = 0;
	unsigned int num_of_non_essential_lists = lps.size();
	float total_maxscore_sum = 0.0f;
	float current_threshold = threshold;
	float current_prefix_sum_max_score = 0.0f;
	const int threshold_unknown = (Fcompare(threshold, 0.0f) == 0) ? 1 : 0; // treat as bool variable
	const int comparison = -1;
	bool check_for_new_essential_list = false;
	bool essential_list_added = false;
	std::vector <float> lists_maxscore;
    std::vector <float> lists_block_max_score;
	std::vector <float> prefix_sum_max_score;
	lists_maxscore.reserve(lps.size());
	lists_block_max_score.reserve(lps.size());
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
		lists_block_max_score.push_back(0.0f);
		total_maxscore_sum += lps.getListMaxScore(i);
	}

	// initialize prefix sum max score array that is used in find_essential_lists
	prefix_sum_max_score.push_back(lps.getListMaxScore(0));
	for (int i=1; i<lps.size()-1; i++)
		prefix_sum_max_score.push_back(prefix_sum_max_score.at(i-1) + lps.getListMaxScore(i));

	// if threshold is known find essential lists
	if (threshold_unknown==0) {
		find_new_essential_lists(num_essential_lists, num_of_non_essential_lists, prefix_sum_max_score, lps, current_threshold, essential_list_added);
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
			find_new_essential_lists(num_essential_lists, num_of_non_essential_lists, prefix_sum_max_score, lps, current_threshold, essential_list_added);
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
		evaluate(lps, res, current_threshold, topK, num_essential_lists, smallest_did, lists_block_max_score, current_prefix_sum_max_score, check_for_new_essential_list, threshold_unknown, comparison, next_smallest_did);

		// check if next smallest did is out of docid range
		if (next_smallest_did >= CONSTS::MAXD)
			break;
	}
}
