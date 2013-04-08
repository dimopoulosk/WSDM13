/*
 * PostingOriented_BMM_NLB.h
 *
 *  Created on: Jul 28, 2012
 *      Author: constantinos
 */

#ifndef POSTING_ORIENTED_BMM_NLB_H_
#define POSTING_ORIENTED_BMM_NLB_H_

#include "PostingOriented_BMW.h"

class PostingOriented_BMM_NLB{
private:
	unsigned int* pages;
	inline void evaluate(lptrArray& lps, QpResult* res, const float& threshold, const int topK, const int& num_of_essential_lists, int& smallest_did, std::vector <float>& lists_block_max_score, const float& current_prefix_sum_max_score, bool& check_for_new_essential_list, const int& /*threshold_unknown*/comparison_heapify, const int& comparison, std::vector <float>& lists_maxscore, const bool& skip_to_next_live_block_enabled) MYNOINLINE;
	inline bool find_smallest_did_in_essential_set(lptrArray& lps, const int& num_of_essential_lists, int& smallest_did) MYNOINLINE;
	inline void find_new_non_essential_lists(unsigned int& num_essential_lists, unsigned int& num_of_non_essential_lists, std::vector <float>& prefix_sum_max_score, lptrArray& lps, const float& threshold, bool& essential_list_added);

	// Finds next live block
	// threshold can be either current threshold (tighter) or dif(threshold - sum(non-essential maxscores))
	// smallest_did is the current did that failed in Filter that allows safe block skipping
	inline void skip_to_next_live_block(lptrArray& lps, const float& threshold, const int& num_essential_lists, const int& check_up_to_list, const int& find_next_live_block_from_list, int& smallest_did, const bool& skip_to_next_live_block_enabled);
public:
	PostingOriented_BMM_NLB(unsigned int* pgs) : pages(pgs) {}
	void operator()(lptrArray& lps, const int& k, QpResult* res, const float& threshold, const int& threshold_knowledge);
};
#endif /* POSTING_ORIENTED_BMM_NLB_H_ */
