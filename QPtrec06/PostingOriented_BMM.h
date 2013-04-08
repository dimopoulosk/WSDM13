/*
 * PostingOriented_BMM.h
 *
 *  Created on: May 22, 2012
 *      Author: constantinos
 */

#ifndef POSTING_ORIENTED_BMM_H_
#define POSTING_ORIENTED_BMM_H_

#include "PostingOriented_BMW.h"

class PostingOriented_BMM{
private:
	unsigned int* pages;
	inline void evaluate(lptrArray& lps, QpResult* res, const float& threshold, const int topK, const int& num_of_essential_lists, int& smallest_did, std::vector <float>& lists_block_max_score, const float& current_prefix_sum_max_score, bool& check_for_new_essential_list, const int& threshold_unknown, const int& comparison, int& next_smallest_did) MYNOINLINE;
	inline bool find_smallest_did_in_essential_set(lptrArray& lps, const int& num_of_essential_lists, int& smallest_did) MYNOINLINE;
	inline void find_new_essential_lists(unsigned int& num_essential_lists, unsigned int& num_of_non_essential_lists, std::vector <float>& prefix_sum_max_score, lptrArray& lps, const float& threshold, bool& essential_list_added);
public:
	PostingOriented_BMM(unsigned int* pgs) : pages(pgs) {}
	void operator()(lptrArray& lps, const int k, QpResult* r, const float& threshold);
};

#endif /* POSTING_ORIENTED_BMM_H_ */
