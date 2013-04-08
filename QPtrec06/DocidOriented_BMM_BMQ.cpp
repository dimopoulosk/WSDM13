/*
 * DocidOriented_BMM_BMQ.cpp
 *
 *  Created on: Jul 24, 2012
 *      Author: constantinos
 */

#include "DocidOriented_BMM_BMQ.h"
#include "globals.h"
#include "profiling.h"
#include "utils.h"
#include <vector>

inline bool DocidOriented_BMM_BMQ::find_smallest_did_in_essential_set(lptrArray& lps, const int& num_of_essential_lists, int& smallest_did) {
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

inline void DocidOriented_BMM_BMQ::evaluate(lptrArray& lps, QpResult* res, const float& threshold, const int topK, const int& num_of_essential_lists, int& smallest_did, std::vector <float>& lists_block_max_score, const float& current_prefix_sum_max_score, bool& check_for_new_essential_list, const int& /*threshold_unknown*/comparison_heapify, const int& comparison, std::vector <float>& lists_maxscore, const bool& skip_to_next_live_block_enabled/*, std::vector< std::vector<float> >& Quantization_lookup_table*//*, std::vector<float>& Quantization_lookup_table*/) {
	float final_score = 0.0f;
    float block_score_sum = 0.0f;
    float essential_block_sum = 0.0f;
    float non_essential_block_sum = 0.0f;
    float early_termination = 0.0f;
    int check_up_to_list = 0;   // needed to keep track of list that block refinement failed and hence for skipping
    bool fail = false; //togo

    // candidate population
    //PROFILER(CONSTS::CAND);

    // sum the block maxes of essential set
    for (int i = lps.size()-1; i>=num_of_essential_lists; --i) {
//    	lists_block_max_score.at(i) = lps[i]->gen.getMaxScore(smallest_did);
//    	essential_block_sum += lists_block_max_score.at(i);
    	//essential_block_sum += Quantization_lookup_table[i*(CONSTS::Quantization+1)+(unsigned int) lps[i]->gen.get_Quantized_Score(smallest_did)];
    	//essential_block_sum += Quantization_lookup_table[i][(unsigned int) lps[i]->gen.get_Quantized_Score(smallest_did)];
    	essential_block_sum += lps[i]->quantile*(float) lps[i]->gen.get_Quantized_Score(smallest_did);
    }

    // Filter 1: sum(essential block maxes) + sum(non-essential maxscores)
    if (! (Fcompare(essential_block_sum + current_prefix_sum_max_score, threshold) <= comparison)) {
    	//PROFILER(CONSTS::EARLYTERMINATION0);

    	// used temporarily
    	early_termination = essential_block_sum + current_prefix_sum_max_score;

    	// get block maxes of non-essential set
		for (int i = num_of_essential_lists-1; i>=0; --i) {
			//lists_block_max_score.at(i) += Quantization_lookup_table[i*(CONSTS::Quantization+1)+(unsigned int) lps[i]->gen.get_Quantized_Score(smallest_did)];
			lists_block_max_score[i] = lps[i]->quantile*(float) lps[i]->gen.get_Quantized_Score(smallest_did);
			//lists_block_max_score.at(i) = Quantization_lookup_table[i][(unsigned int) lps[i]->gen.get_Quantized_Score(smallest_did)];
			non_essential_block_sum += lists_block_max_score[i];
			early_termination -= (lists_maxscore[i] - lists_block_max_score[i]);

			// early termination
			if (Fcompare(early_termination, threshold) <= comparison) {
				check_up_to_list = i; // will be used for skipping
				break;
			}
		}

		// Filter 2: sum(essential set block maxes) + sum(non-essential set block maxes)
		if (! (Fcompare(essential_block_sum + non_essential_block_sum, threshold) <= comparison)) {
			//PROFILER(CONSTS::EARLYTERMINATION1);

			// early termination
//			early_termination = essential_block_sum + non_essential_block_sum;

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
					// togo
/*					early_termination -= (lists_block_max_score.at(i) - score);
				} else
					early_termination -= lists_block_max_score.at(i);

				// early termination
				if (Fcompare(early_termination, threshold) <= comparison) {
					check_up_to_list = i;
					break;
				}
*/				// CAN HAVE EARLY TERMINATION
			}

			block_score_sum = final_score + non_essential_block_sum;
			// Early termination 2: if actual score of essential set + sum of list maxscore of non-essential, less than threshold
			// Filter 3: essential score + sum of block maxes of non-essential set
			if (! (Fcompare(block_score_sum, threshold) <= comparison)) {
				//PROFILER(CONSTS::EARLYTERMINATION3);
				bool fail = false;
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
							block_score_sum -= (lists_block_max_score[i] - score);
						} else
							block_score_sum -= lists_block_max_score[i];

						// Early terminate during non-essential set max block scores refinement
						if (Fcompare(block_score_sum, threshold) <= comparison) {
							//PROFILER(CONSTS::EARLYTERMINATION5);
							fail = true;
							break;
						}
					}

					// if smallest did passed all filters, heapify new result
					if (Fcompare(final_score, threshold) >= comparison_heapify) {
						//PROFILER(CONSTS::HEAPIFY);
						check_for_new_essential_list = true;
						int j;
						for (j = topK-2; (j >= 0) && (Fcompare(final_score, res[j].score)==1); j--)
							res[j+1]=res[j];
						res[j+1].setR(smallest_did,final_score);
					} else {
						// find next smallest did
						smallest_did = CONSTS::MAXD;
						for (int i = lps.size()-1; i>=num_of_essential_lists; --i) {
							// pick the next smallest did
							if (lps[i]->did < smallest_did)
								smallest_did = lps[i]->did;
						}

						//togo
						//if (!fail)
							//PROFILER(CONSTS::GETBLOCKSCORE);
					}
			} else {
				//PROFILER(CONSTS::EARLYTERMINATION4);

				// temp value to help find the next smallest did
				int temp_next_smallest_did = CONSTS::MAXD;

				// skip safely in the essential set
				for (int i = lps.size()-1; i>=num_of_essential_lists; --i) {
					if (smallest_did == lps[i]->did) {
						//PROFILER(CONSTS::NEXTGEQ);
						lps[i]->did = lps[i]->nextGEQ( smallest_did + 1 );
					}

					// pick the next smallest did
					if (lps[i]->did < temp_next_smallest_did)
						temp_next_smallest_did = lps[i]->did;
				}

				// next smallest did has already been found
				smallest_did = temp_next_smallest_did;
			}
		} else {
			//PROFILER(CONSTS::EARLYTERMINATION2);
			skip_to_next_live_block(lps, threshold, num_of_essential_lists, check_up_to_list, 0, smallest_did, skip_to_next_live_block_enabled/*, Quantization_lookup_table*/);
		}
    } else {
    	//PROFILER(CONSTS::EARLYTERMINATION01);
    	skip_to_next_live_block(lps, threshold, num_of_essential_lists, num_of_essential_lists, 0, smallest_did, skip_to_next_live_block_enabled/*, Quantization_lookup_table*/);
    }
}

void DocidOriented_BMM_BMQ::find_new_non_essential_lists(unsigned int& num_essential_lists, unsigned int& num_of_non_essential_lists, std::vector <float>& prefix_sum_max_score, lptrArray& lps, const float& threshold, bool& essential_list_added) {
		while ((num_essential_lists<lps.size()-1) && (Fcompare(prefix_sum_max_score[num_essential_lists], threshold) == -1)) {   // loop until more than threshold
			++num_essential_lists;        // increase number of essential lists
			--num_of_non_essential_lists; // decrease number of non-essential lists
			essential_list_added = true;
		}
}

void DocidOriented_BMM_BMQ::skip_to_next_live_block(lptrArray& lps, const float& threshold, const int& num_essential_lists, const int& check_up_to_list, const int& find_next_live_block_from_list, int& smallest_did, const bool& skip_to_next_live_block_enabled/*, std::vector< std::vector<float> >& Quantization_lookup_table*//*,std::vector<float>& Quantization_lookup_table*/) {
	float sum_of_block_maxes = 0.0f;
	int temp_smallest_max_did;
	unsigned int smallest_max_did = lps[num_essential_lists]->gen.skipToNextBlock(smallest_did)-1;

	// check up to list is used because we may want to find the next live block in
	// a) essential set after failing test sum(ess. block maxes) + sum(non-ess. maxscores) or
	// b) essential set plus the lists from non-essential set, caused filter to skip block
	//    in other words, during non-essential block maxes refinement, we keep track of the
	//    list that caused to break and we need to check the block maxes from lps.size()-1 up to this list

	// find smallest max did in blocks in some lists
	for (int i = lps.size()-1; i>=check_up_to_list; --i) {
		temp_smallest_max_did = lps[i]->gen.skipToNextBlock(smallest_did)-1;
		if (smallest_max_did > temp_smallest_max_did)
			smallest_max_did = temp_smallest_max_did;
	}

   	// the did found is the min max did in blocks of checked lists, so we need to add 1 to get to the next block
	smallest_max_did += 1;

   	while (skip_to_next_live_block_enabled) {
   		if (smallest_max_did >= CONSTS::MAXD)
   			break;

   		// obtain block maxes for the did found
   		for (int i = lps.size()-1; i>=find_next_live_block_from_list; --i)
   			//sum_of_block_maxes += Quantization_lookup_table[i*(CONSTS::Quantization+1)+(unsigned int) lps[i]->gen.get_Quantized_Score(smallest_max_did)];
   			sum_of_block_maxes += lps[i]->quantile*(float) lps[i]->gen.get_Quantized_Score(smallest_max_did);
   			//sum_of_block_maxes += Quantization_lookup_table[i][(unsigned int) lps[i]->gen.get_Quantized_Score(smallest_max_did)];

   		// check if it can make it to the topk, if so exit (we have already computed next safe smallest max did)
   		if (Fcompare(sum_of_block_maxes, threshold) > -1)
   			break;
   		else {
   			//PROFILER(CONSTS::SORT);
   			// find next smallest max_did
   	    	unsigned int next_smallest_max_did = lps[num_essential_lists]->gen.skipToNextBlock(smallest_max_did)-1;
   	    	for (int i = lps.size()-1; i>=find_next_live_block_from_list; --i) {
   				temp_smallest_max_did = lps[i]->gen.skipToNextBlock(smallest_max_did)-1;
   				if (next_smallest_max_did > temp_smallest_max_did)
   					next_smallest_max_did = temp_smallest_max_did;
   			}
   	    	// update next smallest did
   	    	smallest_max_did = next_smallest_max_did + 1;

   	    	// corner case
  	    	if (smallest_max_did >= CONSTS::MAXD)
  	    		break;
   		}
   		// empty sum of block maxes
   		sum_of_block_maxes = 0.0f;
   	}

   	// if smallest_max_did in range find next smallest did
   	if (smallest_max_did < CONSTS::MAXD) {
   		// initialize the next smallest did with a large value
  		smallest_did = CONSTS::MAXD;

   		// perform skipping only on essential set
   		for (int i = lps.size()-1; i>=num_essential_lists; --i) {
  			//PROFILER(CONSTS::NEXTGEQ);
   			lps[i]->did = lps[i]->nextGEQ( smallest_max_did );

  			// pick the next smallest did
   			if (lps[i]->did < smallest_did)
   				smallest_did = lps[i]->did;
   		}
  	} else
  		smallest_did = CONSTS::MAXD;
}

void DocidOriented_BMM_BMQ::operator() (lptrArray& lps, const int& topK, QpResult* res, const float& threshold, const int& threshold_knowledge){
	unsigned int num_essential_lists = 0;
	unsigned int num_of_non_essential_lists = lps.size();
	float total_maxscore_sum = 0.0f;
	float current_threshold = threshold;
	float current_prefix_sum_max_score = 0.0f;
	const int threshold_unknown = (threshold_knowledge == 0) ? 1 : 0;//(Fcompare(threshold, 0.0f) == 0) ? 1 : 0; // treat as bool variable, use 1 when approximate score
	const int comparison = -1;
	const int comparison_heapify = (threshold_knowledge > 0) ? 0 : 1;
	bool check_for_new_essential_list = false;
	bool skip_to_next_live_block_enabled = true;	// set for BMM-NLB and unset for BMM
	bool essential_list_added = true;
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

	// initialize maxscore array
	for (int i=0; i<lps.size(); i++) {
		lists_maxscore.push_back(lps.getListMaxScore(i));
		lists_block_max_score.push_back(0.0f);
		total_maxscore_sum += lps.getListMaxScore(i);
	}

/*
	// Quantization lookup table - version 1. vector of vectors
	std::vector< std::vector<float> > Quantization_lookup_table (lps.size(), std::vector<float>((CONSTS::Quantization+1), 0.0f));

	for (int i=0; i<lps.size(); i++) {
		float temp_quantile = lps[i]->quantile;
		for (int j=0; j<=CONSTS::Quantization; j++)
			Quantization_lookup_table[i][j] = temp_quantile*(float) j ;
	}

	// Quantization lookup table - version 2. vector with offset
	std::vector<float> Quantization_lookup_table(lps.size()*(CONSTS::Quantization+1), 0.0f);

	for (int i=0; i<lps.size(); i++) {
		float temp_quantile = lps[i]->quantile;
		for (int j=0; j<=CONSTS::Quantization; j++)
			Quantization_lookup_table[i*(CONSTS::Quantization+1)+j] = temp_quantile*(float) j ;
	}
*/
	// initialize prefix sum max score array that is used in find_essential_lists
	prefix_sum_max_score.push_back(lps.getListMaxScore(0));
	for (int i=1; i<lps.size()-1; i++)
		prefix_sum_max_score.push_back(prefix_sum_max_score[i-1] + lps.getListMaxScore(i));

	// if threshold is known find essential lists
	if (threshold_knowledge > 0) { //(threshold_unknown==0) {
		find_new_non_essential_lists(num_essential_lists, num_of_non_essential_lists, prefix_sum_max_score, lps, current_threshold, essential_list_added);
		essential_list_added = false;
	}

	// find first smallest did
	find_smallest_did_in_essential_set(lps, num_essential_lists, smallest_did);

	while(true) {
		// if threshold is not known and heapify occured, update the current_threshold, find if a new list can be added to the essential set
		if ((/*threshold_unknown == 1*/threshold_knowledge<2) && (check_for_new_essential_list)) {
			// update threshold
			current_threshold = res[topK-1].score;
			// find essential lists
			find_new_non_essential_lists(num_essential_lists, num_of_non_essential_lists, prefix_sum_max_score, lps, current_threshold, essential_list_added);
			check_for_new_essential_list = false;

			// if essential list added, find new smallest did
			if (! essential_list_added) { // for one corner case
				// if new smallest did out of range break, otherwise return smallest_did in this round
				if (! find_smallest_did_in_essential_set(lps, num_essential_lists, smallest_did))
					break;
		    } else {
				essential_list_added = false;
				// if new smallest did out of range break, otherwise return smallest_did in this round
				if (! find_smallest_did_in_essential_set(lps, num_essential_lists, smallest_did))
					break;
			}
		} else {
			// if threshold is unknown, we have already found the smallest did from the previous round
			if (threshold_knowledge<1/*threshold_unknown == 1*/) {
				// check if next smallest did is out of docid range
				if (smallest_did >= CONSTS::MAXD)
					break;
			} else { // if threshold known AND no check_for_new_essential_list, then we need to find the next smallest did
				if (check_for_new_essential_list) {
					check_for_new_essential_list = false;

					// if smallest did out of range break, otherwise return smallest_did in this round
					if (! find_smallest_did_in_essential_set(lps, num_essential_lists, smallest_did))
						break;
				} else {
					// check if next smallest did is out of docid range
					if (smallest_did >= CONSTS::MAXD)
						break;
				}
			}
		}

		// set correct current prefix sum max score
		current_prefix_sum_max_score = num_essential_lists > 0 ? prefix_sum_max_score[num_essential_lists-1] : 0.0f;

		// evaluate
		evaluate(lps, res, current_threshold, topK, num_essential_lists, smallest_did, lists_block_max_score, current_prefix_sum_max_score, check_for_new_essential_list, /*threshold_unknown*/comparison_heapify, comparison, lists_maxscore, skip_to_next_live_block_enabled/*, Quantization_lookup_table*/);
	}
}

