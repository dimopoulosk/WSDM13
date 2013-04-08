/*
 * exhaustiveOR.h
 *
 *  Created on: May 17, 2012
 *      Author: constantinos
 */

#ifndef EXHAUSTIVEOR_H_
#define EXHAUSTIVEOR_H_

#include "PostingOriented_BMW.h"

class ExhaustiveOR{
private:
	unsigned int* pages;

public:
	ExhaustiveOR(unsigned int* pgs) : pages(pgs) {}
	void operator()(lptrArray& lps, const int topk, QpResult* r);
};

#endif /* EXHAUSTIVEOR_H_ */
