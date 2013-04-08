/*
 * ListIterator.h
 *
 *  Created on: Jan 13, 2012
 *      Author: sergeyn
 */

#ifndef LISTITERATOR_H_
#define LISTITERATOR_H_

#include <vector>
#include <algorithm>
#include <iterator>

#include "SqlProxy.h"
#include "pfor.h"
#include "globals.h"
#include "profiling.h"
#include "sql/sqlite3.h"
#include "BlockGens.h"

class onDemandCpool {
	vecUInt data;
	std::string term;

	static size_t capacity;
	static size_t counter;
	static std::vector<onDemandCpool*> storage;
	static termsMap termsToInd;
	inline void setData(unsigned int* d, size_t sz) {
		data.clear();
		data.reserve(sz);
		for(size_t i=0; i<sz; ++i)
			data.push_back(d[i]);
	}

public:
	onDemandCpool()  {}
	~onDemandCpool() {evictData(); }

	void evictData();
	static void initPool(size_t cap);
	size_t errorCheck();
	inline void setName(const std::string& off) { term = off; }
	inline const std::string& getName() { return term ;}
	//the data must be materialized before using this operator!
	inline unsigned int& operator[](size_t i) { return data[i]; }
	inline size_t getSize() { return data.size(); }
	unsigned int materialize( const std::vector<std::string>& peers);
};

class BasicList {
public:
	std::string term;
	size_t termId;
	unsigned int lengthOfList;				// the padded list length
	unsigned int unpadded_list_length;		// needed only for layering (unpadded list length)
	unsigned int original_unpadded_list_length; // for layering, useless but we need to remove it from prepare_list - trecreader
	float maxScoreOfList; //optional
	BasicList(const std::string& t="", size_t id=0) : term(t), termId(id) {}
};

class RawIndexList : public BasicList{
public:
	vecUInt doc_ids;
	vecUInt freq_s;
	std::vector<float> scores; //optional
	explicit RawIndexList(const BasicList& rhs=BasicList()) :BasicList(rhs) { }
	void rankWithBM25(unsigned int* doclen); //will fill the scores
	void rankWithBM25(unsigned int* doclen, unsigned int& unpadded_list_length); // compute the correct max score of list
};

class CompressedList : public BasicList {
public:
	size_t BS;
	vecUInt maxDidPerBlock;		// the max value for each chunk
	std::vector<float> DEPmaxScorePerBlock;
	vecUInt sizePerBlock;		// the size for each chunk
	vecUInt flag;			// the flag for each chunk
	vecUInt compressedVector;	// the compressed data
	onDemandCpool cpool;
	size_t cpool_leng;	// the length of cpool
	size_t info_leng;		// the length of info

	CompressedList() :BS(0) {}
	CompressedList(const RawIndexList& rl, const size_t bs=CONSTS::BS);
	void compressDocsAndFreqs(const RawIndexList& rl);
	void serializeToFS(const std::string& rootPath);
	void serializeToDb(const std::string& rootPath, SqlProxy&);
};

class BaseLptr : public CompressedList{
public:
	int checkblock;
	int currMax;
	unsigned int *off;		// the running pointer
	unsigned int *decmpPtr;		// another running pointer
	unsigned int block; 	// the current block number
	unsigned int elem;
	int did;
	int dflag;
	int fflag;
	unsigned int dbuf[CONSTS::BS];		// buffer for storing the decoded docid
	unsigned int fbuf[CONSTS::BS];		// buffer for storing the decoded freq
	float pre;

  	// Quantization
	float quantile;

///*
	// Layering
	bool is_essential;		  			// if term is essential
	bool has_layers;		  			// if term has layers true
	unsigned int list_id;				// if term has layers, then list_id
	int layer_index;					// their layered id, otherwise
	int layer_status;
	//int metainfo;						// metainfo needed during early termination

	inline void reset() { block = 0; dflag = 1; fflag = 0; elem = 0; checkblock = 0; }

	inline BaseLptr() { reset(); }
	void open(const std::vector<std::string>& peers, const unsigned int lengthOfLst);

	// Usage: reset all pointers so that we can re-use list
	inline void reset_list() {
		block = 0;
		dflag = 1;
		fflag = 0;
		elem = 0;
		checkblock = 0;
		off = &(cpool[0]);
		currMax = maxDidPerBlock[0];
		decmpPtr = pack_decode(dbuf, off, flag[block]&65535, 0);
		did = dbuf[0];
	}
};

// structure for a pointer within an inverted list
class lptr : public BaseLptr {
public:
	class comparator {
	public:
		inline bool operator()(lptr* a, lptr* b) { return a->did<b->did;} //to sort by list length -- aka comparelps
	};

	class scorecomparator { // ADDON
		public:
			inline bool operator()(lptr* a, lptr* b) { return a->maxScoreOfList<b->maxScoreOfList;} // sort by score (function for sort())
		};

	#define CALC_SCORE(f,did) pre * (float)(f) / ((float)(f + 2 * ( 0.25 + 0.75 * (float)(did)/(float)(130) ) ) )
	SingleHashBlocker gen;
	unsigned int maskOn;
	unsigned int maskOff;
	unsigned int indexInSorted;
	lptr() {}

	inline float calcScore(float f, int did) { return CALC_SCORE(f,did);}

	inline float DEPgetMaxScoreOfDeepBlock() { return DEPmaxScorePerBlock[block]; }

	// used for the MaxscoreBlocks
	inline float getMaxScoreOfBlock() {	return DEPmaxScorePerBlock[checkblock]; }

	//move the pointer to the right block
	inline int DEPshallowMove(int did) MYNOINLINE {
		while (did > (int)maxDidPerBlock[checkblock]) 		//++total_move;
			++checkblock;
		return maxDidPerBlock[checkblock] + 1 ;
	}

	// get the frequency of the current posting
	inline int getFreq() MYNOINLINE	{
		if (fflag < 1)  {
			decmpPtr = pack_decode(fbuf,decmpPtr, flag[block]>>16);
			fflag = 1;
		}
		return(fbuf[elem]+1);
	}

	// skip blocks until you come to the right one
	inline void skipToDidBlockAndDecode(int d) MYNOINLINE	{
		//profiler::getInstance().stepCounter(CONSTS::SKIPS);
		off += (sizePerBlock[block++]+1);
		while (d > (int)maxDidPerBlock[block]) {
			off += (sizePerBlock[block++]+1);
			//profiler::getInstance().stepCounter(CONSTS::SKIPSB);
		}

		currMax=maxDidPerBlock[block]; //we keep track of current max to avoid max[block] in nextGEQ -- it is not cache friendly...
		decmpPtr = pack_decode(dbuf, off, flag[block]&65535, maxDidPerBlock[block-1]); //decode new block
		elem = 0; //point to first element in that block
		did = dbuf[0]; //get the first did IN THE DECODED block
		dflag = 1;
		fflag = 0;

	}

	// move list pointer to first posting with docID >= d */ //
	inline unsigned int nextGEQ(int d) MYNOINLINE	{
		PROFILER(CONSTS::NEXTGEQ);
		if (d > currMax) // check if we need to skip to a new block of size BS
			skipToDidBlockAndDecode(d);

		for(elem=elem+1; d>did; ++elem) //try to go to the right element
			did +=dbuf[elem]; //prefix summing
		--elem;

		return did;
	}
};


class lptrArray : public std::vector<lptr*> {
public:
	lptrArray() {}
	inline void sort() {
		std::sort(begin(),end(),lptr::comparator()); //TODO: speedup
		for(unsigned int i=0; i<size(); ++i)
			this->operator[](i)->indexInSorted = i;
	}

	inline void sort(const int start, const int end) { // TODO handle error if put pivot out of range
		std::sort(begin() + start, begin() + end + 1, lptr::comparator());
		for (unsigned int i=start; i<end+1; ++i)
			this->operator[](i)->indexInSorted = i;
	}

	// ADDON
	inline void sortbyscore() {
		std::sort(begin(),end(),lptr::scorecomparator());
		for(unsigned int i=0; i<size(); ++i)
			this->operator[](i)->indexInSorted = i;
	}

	inline void sortbyscore(const int pivot) { // TODO handle error if put pivot out of range
		std::sort(begin(),begin() + pivot + 1,lptr::scorecomparator());
		for(unsigned int i=0; i<pivot + 1; ++i)
			this->operator[](i)->indexInSorted = i;
	}

	inline float getListMaxScore(size_t i) { return this->operator [](i)->maxScoreOfList; }

	inline void popdown(int i) 	{
		lptr *temp = this->operator [](i);
		int j = i + 1;
		while( j < this->size() && temp->did >= ( (this->operator [](j))->did) )	{
			this->operator [](j-1) = this->operator [](j); //swap
			this->operator [](j)->indexInSorted = j-1;
			++j;
		}
		this->operator [](j - 1) = temp;
		this->operator [](j-1)->indexInSorted = j-1;
	}
};

class termsCache {
	std::vector<lptr> data;

public:
	termsCache() {}
	void addSingleToCache(const std::string& term);
	void fill_cache(const char* mappingPath);
	inline lptr& operator[] (size_t i) { return data[i]; }
	inline size_t size() { return data.size(); }
};

template<typename T>
void injectBlocker(T& blocker, RawIndexList& rilist ) {
	T gen(rilist.doc_ids);
	gen.generateBlocks();
	gen.generateMaxes(rilist.scores);
	blocker = gen;
}

//Quantization
template<typename T>
void Quantized_injectBlocker(T& blocker, RawIndexList& rilist , float& quantile) {
	T gen(rilist.doc_ids);
	gen.generateBlocks();
	gen.generate_Quantized_Maxes(rilist.scores, quantile);
	blocker = gen;
}

template<typename T>
void clearBlocks(T& blocker) {
	blocker = T();
}

#endif /* LISTITERATOR_H_ */
