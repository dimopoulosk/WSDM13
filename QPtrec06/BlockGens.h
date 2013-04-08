/*
 * BlockGens.h
 *
 *  Created on: Feb 27, 2012
 *      Author: sergeyn
 */

#ifndef BLOCKGENS_H_
#define BLOCKGENS_H_

#include "globals.h"
#include "utils.h"
#include <cmath>

typedef vecUInt SingleBlock;
typedef std::vector<SingleBlock>::const_iterator BlocksIterator;

template<typename T>
T getMaxIndirect(const SingleBlock& blockInds, const std::vector<T>& vals) {
	T m = vals[ *(blockInds.cbegin())];
	for( auto it = blockInds.cbegin(); it != blockInds.cend(); ++it)
		if(m < vals[(*it)] )
			m = vals[(*it)];
	return m;
}

template<typename T>
std::vector<T> getVaues(const SingleBlock& blockInds, const std::vector<T>& vals) {
	std::vector<T> tmp;
	tmp.reserve(blockInds.size());
	for(SingleBlock::const_iterator it = blockInds.cbegin(); it != blockInds.end(); ++it)
		tmp.push_back(vals[*it]);
	return tmp;
}

class AbsBlocksGen : public std::vector<SingleBlock> {
public:
	AbsBlocksGen() {}
	virtual ~AbsBlocksGen() {}
	BlocksIterator getIterator() { return this->cbegin(); }
	BlocksIterator getEnd() { return this->cend(); }
	virtual void generateBlocks() {}

	static unsigned int getMaxDidOfBlock(const SingleBlock& block, const vecUInt& dids ) { return getMaxIndirect(block,dids); }
	static float getMaxScoreOfBlock(const SingleBlock& block, const std::vector<float>& scores) { return getMaxIndirect(block,scores); }
	static unsigned int getMaxDidOfBlock(const BlocksIterator& block, const vecUInt& dids ) { return getMaxIndirect(*block,dids); }
	static float getMaxScoreOfBlock(const BlocksIterator& block, const std::vector<float>& scores) { return getMaxIndirect(*block,scores); }
};

class PostingsBlocker : public AbsBlocksGen {
	size_t BS;
	size_t length;

public:
	PostingsBlocker(size_t len, size_t bs = CONSTS::BS ) : BS(bs), length(len){}
	virtual void generateBlocks();
};

class PostingsBlockerWithScores : public PostingsBlocker {
	const vecUInt* docIds;
	hashMapStd<unsigned int, unsigned int> didToBlock;
	std::vector<float> blockIdToMaxScore;

	inline unsigned int whereDidMapped(unsigned int did) { return didToBlock[did];}
public:
	PostingsBlockerWithScores() :PostingsBlocker(0,0),docIds(0) {}
	PostingsBlockerWithScores(vecUInt& dids, size_t bs = CONSTS::BS ) : PostingsBlocker(dids.size(),bs), docIds(&dids){}
	void generateMaxes(const std::vector<float>& scores);
	inline float getMaxScore(unsigned int did) { return blockIdToMaxScore[whereDidMapped(did)]; }
	inline unsigned int skipToNextBlock(unsigned int did) {
		unsigned int block = whereDidMapped(did);
		while(whereDidMapped(++did) == block) ;
		return did;
	}
};

class FixedBitsHashBlocker : public AbsBlocksGen {
protected:
	unsigned char bits;
	const vecUInt* docIds;
	std::vector<float> blockIdToMaxScore;
	//std::vector<unsigned char> QuantizedblockIdToMaxScore;
	inline unsigned int whereDidMapped(unsigned int did) { return did >> bits;}
public:
	static unsigned int expectedBlockSize;
	std::vector<unsigned char> QuantizedblockIdToMaxScore; // in protected ?
	FixedBitsHashBlocker(const vecUInt& dids, unsigned char _bits, size_t sz=0) : bits(_bits ? _bits : expectedBlockSize), docIds(&dids) {}

	unsigned int getBits() const { return bits; }
	const std::vector<float>& getScores() const { return blockIdToMaxScore; }

	void generateMaxes(const std::vector<float>& scores);
	//void generate_Quantized_Maxes(const std::vector<float>& scores, float& quantile); // Quantization
	//void generate_Float_Quantized_Maxes(const std::vector<float>& scores, float& quantile); // Float Quantization
	inline virtual float getMaxScore(unsigned int did) { return blockIdToMaxScore[whereDidMapped(did)]; }
	//inline virtual float get_Quantized_Score(unsigned int did) { return QuantizedblockIdToMaxScore[whereDidMapped(did)]; }
	//inline virtual float get_Quantized_Block_Score(int block_number) { return QuantizedblockIdToMaxScore[block_number]; }
	inline virtual unsigned int skipToNextBlock(unsigned int did) { return (whereDidMapped(did)+1)<<bits; }
};

class SingleHashBlocker : public AbsBlocksGen {
protected:
	unsigned char bits;
	const vecUInt* docIds;
	std::vector<float> blockIdToMaxScore;
	//std::vector<unsigned char> QuantizedblockIdToMaxScore;
	inline unsigned int whereDidMapped(unsigned int did) { return did >> bits;}
public:
	static unsigned int expectedBlockSize;
	std::vector<unsigned char> QuantizedblockIdToMaxScore; // in protected ?
	SingleHashBlocker(size_t sz=0) : bits(10), docIds(0) {}

	static unsigned int getLogBlockSize(int length) {
		// Pick either (a) or (b) and comment accordingly
		// (a) Expected Block Selection scheme - In other words, expected postings per docID oriented block
		// First, set the number of expected postings you want per block in BITS, for example if you want 4 expected postings per block, set expectedBits = 2 (because 2^2=4)
/*		unsigned int expectedBits = expectedBlockSize;
		unsigned int maxdBits = CONSTS::MAXDBITS;
		unsigned int lenBits = intlog2(length);
		unsigned int bbits = (maxdBits/lenBits)*expectedBits; // assuming uniform distribution of docIDs, maxdBits/lenBits gives the size of the block in order to have 1 expected posting per block, so multiply by expectedBits to get the right # of expected postins per block
*/
		// (b) Fixed Block Selection scheme - fixed block size for all terms
		// Just change the return value to the number of fixed bits you want for all terms, i.e. return 6 (means 2^6 = 64 consecutive docIDs per block and so on ...)
		unsigned int bbits = 7;
		return bbits;
	}

	static unsigned int bitOracle(unsigned int length) {
		unsigned int lenBits = intlog2(length);
  		std::vector<int> Bucket (18, 0);
  		// Parameters
		Bucket.at(0) = 10; // 2^7 (list length less than 2^7) ...
		Bucket.at(1) = 10; // 2^8
		Bucket.at(2) = 10; // 2^9
		Bucket.at(3) = 10; // 2^10
		Bucket.at(4) = 6; // 2^11
		Bucket.at(5) = 6; // 2^12
		Bucket.at(6) = 7; // 2^13
		Bucket.at(7) = 7; // 2^14
		Bucket.at(8) = 7; // 2^15
		Bucket.at(9) = 8; // 2^16
		Bucket.at(10) = 8; // 2^17
		Bucket.at(11) = 7; // 2^18
		Bucket.at(12) = 6; // 2^19
		Bucket.at(13) = 6; // 2^20
		Bucket.at(14) = 6; // 2^21
		Bucket.at(15) = 6; // 2^22
		Bucket.at(16) = 6; // 2^23
		Bucket.at(17) = 6; // 2^24

  		// pick the right bucket
		int B_counter = 0;
		for (int i=7; i<25; i++) {
			if (lenBits<i)
				return Bucket.at(B_counter);
			else
				++B_counter;
		}

/*
// maxscore DEFAULT BEST version
		if (lenBits<8)
			return 10;
		else
			return 6;
*/

/* hashbmw version
		if (lenBits<11)
			return 10;
		else //if (lenBits<20)
			return 7;
		//else
		//	return 7;
*/
	}

	// BlockGenCost: given lenbits and bbits give the extra bits hash needed as cost
	static unsigned int BlockGenCost(unsigned int lenBits, unsigned int bitshift) {
		unsigned int PFDelta_num_blocks = lenBits - 6; // dont hardcode 6 => 2^6 = 64 // we are sure lenBits always more than 6
		unsigned int Hash_num_blocks = lenBits - bitshift;
		return (Hash_num_blocks - PFDelta_num_blocks);
	}

	unsigned int getBits() const { return bits; }
	const std::vector<float>& getScores() const { return blockIdToMaxScore; }
	SingleHashBlocker(vecUInt& dids) : bits(bitOracle(dids.size())), docIds(&dids) { // CHANGE bitOracle()to getLogBlockSize(dids.size()) for fixed and expected block selection scheme AND SEE THE code of getLogBlockSize() and change the code accordingly
	}

	void generateBlocks() {}
	void generateMaxes(const std::vector<float>& scores);
	void generate_Quantized_Maxes(const std::vector<float>& scores, float& quantile); // Quantization
	void generate_Float_Quantized_Maxes(const std::vector<float>& scores, float& quantile); // Float Quantization
	inline virtual float getMaxScore(unsigned int did) { return blockIdToMaxScore[whereDidMapped(did)]; }
	inline virtual float get_Quantized_Score(unsigned int did) { return QuantizedblockIdToMaxScore[whereDidMapped(did)]; }
	inline virtual float get_Quantized_Block_Score(int block_number) { return QuantizedblockIdToMaxScore[block_number]; }
	inline virtual unsigned int skipToNextBlock(unsigned int did) { return (whereDidMapped(did)+1)<<bits; }
};

class OnTheFlyBlocker : public SingleHashBlocker {
public:
	explicit OnTheFlyBlocker(const SingleHashBlocker& largerBlocker, const unsigned int shrinkFactor) {
		const std::vector<float>& scores(largerBlocker.getScores());
		const size_t sz(scores.size());
		blockIdToMaxScore.reserve(sz/2);
		for(size_t i=0; i<sz; i+=shrinkFactor)
			blockIdToMaxScore.push_back(*std::max_element(scores.begin()+i,scores.begin()+i+shrinkFactor));
	}
};

#endif /* BLOCKGENS_H_ */

