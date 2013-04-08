/*
 * loadIndex.cpp
 *
 *  d on: Jan 12, 2012
 *      Author: sergeyn
 */

#include "ListIterator.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sql/sqlite3.h"
#include "math.h"
#include "BlockGens.h"

using namespace std;

float temp_score[CONSTS::MAXD + 128];
unsigned int temp_store[CONSTS::MAXD];

size_t onDemandCpool::counter;
size_t onDemandCpool::capacity;
std::vector<onDemandCpool*> onDemandCpool::storage;
termsMap onDemandCpool::termsToInd;

void onDemandCpool::initPool(size_t cap) {
	counter = 0;
	capacity = cap;

	storage.reserve(cap);
	for(size_t i=0; i<cap; ++i)
		storage.push_back(NULL);
}

unsigned int readEntireFileToTmpBuffer(const char* fname, unsigned int* buffer){
	FILE* fin = fopen(fname,"r");
	if(fin == NULL) {
		CERR << "couldn't open " << fname << Log::endl;
		return 0;
	}

	unsigned int r1 = 0;
	while( fread(&(buffer[r1]), sizeof(unsigned int), 1, fin) == 1 )
		r1 ++;

	fclose(fin);
	return r1;
}

inline bool isCollision(const std::string& victim, const std::vector<std::string>& peers) {
	for(size_t i=0; i<peers.size(); ++i)
		if(victim == peers[i])
			return true;
	return false;
}

void onDemandCpool::evictData() {
	data.clear();
	data = vecUInt(); //clearing and resizing was causing mem. leaks on my system - hence this
}

unsigned int onDemandCpool::materialize( const std::vector<std::string>& peers) {
	//test if hit
	termsMap::const_iterator ind = termsToInd.find(term);
	if(ind != termsToInd.end() && (*ind).second >= 0)
		return data.size(); //this one is materialized already

    //no hit, so we need to place it in storage
	while(storage[counter] && isCollision(storage[counter]->term,peers)) //pick non-colliding victim for eviction
		++counter;

	if(storage[counter]) { //the bucket is occupied
		storage[counter]->evictData();
		termsToInd[storage[counter]->term]=-1; //evict
		//profiler::getInstance().stepCounter(CONSTS::CNT_EVICTIONS);
	} // when query has duplicate terms this scheme should still work

	termsToInd[term]=counter;
	storage[counter]=this;
	counter = (counter+1)%capacity;

	//now actually get the data, update the buffer pointer and the size
	size_t size = readEntireFileToTmpBuffer(term.c_str(),temp_store);

	data.reserve(size);
	for(unsigned int i=0;i<size;i++)
			data.push_back(temp_store[i]);

	return size;
}

template <typename T>
void getBlob(sqlite3_stmt* ppStmt, int cid, std::vector<T>& v) {
	if(sqlite3_column_type(ppStmt, cid) != SQLITE_BLOB)
		return;
	const int size = sqlite3_column_bytes(ppStmt, cid) / sizeof(T);
    const T* src = (const T*) sqlite3_column_blob(ppStmt, cid);
    v.clear();
    v.reserve(size);
    for(int i=0; i<size; ++i)
    	v.push_back(src[i]);
}

std::vector<lptr> load(const std::string& where=""){
	SqlProxy ld;
	sqlite3_stmt* ppStmt = ld.prepare("select * from terms "+where); //tODO: make iterator
	if(! ppStmt)
		CERR << "select failed" << EFATAL;


  std::vector<lptr> centry;
  centry.reserve(11000);

  //'/flag/','/max/','/min/','/score/','/size/'
  while (sqlite3_step(ppStmt) == SQLITE_ROW) {  // For each row returned
	  lptr c;
	  c.term = std::string((const char*)sqlite3_column_text(ppStmt, 0));
	  c.cpool.setName(CONSTS::TERM_PATH + "/" + c.term);
	  c.info_leng = sqlite3_column_bytes(ppStmt, 1) / sizeof(unsigned int);
	  getBlob(ppStmt,1,c.flag);
	 // c.maxDidPerBlock = (unsigned int*) getBlob(ppStmt,2);
	  getBlob(ppStmt,2,c.maxDidPerBlock);
	  getBlob(ppStmt,4,c.DEPmaxScorePerBlock);
	  getBlob(ppStmt,5,c.sizePerBlock);
	  c.maxScoreOfList = *(std::max_element(c.DEPmaxScorePerBlock.begin(),c.DEPmaxScorePerBlock.end()));
	  centry.push_back(c);
  } //end while
  sqlite3_finalize(ppStmt);
  return centry;
}

void termsCache::addSingleToCache(const std::string& term) {
	std::vector<lptr> arr = load("where term='"+term+"'");
	data.push_back(arr[0]);
}

void termsCache::fill_cache(const char* mappingPath){
	onDemandCpool::initPool(CONSTS::STORAGE_CAPACITY);
	data = load();
	COUT2<< "Loaded: " << data.size() << " index info items" << Log::endl;
}

void RawIndexList::rankWithBM25(unsigned int* doclen) {
	float idf;
	scores.reserve(lengthOfList);
	for(size_t i = 0; i<lengthOfList; ++i)	{
		unsigned int did = doc_ids[i];
		float f  = freq_s[i];
		unsigned int dlen = doclen[ did ];
		idf = log( 1+(float)(CONSTS::MAXD - lengthOfList + 0.5)/(float)(lengthOfList + 0.5))/log(2) ;
		scores.push_back(idf * ( (float)(f * 3)/(float)(f + 2*( 0.25 + 0.75 * (float)(dlen)/(float)(130) ) ) ));
	}
	maxScoreOfList = *(std::max_element(scores.begin(),scores.end()));
}

void RawIndexList::rankWithBM25(unsigned int* doclen, unsigned int& unpadded_list_length) {
	float idf;
	scores.reserve(unpadded_list_length);
	for(size_t i = 0; i<unpadded_list_length; ++i)	{
		unsigned int did = doc_ids[i];
		float f  = freq_s[i];
		unsigned int dlen = doclen[ did ];
		idf = log( 1+(float)(CONSTS::MAXD - unpadded_list_length + 0.5)/(float)(unpadded_list_length + 0.5))/log(2) ;
		scores.push_back(idf * ( (float)(f * 3)/(float)(f + 2*( 0.25 + 0.75 * (float)(dlen)/(float)(130) ) ) ));
	}
	maxScoreOfList = *(std::max_element(scores.begin(),scores.begin() + unpadded_list_length));
}

template <typename T>
void dumpArrayToFile(const std::string& path, T* data, size_t sz) {
	FILE *fdo = fopen(path.c_str(),"w");
	for (size_t i = 0; i < sz; ++i)
		if (fwrite(&(data[i]), sizeof(T), 1, fdo) != 1)
			CERR << "Failed writing to " << path << Log::endl;
	fclose(fdo);
}

template <typename Iterator>
void dumpToFile(const std::string& path, Iterator start, Iterator end) {
	FILE *fdo = fopen(path.c_str(),"w");
	if(! fdo) {
		CERR << "Failed writing to " << path << " (Hint: create folder for path?)" << Log::endl;
		return;
	}
	assert(fdo);
	for (; start !=end; ++start)
		if (fwrite(&(*start), sizeof(typename Iterator::value_type), 1, fdo) != 1)
			CERR << "Failed writing to " << path << " (Hint: create folder for path?)" << Log::endl;
	fclose(fdo);
}

void CompressedList::serializeToFS(const std::string& rootPath ) {
	const std::string TERM_PATH = rootPath+"pool/";
	const std::string FLAG_PATH =  rootPath+"flag/";
	const std::string MAX_PATH = rootPath+"max/";
	const std::string SCORE_PATH = rootPath+"score/";
	const std::string SIZE_PATH = rootPath+"size/";

	dumpToFile(TERM_PATH+term,compressedVector.begin(),compressedVector.end()); // dump the compressed list
	dumpToFile(FLAG_PATH+term,flag.begin(),flag.end()); // dump the flag

	dumpToFile(SIZE_PATH+term,sizePerBlock.begin(),sizePerBlock.end()); // dump the size for each chunk
	dumpToFile(MAX_PATH+term,maxDidPerBlock.begin(),maxDidPerBlock.end());
	dumpToFile(SCORE_PATH+term,DEPmaxScorePerBlock.begin(), DEPmaxScorePerBlock.end());// dump the max score for each chunk
	// dump the min did for each chunk -- nobody used it
}

void CompressedList::serializeToDb(const std::string& rootPath, SqlProxy& sql) {
	  const std::string TERM_PATH = rootPath+"pool/";

//	  const char *zSql = "insert or ignore into terms values (?,?,?,?,?,?)";
	  sqlite3_stmt *pStmt =  sql.prepare("insert or replace into terms values (?,?,?,?,?,?)");
	  const size_t numb = lengthOfList/BS;

/*
	    rc = sqlite3_prepare(sql.getDB(), zSql, -1, &pStmt, 0); // Compile the INSERT statement into a virtual machine.
	    while(rc == SQLITE_BUSY || rc == SQLITE_LOCKED) { //when running in a MT env. tables could be locked
	    	    	sleep(1);
	    	    	rc = sqlite3_prepare(sql.getDB(), zSql, -1, &pStmt, 0);
   	    }
*/
	    sqlite3_bind_text(pStmt, 1, term.c_str(), -1, SQLITE_STATIC);
	    sqlite3_bind_blob(pStmt, 2, &(flag[0]), numb*sizeof(unsigned int), SQLITE_STATIC);
	    sqlite3_bind_blob(pStmt, 3, &(maxDidPerBlock[0]), numb*sizeof(unsigned int), SQLITE_STATIC);
	    sqlite3_bind_blob(pStmt, 4, &(maxDidPerBlock[0]), 1, SQLITE_STATIC); //the min dummy
	    sqlite3_bind_blob(pStmt, 5, &(DEPmaxScorePerBlock[0]), numb*sizeof(float), SQLITE_STATIC);
	    sqlite3_bind_blob(pStmt, 6, &(sizePerBlock[0]), numb*sizeof(unsigned int), SQLITE_STATIC);
		COUT1 << "dump to file: " << term <<  " " <<  compressedVector.size() << " " << numb << " comp. ints" << Log::endl;
	    dumpToFile(TERM_PATH+term,compressedVector.begin(),compressedVector.end()); // dump the compressed list
	    sql.executeStoredStatement();

/*
	    rc = sqlite3_step(pStmt);
	    while(rc == SQLITE_BUSY || rc == SQLITE_LOCKED) { //when running in a MT env. tables could be locked
	    	if(rc == SQLITE_LOCKED)
	    		sqlite3_reset(pStmt);
	    	sleep(1);
	    	rc = sqlite3_step(pStmt);
	    }

	    assert( rc!=SQLITE_ROW );
	    rc = sqlite3_finalize(pStmt);
	  } while( rc==SQLITE_SCHEMA );
	  assert(rc == SQLITE_OK);
	  sqlite3_finalize(pStmt);
	  */
}

void CompressedList::compressDocsAndFreqs(const RawIndexList& rl) {
	size_t ind;

	DEPmaxScorePerBlock.reserve(lengthOfList/BS);
	maxDidPerBlock.reserve(lengthOfList/BS);

    //build delta list for doc ids
	vecUInt deltaList = rl.doc_ids;
	for (ind = lengthOfList-1; ind > 0; ind--)
		deltaList[ind] -= deltaList[ind-1];

	//build delta list for freqs
	vecUInt freq = rl.freq_s;
	for (ind = 0; ind < lengthOfList; ++ind)
		freq[ind]--;

	assert(lengthOfList == deltaList.size() && lengthOfList == freq.size());

	int j,flag1, flag2, numb;
	unsigned int *start, *ptrOfDelta, *ptrOfFreq;

	vecUInt buffer(lengthOfList*2);
	buffer.resize(lengthOfList*2);

	flag.clear();
	flag.reserve(lengthOfList/BS);
	sizePerBlock.clear();
	sizePerBlock.reserve(lengthOfList/BS);

	start = &buffer[0];
	PostingsBlocker blocker(deltaList.size(),BS);
	blocker.generateBlocks();

	//blocks iterator:
	//for (size_t i = 0; i < lengthOfList; i += BS)	{
	for(BlocksIterator it = blocker.getIterator(); it!=blocker.getEnd(); ++it) {
		vecUInt deltaTmp = getVaues(*it,deltaList);
		vecUInt freqTmp = getVaues(*it,freq);
		std::vector<float> scoresTmp = getVaues(*it,rl.scores);
		unsigned int* deltaPtr = &(deltaTmp[0]);
		unsigned int* freqPtr = &(freqTmp[0]);

		maxDidPerBlock.push_back(AbsBlocksGen::getMaxDidOfBlock(it,rl.doc_ids)); //(rl.doc_ids[i+BS-1]);
		DEPmaxScorePerBlock.push_back(AbsBlocksGen::getMaxScoreOfBlock(it,rl.scores)); //( *(std::max_element(rl.scores.begin()+i,rl.scores.begin()+i+BS)) );

		flag1 = -1;
		for (j = 0; flag1 < 0; j++)	{
			ptrOfDelta = start;
			flag1 = pack_encode(&ptrOfDelta,deltaPtr, j);
		}

		flag2 = -1;
		for (j = 0; flag2 < 0; j++)	{
			ptrOfFreq = ptrOfDelta;
			flag2 = pack_encode(&ptrOfFreq, freqPtr, j);
		}

		flag.push_back(((unsigned int)(flag1))|(((unsigned int)(flag2))<<16));

		if (ptrOfFreq - start > 256)
			CERR << "One block is " << ptrOfFreq - start << " bytes" << EFATAL;

		sizePerBlock.push_back( (unsigned char)((ptrOfFreq - start) - 1));
		numb++;
		start = ptrOfFreq;
	} //end for block

	cpool_leng = start-&(buffer[0]);
	const size_t cells = cpool_leng;//sizeof(unsigned int);
	compressedVector.reserve(cells);
	for(size_t i=0; i<cells; ++i)
		compressedVector.push_back(buffer[i]);

	COUT1 << term << "lengthOfList: " << lengthOfList << " deltaList.size: " << deltaList.size() << "cpool_leng: " << cpool_leng << Log::endl;

}

CompressedList::CompressedList(const RawIndexList& rl, const size_t bs) : BS(bs) {
	lengthOfList = rl.lengthOfList;
	maxScoreOfList = rl.maxScoreOfList;
	term = rl.term;
	termId = rl.termId;

	compressDocsAndFreqs(rl);
}

void BaseLptr::open(const std::vector<std::string>& peers, const unsigned int lengthOfLst) {
	reset();
	//materialize compressed data if not present and save its size
	cpool_leng = cpool.materialize(peers);
	off = &(cpool[0]);
	currMax = maxDidPerBlock[0];
	decmpPtr = pack_decode(this->dbuf, this->off, this->flag[this->block]&65535, 0);
	did = dbuf[0];
	unpadded_list_length = lengthOfLst;
	// compute pre using the unpadded list length
	pre = (float)3 * (float)(  log( 1 + (float)(CONSTS::MAXD - unpadded_list_length + 0.5)/(float)( unpadded_list_length + 0.5))/log(2));

	// set padded list length without storing anything
	// Note: lengthofList is unsigned int, but the expression can not be less than 0, so it's safe
	lengthOfList = lengthOfLst + (CONSTS::BS - (unpadded_list_length%CONSTS::BS));
}

void PostingsBlocker::generateBlocks() {
	unsigned int i = 0;
	assert(0==length%BS); //we currently assume that the postings count is block aligned
	vecUInt indices;
	for(unsigned int b = 0; b < length; ++b ) {
		indices.push_back(b);
		i = (i+1)%BS;
		if(i==0) {
			this->push_back(indices);
			indices.clear();
		}
	}
}

//#if __GNUC__ == 4 && __GNUC_MINOR__ > 4
unsigned int SingleHashBlocker::expectedBlockSize(8);
unsigned int FixedBitsHashBlocker::expectedBlockSize(8);

void FixedBitsHashBlocker::generateMaxes(const std::vector<float>& scores) {
	size_t block = 0;
	std::vector<float> localScores;
	const size_t numOfBlocks(1+whereDidMapped(CONSTS::MAXD));

	assert(scores.size() == docIds->size());

	for(size_t i = 0; i<numOfBlocks;  ++i)
		blockIdToMaxScore.push_back(0.0);

//	if(CONSTS::MAXD < docIds->at( docIds->size()-1) )
	//	COUT << docIds->at(docIds->size()-1) << Log::endl;

	for(size_t i=0; i<docIds->size() && docIds->at(i)<CONSTS::MAXD; ++i) {
		const unsigned int cblock =  whereDidMapped((*docIds)[i]);
		if(cblock != block) {
			if(localScores.size())
				blockIdToMaxScore[block] = *(std::max_element(localScores.cbegin(), localScores.cend()));
			localScores.clear();
			block = cblock;
		}
		localScores.push_back(scores[i]);
	}

	if(block == numOfBlocks)
		CERR << "ha?" << Log::endl;
	if(localScores.size())
		blockIdToMaxScore[block] = *(std::max_element(localScores.cbegin(), localScores.cend()));
}

void SingleHashBlocker::generateMaxes(const std::vector<float>& scores) {
	size_t block = 0;
	std::vector<float> localScores;
	const size_t numOfBlocks(1+whereDidMapped(CONSTS::MAXD));

	assert(scores.size() == docIds->size());

	for(size_t i = 0; i<numOfBlocks;  ++i)
		blockIdToMaxScore.push_back(0.0);

//	if(CONSTS::MAXD < docIds->at( docIds->size()-1) )
	//	COUT << docIds->at(docIds->size()-1) << Log::endl;

	for(size_t i=0; i<docIds->size() && docIds->at(i)<CONSTS::MAXD; ++i) {
		const unsigned int cblock =  whereDidMapped((*docIds)[i]);
		if(cblock != block) {
			if(localScores.size())
				blockIdToMaxScore[block] = *(std::max_element(localScores.cbegin(), localScores.cend()));
			localScores.clear();
			block = cblock;
		}
		localScores.push_back(scores[i]);
	}

	if(block == numOfBlocks)
		CERR << "ha?" << Log::endl;
	if(localScores.size())
		blockIdToMaxScore[block] = *(std::max_element(localScores.cbegin(), localScores.cend()));
}

void SingleHashBlocker::generate_Float_Quantized_Maxes(const std::vector<float>& scores, float& quantile) {
	size_t block = 0;
	std::vector<float> localScores;
	const size_t numOfBlocks(1+whereDidMapped(CONSTS::MAXD));

	assert(scores.size() == docIds->size());

	for(size_t i = 0; i<numOfBlocks;  ++i)
		blockIdToMaxScore.push_back(0.0);

//	if(CONSTS::MAXD < docIds->at( docIds->size()-1) )
	//	COUT << docIds->at(docIds->size()-1) << Log::endl;

	for(size_t i=0; i<docIds->size() && docIds->at(i)<CONSTS::MAXD; ++i) {
		const unsigned int cblock =  whereDidMapped((*docIds)[i]);
		if(cblock != block) {
			if(localScores.size()) {
				// get float block max
				float blockmax_float = *(std::max_element(localScores.cbegin(), localScores.cend()));
				// divide 4-byte float block max with quantile of term
				float temp = blockmax_float/quantile;
				// get quantized value (0-255)
				unsigned int quantized_value = (unsigned int) temp;
				// take the remainder from the division
				float remainder = temp - quantized_value;
				blockIdToMaxScore[block] = ( (Fcompare(remainder, 0.0f) != 0)&&(quantized_value!=255) ) ? (float) quantile*(quantized_value + 1) : (float) quantile*quantized_value;
			}

			localScores.clear();
			block = cblock;
		}
		localScores.push_back(scores[i]);
	}

	if(block == numOfBlocks)
		CERR << "ha?" << Log::endl;
	if(localScores.size()) {
		// get float block max
		float blockmax_float = *(std::max_element(localScores.cbegin(), localScores.cend()));
		// divide 4-byte float block max with quantile of term
		float temp = blockmax_float/quantile;
		// get quantized value (0-255)
		unsigned int quantized_value = (unsigned int) temp;
		// take the remainder from the division
		float remainder = temp - quantized_value;
		blockIdToMaxScore[block] = ( (Fcompare(remainder, 0.0f) != 0)&&(quantized_value!=255) ) ? (float) quantile*(quantized_value + 1) : (float) quantile*quantized_value;
	}
}

// Quantization
void SingleHashBlocker::generate_Quantized_Maxes(const std::vector<float>& scores, float& quantile) {
	size_t block = 0;
	std::vector<float> localScores;
	const size_t numOfBlocks(1+whereDidMapped(CONSTS::MAXD));
	assert(scores.size() == docIds->size());

	for(size_t i = 0; i<numOfBlocks;  ++i)
		QuantizedblockIdToMaxScore.push_back((unsigned char) 0);

	for(size_t i=0; i<docIds->size() && docIds->at(i)<CONSTS::MAXD; ++i) {
		const unsigned int cblock =  whereDidMapped((*docIds)[i]);
		if(cblock != block) {
			if(localScores.size()) {
				// get float block max
				float blockmax_float = *(std::max_element(localScores.cbegin(), localScores.cend()));
				// divide 4-byte float block max with quantile of term
				float temp = blockmax_float/quantile;
				// get quantized value (0-255)
				unsigned int quantized_value = (unsigned int) temp;
				// take the remainder from the division
				float remainder = temp - quantized_value;
				QuantizedblockIdToMaxScore[block] = ( (Fcompare(remainder, 0.0f) != 0)&&(quantized_value!=255) ) ? (unsigned char) quantized_value + 1 : (unsigned char) quantized_value;
				//debug
				//std::cout << block << " block has blockmax float: " << blockmax_float << " and quantized value: " << (unsigned int) QuantizedblockIdToMaxScore[block] << std::endl;
			}
			localScores.clear();
			block = cblock;
		}
		localScores.push_back(scores[i]);
	}

	if(block == numOfBlocks)
		CERR << "ha?" << Log::endl;
	if(localScores.size()) {
		// get float block max
		float blockmax_float = *(std::max_element(localScores.cbegin(), localScores.cend()));
		// divide 4-byte float block max with quantile of term
		float temp = blockmax_float/quantile;
		// get quantized value (0-255)
		unsigned int quantized_value = (unsigned int) temp;
		// take the remainder from the division
		float remainder = temp - quantized_value;
		QuantizedblockIdToMaxScore[block] = ( (Fcompare(remainder, 0.0f) != 0)&&(quantized_value!=255) ) ? (unsigned char) quantized_value + 1 : (unsigned char) quantized_value;
	}
}

void PostingsBlockerWithScores::generateMaxes(const std::vector<float>& scores) {
	for(BlocksIterator it=this->cbegin(); it!=this->cend(); ++it) {
		float max = AbsBlocksGen::getMaxScoreOfBlock(it,scores);
		for(size_t i=0; i<(*it).size(); ++i) {
		//	blockIdToMaxScore1[docIds->at((*it)[i])] = max;
		//	COUT << docIds->at((*it)[i]) << "->" << blockIdToMaxScore1[docIds->at((*it)[i])] << Log::endl;
		}
	}
}
