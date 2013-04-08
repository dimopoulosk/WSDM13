/*
 * globals.h
 *
 *  Created on: Jan 15, 2012
 *      Author: sergeyn
 */

#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <ostream>
#include <vector>
#include <assert.h>
#include <stdint.h>

//#define PROFILING
#ifdef GPROFILING
#define MYNOINLINE __attribute__((noinline))
#define PROFILER(a) //
#else
#define MYNOINLINE
#define PROFILER(a) profilerC::getInstance().stepCounter(a);
#endif

#ifdef CPP0X
	#include <unordered_map>
	#define hashMapStd std::unordered_map
	typedef std::unordered_map<std::string, int> termsMap;
#else
	#include <map>
	typedef std::map<std::string, int> termsMap;
	#define hashMapStd std::map
#endif




namespace Log {
	enum verbLevels { VALL, VDEBUG, VPROFILING, VOUTPUT };
#ifdef GPROFILING
#define COUT1 //
#define COUT2 //
#else
	#define COUT Log::logger()
	#define COUT1 Log::logger() << Log::verb<Log::VDEBUG>
	#define COUT2 Log::logger() << Log::verb<Log::VPROFILING>
	#define COUT3 Log::logger() << Log::verb<Log::VOUTPUT>
	#define COUT4 Log::logger(std::cout,4)
	#define COUT5 Log::logger(std::cout,5)
	#define CERR Log::logger(std::cerr,100)
	#define EFATAL Log::endl << " " << __FILE__ << " " << __LINE__ << Log::fatal
#endif
	void setGlobalVerbosityForAllLoggers(int v); //{ logger::GLOBAL_VERBOSITY = v; }

	typedef std::ostream outstrType; //just in case

	class logger {
		outstrType& out; //the stream of this logger
		int verbosity; //current verbosity
		static int GLOBAL_VERBOSITY;
	public:

		logger(outstrType& s=std::cout, int v=0) : out(s), verbosity(v) {}
		static void setGlobalV(int v) { GLOBAL_VERBOSITY = v; }
		template <typename T>
		logger& operator<<(const T& rhs) {
			if(verbosity >= GLOBAL_VERBOSITY)
			  out <<  rhs;
			return *this;
		 }

		//manipulators of the stream types
		typedef logger& (*endlType)(logger&);
		typedef logger& (*setVerbType)(logger&, int& v);
		logger& operator<<(endlType manip){ return manip(*this); }
		logger& operator<<(setVerbType manip){ return manip(*this,verbosity); }

		void setV(int v) { verbosity = v;}
		void flush() { if(verbosity >= GLOBAL_VERBOSITY) out << std::endl; }
	};

	logger& endl(logger& stream);  //stream flush
	logger& fatal(logger& stream);  //print and throw

	template<int val> //verbosity setter
	logger& verb(logger& stream, int& v) {
		v=val;
		return stream;
	}


/*

  usage
	COUT << "1 " << 2 << " " << 3.1 << log::endl << "new line" << log::endl;
	Log::setGlobalVerbosityForAllLoggers(1);
	CERR << "1 " << 2 << " " << 3.1 << log::endl << "new line" << log::endl;
	log::logger() << log::verb<5> << "1 " << 2 << " " << 3.1 << log::endl << "new line" << log::endl;
	COUT5 << "1 " << 2 << " " << 3.1 << log::endl << "new line" << log::endl;
*/

};

namespace CONSTS {

	enum counters { CNT_EVICTIONS,
		WAND,GETBLOCKSCORE,GETBLOCKSCORE1, BLOCKFAIL,ALIGN,MISALIGNED, GETFREQ, EVAL,HEAPIFY,HEAPIFYOK,NEXTGEQ,NEXTGEQ1,ALLQS,OTFBMG,
		SKIPS, SKIPSB, DECODE, CAND, MISTRUST, ESSENTIAL, SHALLOWPOINTER, SHALLOWPOINTER1, SORT, DOCIDFILTER1, DOCIDFILTER2,
		EARLYTERMINATION0, EARLYTERMINATION01, EARLYTERMINATION1, EARLYTERMINATION2, EARLYTERMINATION3, EARLYTERMINATION4, EARLYTERMINATION5,
		NONESSENTIAL,STEP1, STEP2, STEP3, NOP };

	const std::string QLogPath("../QueryLog/");
	const std::string fResultLog("../result_log");
	const std::string fQueryScore("../query_score");
	const std::string Golden_Threshold("../QueryLog/Golden_Threshold_trec06"); // Golden_Threshold_trec06 // Golden_Threshold_trec05
	const std::string zeroMapping("../QueryLog/10000q_terms.mapping"); //10000q_terms.mapping CHANGED !!
	const std::string termsMapping("../QueryLog/10000q_terms.mapping");
	// 10000q_terms.mapping // normal and reordered
	// 10000q_terms.mapping_layered // layered
	// 10000q_terms.mapping_layered_sorted // reordered + layered
	// version2: 10000q_terms_mapping_layered_version2 // layered (version2)
	// TREC05
	// trec05_mapping
	const std::string termsMapping_Layer("../QueryLog/10000q_terms.mapping_layered");

	// not used
	const std::string Index_term_mapping("word_file");
	const std::string entireIndexListLengths("../QueryLog/Entire_index_list_lengths");
	const std::string Unpadded_10k_postings("../QueryLog/10k_pair_term_unpadded_postings");
	const std::string Entire_Index_Maxscore_Unpadded_List_length("../QueryLog/Entire_Index_Maxscore_Unpadded_List_length"); //../QueryLog/
	const std::string Starting_Threshold("../QueryLog/10th_Starting_Threshold");
	const std::string Top1_Threshold("../QueryLog/Top-1_Threshold");
	const std::string Top10_Threshold("../QueryLog/Top-10_Threshold");

	const std::string ikQueryScores("../QueryLog/1k_scores_float");  // 1k_scores
	const std::string ikQuery("../QueryLog/1000query"); //1k_trec05; //1000query
	const std::string trecRoot("/research/index/trec06/");
	// indexes
	// /data5/constantinos/trec06/
	// /data5/constantinos/trec06_sorted/
    // /data5/constantinos/trec06_layer/
	// /data5/constantinos/trec06_sorted_layer/
	// /data5/constantinos/trec06_layer2/

	// Other options for index
	// Classical Index
	// /home/sergeyn/BMW/index64trec06/
	// /data2/BMW/BACKUPindex64trec06/ --
	// home/constantinos/Desktop/BMW_next/index/index64trec06/ (external disk)
	// correct one:
	// /data2/BMW/trec06/
	// SORTED
	// /data2/BMW/index64trec06/
	// Layered
	// /data2/BMW/trec06_Layers/

	#define EXPECTED_BLOCKSIZE "64"
	const int EXPECTED_BLOCK_SIZE(atoi(EXPECTED_BLOCKSIZE));
//const std::string trecRoot("/data1/Indexing/index" + std::string(EXPECTED_BLOCKSIZE) + "/"); // index location
// for BMW const std::string trecRoot("/data2/BMW/BMW_Code_Backup/index64trec06/"); or //const std::string trecRoot("/data/sding/BMW_Code_Backup/index64trec06/");
// ("/data1/Indexing/index" + std::string(EXPECTED_BLOCKSIZE) + "/");

	const std::string trecRawRoot("/data/sding/TwoLevelTrec/index_result/");///home/sergeyn/BMW/index_result/");
	// /data/sding/TwoLevelTrec/index_result/
	// media/Book_/Indexing Backup/index_result/
	// /home/sergeyn/BMW/index_result/

	const std::string doclenFileName(trecRoot+"doclen_file"); // /data/sding/TwoLevelTrec/index_result/doclen_file or doclen_file_sorted
	const std::string sqlPath(trecRoot+"indx.sqlite");
	const std::string basic_table(trecRoot+"basic_table");
	const std::string TERM_PATH(trecRoot+"pool");
	const std::string FLAG_PATH(trecRoot+"flag");
	const std::string MAX_PATH(trecRoot+"max");
	const std::string SCORE_PATH(trecRoot+"score");
	const std::string SIZE_PATH(trecRoot+"size");

	const std::string GOOD_TERM_SUFFIX("_good");
	const std::string BAD_TERM_SUFFIX("_bad");

	const std::string INFO_INDEX("8_1.inf");
	const std::string INDEX("8_1.dat");
	const std::string SORTED("_sorted");
	const std::string DOCUMENT_LENGTH("doclen_file");
	const std::string MERGED_BOOL_PATH("merged_bool/");
	const std::string WORD_FILE("/data/sding/TwoLevelTrec/index_result/word_file");

	const int MAXD(25205179);
	const int MAXTERMS(32818983);
	const int MAXDBITS(25);
	const int TOPK(10);
	const unsigned int STORAGE_CAPACITY(218); //2187 -- the lowest value to have 0 evictions on our 1K queries sample
	const int BS(64);

	const unsigned int bitvector_size = 1<<15;

	// Quantization
	const int Quantization(255);

	// Layering
	const int LAYER_SPLIT_POSTINGS_THRESHOLD(4096); // 65536 : Version1 (BMW): 50000 // Version2: 4096
	const float LAYER_SCORE_PARAMETER(0.25); // Version1 (BMW): 0.02  // Version2: 0.25  //0.5f // 0.25f

	const double FRAC(0.10); //no idea what those two are...
	const int S(16);

};

typedef	std::vector<unsigned int> vecUInt;
typedef	std::vector<unsigned int>::iterator vecUIntItr;
typedef	std::vector<unsigned int>::const_iterator vecUIntCItr;


#endif /* GLOBALS_H_ */
