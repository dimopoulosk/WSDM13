#include "globals.h"
#include "ListIterator.h"
#include "PostingOriented_BMW.h"
#include "TrecReader.h"
#include <bitset>

#ifdef CPP0X
	typedef std::unordered_map<std::string, float> strToFloatMap;
#else
	typedef std::map<std::string, float> strToFloatMap;
#endif

class QueryLogManager {
	std::vector<std::vector<std::string> > queriesD;
	strToFloatMap mapQueryToScore;
	termsMap *lex;

public:
	class queriesFileIterator : public std::iterator<std::forward_iterator_tag, std::vector<std::string> > 	{
		const QueryLogManager& queriesP;
		int curr;
		size_t count;
		size_t limit;
		size_t bucket;

		bool isSkip() const;
	public:
	  queriesFileIterator(const QueryLogManager& lm, int cur=-1) :queriesP(lm), curr(cur), count(0), limit(0), bucket(0) {}
	  queriesFileIterator(const queriesFileIterator& it) : queriesP(it.queriesP), curr(it.curr){}
	  bool operator==(const queriesFileIterator& rhs) {return curr==rhs.curr;}
	  bool operator!=(const queriesFileIterator& rhs) {return curr!=rhs.curr;}
	  const std::vector<std::string>& operator*() {return queriesP[curr];}
	  float score() const { return curr<queriesP.size() ? queriesP.score(queriesP[curr]) : -1.0; }
	  void changeSkipPolicy(size_t lim, size_t buck) { limit=lim; bucket=buck;}
	  queriesFileIterator& operator++();

	};

	QueryLogManager(const char* fname, termsMap *l); //will add more loaders in future?
	QueryLogManager(const char* fname, termsMap *l, bool& nothing);
	void setScoreForQuery(const std::vector<std::string>& terms, float v);
	size_t size() const { return queriesD.size(); }
	const queriesFileIterator end() { return queriesFileIterator(*this,queriesD.size()); }
	const std::vector<std::string>& operator[] (size_t q) const { return queriesD[q]; }
	float score(const std::vector<std::string>& terms) const;
	void loadScoresFromFile(const char* fname);
};


void dumpResult(int qn, int rlen, const std::vector<std::string>& word_l, QpResult* res, FILE* fresult_log, FILE* fqscore );

class QueryProcessing {
	int qn;
	termsCache& Cache;
	termsMap termsInCache;
	termsMap mappingForOnDemand;
	unsigned int* pages;
	TrecReader* trecReader;
	unsigned int* loaddoclen(const char* fname=CONSTS::doclenFileName.c_str());
	lptrArray openLists(const std::vector<std::string>& word_l, termsMap& lex);

public:
	QueryProcessing(termsCache& cache);
	void onDemandCall(const std::string& term);
	void fillTermsMap(termsMap& lex, const std::string&);
	void operator()(const char* queryLog, const int buckets, const int limit, const int topk, const int layer);
	void printReport() const;

	// newly added (auxiliary) functions
	void Print_Stats(std::vector<std::string>& word_l, std::vector<float>& threshold_increase, int& number_of_essential_lists, double& qp_time, lptrArray& lps, QpResult* res, const int& topk);
	void Write_Inverted_Lists_in_File(std::vector<std::string>& word_l, lptrArray& lps);
	std::vector<int> Get_List_Lengths(termsMap& lex);
	std::vector<int> Load_Entire_Index_List_Lengths();
	std::vector<int> Load_Query_log_List_lengths(termsMap& lex);
	void Space_Comparison(std::vector<int>& List_Lengths);
	void List_Length_Distribution(std::vector<int>& List_Lengths);
	std::unordered_map<std::string, int> Load_unpadded_postings(const std::string& path);
	std::unordered_map<std::string, float> Load_Threshold(const std::string& path);
	void on_the_fly_max_array_generation(std::vector<float>& max_array, RawIndexList& rilist, int& bits);
	void translate_and_assign_ids_to_layered_terms(std::vector<std::string>& word_l, termsMap& lex, std::vector<std::string>& translated_terms, std::vector<int>& list_ids);
	void set_layered_metainfo_to_lps(lptrArray& lps, std::vector<int>& list_ids, termsMap& Unpadded_Length_Map, std::vector<std::string>& layered_terms);
	RawIndexList naive_lps_to_RawIndexList(lptr*& lps);
	RawIndexList lps_to_RawIndexList(lptr*& lps);
	void on_the_fly_max_array_generation(lptr*& lps, std::vector<float>& max_array, int& bits);
	void Term_To_list_length_Map(termsMap& lex, const std::string& path, const bool& unpadded_list_length);
	void Bucketized_Maxscore_Stats();
	float Get_Average_Document_Length();
	void Set_Quantiles(lptrArray& lps);
	void Query_log_Maxscore_Length_Correlation(termsMap& lex);
	void bitOracle(unsigned int& length, int& lenBits);
};
