#ifndef _TREC_READER_
#define _TREC_READER_

#include <stdio.h>
#include <stdlib.h>

#include "globals.h"
#include "ListIterator.h"
#include "SqlProxy.h"

typedef struct resultItem
{
	int doc;
	float score;
} resultItem;

typedef std::pair<std::vector<std::string>, std::vector<size_t> > stringIntVectorsPair;

termsMap getTrecWordsMappingMap(const std::string& qMappingPath);

class TrecReader
{
public:
	TrecReader() {}
	// load the lex into mem, get the urls and doclength
	TrecReader(const std::string& index_path,
			   const std::string& lex_path,
			   const std::string& doclength_path,
			   const std::string&  url_path, int docn);

	void loadRawListIn(RawIndexList& tList);

	RawIndexList getRawList(const std::string& term, size_t wid);
	// Usage: Given the term and its term_id, return a RawIndexList structure (vector of scores, dids, freqs) - padded
	RawIndexList load_raw_list(const std::string& term, size_t wid);

	void dump_To_File_Pairs_of_Maxscore_and_Unpadded_List_Lengths(const stringIntVectorsPair& tmap, const std::string& path, int offset, int limit);

	// Usage: Creating Layered index
	// Input: Given the original rawindexlist construct the good and the bad term
	void build_layered_index(RawIndexList& original_Term, RawIndexList& good_Term, RawIndexList& bad_Term);
	int prepare_list_Layers(const stringIntVectorsPair& tmap, const std::string& rootPath, int offset, int limit);
	float get_split_threshold(std::vector<float>& scores);

	int prepare_list(const stringIntVectorsPair& mapping, const std::string& rootPath, int offset, int limit);
	void load_doclength();
	~TrecReader(void);
private:
	FILE *findex;
	FILE *finf;
	FILE *flex;
	FILE *fdoclength;

	int docn;
	unsigned int* doclen;
	SqlProxy sql;

//	int currentwid;
	unsigned int* inf_buffer;
	unsigned int* hold_buffer;

	size_t* inf_prefix_sum;
#ifdef CPP0X
	std::unordered_map<std::string, int*> lex;
#else
	std::map<std::string, int*> lex;
#endif
};

void TrecFactory(TrecReader& reader,const std::string& resultRootPath, int offset, int limit, stringIntVectorsPair& tmap);
TrecReader* TrecFactory(const std::string& index_res);


#endif
