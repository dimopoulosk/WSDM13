#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "TrecReader.h"
#include "qp.h"
#include "optionparser.h"

 enum  optionIndex { UNKNOWN, HELP, BUILD, ONDEMAND, OUTPATH, QP, QPATH, OFFSET, LIMIT, VERB, TOPK, LAYER, INDEXROOT };
 const option::Descriptor usage[] =
 {
  {UNKNOWN, 	0,"" , ""    ,		option::Arg::None, "USAGE: example [options]\n\nOptions:" },
  {HELP,    	0,"h", "help",		option::Arg::None, "  --help  \tPrint usage and exit." },
  {BUILD,   	0,"b", "buildindex",option::Arg::Optional, "  --buildindex[=<raw index path>], -b [<raw index path>]  \tBuild index with optional path to raw index" },
  {ONDEMAND,   	0,"d", "demandBld", option::Arg::Optional, "  --demandBld[=<raw index path>], -d [<raw index path>]  \tBuild term lists on demand with optional path to raw index" },
  {OUTPATH,   	0,"o", "outputpath",option::Arg::Optional, "  --outputpath[=<index path>], -o [<index path>]  \tWhere to store the index" },
  {QPATH,   	0,"p", "qpath",     option::Arg::Optional, "  --qpath[=<query log path>], -p [<query log path>]  \tWhere to get the queries from" },
  {QP,   		0,"q", "qp",		option::Arg::Optional, "  -q [<buckets>], \t--qp[=<buckets>] \tRun query proc. with optional buckets" },
  {LIMIT,  		0,"l", "limit",		option::Arg::Optional, "  -l [<limit>], \t--limit[=<limit>] \tRun at most limit queries" },
  {OFFSET, 		0,"s", "start",		option::Arg::Optional, "  -s [start], \t--start[=<start>] \tThe first index to start at -> for build, but when run selects the exp. block size" },
  {VERB, 		0,"v", "verbosity",	option::Arg::Optional, "  -v [int], \t--verbosity[=int] \tLogging level -- lower level - more messages" },
  {TOPK, 		0,"k", "topk",	    option::Arg::Optional, "  -k [int], \t-- number of top-k documents to retrieve" },
  {LAYER, 		0,"r", "layer",		option::Arg::Optional, "  -r [int], \t-- Layer off: 0, Layer on: 1" },
  {INDEXROOT, 	0,"i", "indexroot",	option::Arg::Optional, "  -i [str], \t--indexroot[=str] \tThe path to trecroot" },
  {0,0,0,0,0,0}
 };


int main(int argc, char * argv[] ){
	   argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
	   option::Stats  stats(usage, argc, argv);
	   option::Option options[stats.options_max], buffer[stats.buffer_max];
	   option::Parser parse(usage, argc, argv, options, buffer);

	   if (parse.error())
	     return 1;

	   if (options[HELP]) {
	     option::printUsage(std::cout, usage);
	     return 0;
	   }

	   if (options[VERB] && options[VERB].arg)
		   Log::setGlobalVerbosityForAllLoggers(atoi(options[VERB].arg));

	   int buckets = (options[QP] && options[QP].arg) ? atoi(options[QP].arg) : 0;
	   int limit = (options[LIMIT] && options[LIMIT].arg) ? atoi(options[LIMIT].arg) : CONSTS::MAXD;
	   int offset = (options[OFFSET] && options[OFFSET].arg) ? atoi(options[OFFSET].arg) : 0;
	   int topk = (options[TOPK] && options[TOPK].arg) ? atoi(options[TOPK].arg) : 10;
	   int layer = (options[LAYER] && options[LAYER].arg) ? atoi(options[LAYER].arg) : 0;

	   if(options[BUILD]) {
		   const std::string&  trecRawRoot = options[BUILD].arg ? std::string(options[BUILD].arg) : CONSTS::trecRawRoot;
		   std::string outputPath = (options[OUTPATH] && options[OUTPATH].arg) ? std::string(options[OUTPATH].arg): CONSTS::trecRoot;
		   stringIntVectorsPair tmap; //empty one
		   TrecReader* tr =  TrecFactory(trecRawRoot);
		   TrecFactory(*tr,outputPath, offset,limit,tmap);
		   delete tr;
		   return 0;
	   }

	   std::string queryFile = (options[QPATH] && options[QPATH].arg) ? std::string(options[QPATH].arg): CONSTS::ikQuery;

	   SingleHashBlocker::expectedBlockSize = (options[OFFSET] && options[OFFSET].arg) ? atoi(options[OFFSET].arg) : 8;
	   //COUT3 << "Expected block size: " << SingleHashBlocker::expectedBlockSize << Log::endl; // togo

	   //togo
	   //std::cout << "########### expected block size: 2^" << SingleHashBlocker::expectedBlockSize << std::endl;

	   if(nice(-3)) {}

	   termsCache Cache;
	   if(options[ONDEMAND] )
		   onDemandCpool::initPool(CONSTS::STORAGE_CAPACITY); //prepare the cache but load none
	   else {
		   // load different terms_mapping
		   if (layer == 0)
			   Cache.fill_cache(CONSTS::termsMapping.c_str()); //load terms //
		   else
			   Cache.fill_cache(CONSTS::termsMapping_Layer.c_str()); //load terms Layering
	   }

	   QueryProcessing qp(Cache);
	   qp(queryFile.c_str(), buckets, limit, topk, layer); //process the queries
	   qp.printReport();
	   COUT1 << "start cleanup" << Log::endl;
	   return 0;
}

