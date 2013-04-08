/*
 * SqlProxy.h
 *
 *  Created on: Feb 16, 2012
 *      Author: sergeyn
 */

//What a !@#@! mess
#ifndef _SQLPROXY_H_
#define _SQLPROXY_H_

#include <vector>
#include <iostream>
#include <algorithm>
#include "globals.h"
#include "profiling.h"
#include "sql/sqlite3.h"

class SqlProxy {
	sqlite3 *db;
	sqlite3_stmt *ppStmt;
	std::string dbname;

public:
	sqlite3_stmt * prepare(const std::string& select) {
	  const char *zSqlSelect = select.c_str();
	  if( sqlite3_prepare_v2(db, zSqlSelect, -1, &ppStmt, NULL) != SQLITE_OK )  {
		  CERR << "db error: " << sqlite3_errmsg(db) << Log::endl;
	      sqlite3_close(db);
	      ppStmt = 0;
		  return 0;
	  }

	  return ppStmt;
	}


	SqlProxy(std::string path=CONSTS::sqlPath) :db(0), ppStmt(0), dbname(path) {
		  if(sqlite3_open(dbname.c_str(), &db) ) {
			    CERR << "Can't open database: " << dbname << Log::endl;
				throw 1;
		  }
		  doInsertOrCreate("CREATE  TABLE  IF NOT EXISTS 'terms'('term' TEXT PRIMARY KEY  NOT NULL  UNIQUE, 'flagB' BLOB, 'maxB' BLOB, 'minB' BLOB, 'scoreB' BLOB, 'sizeB' BLOB)");
	}
	~SqlProxy() {
	  sqlite3_exec(db, "END", NULL, NULL, NULL);
	  sqlite3_close(db);
	}
	sqlite3 * getDB() { return db; }




	/*
	bool doInsertOrCreate(sqlite3_stmt *p) {
		 int rc =  sqlite3_step(p);
		 sqlite3_finalize(p);
		 return rc == SQLITE_OK;
	}
*/
	bool executeStoredStatement(bool doClose=false) {
		 int rc = sqlite3_step(ppStmt);
		 sqlite3_finalize(ppStmt);

		  if(doClose) {
			  sqlite3_exec(db, "END", NULL, NULL, NULL);
			  sqlite3_close(db);
		  }
		  return rc == SQLITE_OK;

	}

	bool doInsertOrCreate(const std::string& query, bool doClose=false) {

		  const char *zSql = query.c_str();
		  if( sqlite3_prepare_v2(db, zSql, -1, &ppStmt, NULL) != SQLITE_OK )  {
		      CERR << "db error: " << sqlite3_errmsg(db) << Log::endl;
		      sqlite3_close(db);
		      ppStmt = 0;
			  return false;
		  }

			 int rc = sqlite3_step(ppStmt);
			 sqlite3_finalize(ppStmt);


		  if(doClose) {
			  sqlite3_exec(db, "END", NULL, NULL, NULL);
			  sqlite3_close(db);
		  }
		  return rc == SQLITE_OK;
	}

};
#endif
