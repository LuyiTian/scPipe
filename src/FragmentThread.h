#ifndef FRAGMENT_THREAD_H
#define FRAGMENT_THREAD_H


#include <string>
#include <map>
#include <vector>
#include <Rcpp.h>
#include <memory>

#include "bam.h"
#include <htslib/sam.h>
#include "ThreadOutputFile.h"
#include "Fragments.h"


class FragmentThread {
	public:
		FragmentThread(
			std::string,
			std::string,
			int,
			unsigned int,
			std::string,
			bam_header_t *,
			unsigned int,
			std::string,
			std::string,
			std::vector<std::string>,
			unsigned int,
			unsigned int,
			unsigned int
		);

		FragmentThread(const FragmentThread &);

		//////////////////////////////////////////////////////
		////////   Static methods for FragmentThread   ///////
		//////////////////////////////////////////////////////
		
		// fetchCall is given as callback to samfetch to call for every segment
		// the parent FragmentThread is passed in through second arg
		// typedef int (*bam_fetch_f)(const bam1_t *b, void *data)
		// typedef for int function called bam_fetch_f which takes a bam1_t and void *
		static int fetchCall(const bam1_t *, void *);

		static std::vector<FragmentStruct> collapseFragments(FragmentMap &);

		static std::map<std::string, int> collapseOverlapFragments(std::map<std::string, int> &, bool);

		static std::map<std::string, int> CounterMapFragment(
			FragmentMap &, std::function<std::string(FragmentStruct&, bool, bool, bool, bool)>);
		static std::map<std::string, int> CounterMapString(std::vector<std::string> &);

		static std::map<std::string, std::vector<std::string> > createPositionLookup(
			std::vector<std::string> &, bool);
		static std::map<std::string, std::vector<std::string> > createPositionLookup(
			std::map<std::string, int> &, bool);

		//////////////////////////////////////////////////////
		//////// Instance methods for FragmentThread   ///////
		//////////////////////////////////////////////////////
		void operator() ();

		void addToFragments(
				std::string,
				std::string,
				int,
				int,
				std::string,
				bool);

		void updateFragmentDict(const bam1_t *);

		FragmentMap findCompleteFragments(unsigned int);

		bool updateFragmentCount();

		void completeCollapseAndWriteFragments(const unsigned int);
	
	public:
		std::string contig;
		int tid;
		unsigned int end;
		std::string bam;
		bam_header_t *bam_header;
		unsigned int min_mapq;
		std::string cellbarcode;
		std::string readname_barcode;
		std::vector<std::string> cells;
		unsigned int max_distance;
		unsigned int min_distance;
		unsigned int chunksize;
		
		unsigned int fragment_count;

		FragmentMap fragment_dict;
		
	private:

		ThreadOutputFile debug;

		ThreadOutputFile fragfile;
};

#endif
