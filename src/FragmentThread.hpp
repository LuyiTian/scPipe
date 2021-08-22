#include <string>
#include <map>
#include <Rcpp.h>
#include <htslib/sam.h>

#ifndef FRAGMENT_THREAD_H
#define FRAGMENT_THREAD_H
struct FragmentStruct {
	std::string chromosome;
	int start;
	int end;
	std::string cell_barcode;
	bool complete;
};

// Class for a single FragmentThread
// holds input information as well as result information
//  
class FragmentThread {

	private:
		// fetchCall is given as callback to samfetch to call for every segment
		// the parent FragmentThread is passed in through second arg
		// typedef int (*bam_fetch_f)(const bam1_t *b, void *data)
		// typedef for int function called bam_fetch_f which takes a bam1_t and void *
		static int fetchCall(const bam1_t *, void *);

	public:
		FragmentThread(
			std::string,
			std::string,
			unsigned int,
			std::string,
			unsigned int,
			std::string,
			std::string,
			Rcpp::CharacterVector,
			unsigned int,
			unsigned int,
			unsigned int
		);

		std::string outname;
		std::string contig;
		unsigned int end;
		std::string bam;
		unsigned int min_mapq;
		std::string cellbarcode;
		std::string readname_barcode;
		Rcpp::CharacterVector cells;
		unsigned int max_distance;
		unsigned int min_distance;
		unsigned int chunksize;
		unsigned int fragment_count;

		std::map<std::string, FragmentStruct> fragment_dict;
		
		// equivalent to  getFragments in sinto_fragments.py
		void operator() ();

		void addToFragments(
			std::string,
			std::string,
			int,
			int,
			std::string,
			bool);

		//void writeFragments();

		//void collapseFragments();

		//void findCompleteFragments();

		void updateFragmentDict();
};

#endif