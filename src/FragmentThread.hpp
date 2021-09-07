#include <string>
#include <map>
#include <Rcpp.h>
#include "bam.h"
#include <htslib/sam.h>

#ifndef BAM_UTIL_FUNCS_H
#define BAM_UTIL_FUNCS_H

// seqment->data is qname-cigar-seq-qual-aux
// qname is at bam1_qname(seqment) (which is (char*)seqment->data)
			
/*
Get the mapping quality of the alignment
@param b bam1_t pointer to an alignment
@return uint8_t mapping quality of alignment
*/
#define bam_mapping_qual(b) (b->core.qual)

/*
Get the 0-based leftmost coordinate of an alignment
@param b bam1_t pointer to an alignment
@return int32_t pos
*/
#define bam_alignment_start(b) (b->core.pos)

/*
Convert the uint8_t aux value to a string
Shorter version of htslib/sam.h bam_aux2Z
Removed return of 0 (null terminator) for bad input,
as Rcpp could not handle it and crashed
Now, valid input can be checked by the returning string's length
@param s uint8_t * pointer to tag data returned by bam_aux_get()
@return std::string of tag data, or empty string if tag type is not Z
*/
#define bam_aux2string(s) ( (s != NULL) ? std::string((*s++ == 'Z') ? (char *)s : "") : "")

#endif

#ifndef FRAGMENT_THREAD_H
#define FRAGMENT_THREAD_H
struct FragmentStruct {
	std::string chromosome;			// 0
	int32_t start;					// 1
	int32_t end;					// 2
	std::string cell_barcode;		// 3
	bool complete;					// 4
};

// Class for a single FragmentThread
// holds input information as well as result information
//
class FragmentThread {

private:

public:
	FragmentThread(
		std::string,
		std::string,
		int,
		unsigned int,
		std::string,
		bam_header_t *,
		uint8_t,
		std::string,
		std::string,
		Rcpp::CharacterVector,
		unsigned int,
		unsigned int,
		unsigned int
	);

	std::string outname;
	std::string contig;
	int tid;
	unsigned int end;
	std::string bam;
	bam_header_t *bam_header;
	uint8_t min_mapq;
	std::string cellbarcode;
	std::string readname_barcode;
	Rcpp::CharacterVector cells;
	unsigned int max_distance;
	unsigned int min_distance;
	unsigned int chunksize;
	unsigned int fragment_count;

	std::map<std::string, FragmentStruct> fragment_dict;
	
	// fetchCall is given as callback to samfetch to call for every segment
	// the parent FragmentThread is passed in through second arg
	// typedef int (*bam_fetch_f)(const bam1_t *b, void *data)
	// typedef for int function called bam_fetch_f which takes a bam1_t and void *
	static int fetchCall(const bam1_t *, void *);

	static std::map<std::string, FragmentStruct> 
		*findCompleteFragments(std::map<std::string, FragmentStruct> *, unsigned int, int32_t);
	
	static void writeFragmentsToFile(std::map<std::string, FragmentStruct> *, std::string);

	static std::map<std::string, FragmentStruct> *collapseFragments(std::map<std::string, FragmentStruct> *);

	// equivalent to  getFragments in sinto_fragments.py
	void operator() ();

	void addToFragments(
			std::string,
			std::string,
			int,
			int,
			std::string,
			bool);

	void updateFragmentDict(const bam1_t *);

	///////////// Utility Functions for Writing //////////////////
	// update the internal fragment_count and return true if 
	// we need to write to file, otherwise false
	bool updateFragmentCount();

	void writeFragments(const bam1_t *);



};

#endif
