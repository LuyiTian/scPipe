#include <string>
#include <map>
#include <vector>
#include <Rcpp.h>
#include "bam.h"
#include <htslib/sam.h>
#include "Fragments.hpp"

#ifndef FRAGMENT_THREAD_H
#define FRAGMENT_THREAD_H

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

	FragmentMap fragment_dict;

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
	
	static void writeFragmentsToFile(std::vector<FragmentStruct> &, std::string);

	static std::map<std::string, int> CounterMapFragment(
		FragmentMap &, std::function<std::string(FragmentStruct&, bool, bool, bool, bool)>);
	static std::map<std::string, int> CounterMapString(std::vector<std::string> &);

	static std::map<std::string, std::vector<std::string>> createPositionLookup(
		std::vector<std::string> &, bool);
	static std::map<std::string, std::vector<std::string>> createPositionLookup(
		std::map<std::string, int> &, bool);

	//////////////////////////////////////////////////////
	//////// Instance methods for FragmentThread   ///////
	//////////////////////////////////////////////////////
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

	FragmentMap findCompleteFragments(unsigned int);

	///////////// Utility Functions for Writing //////////////////
	// update the internal fragment_count and return true if 
	// we need to write to file, otherwise false
	bool updateFragmentCount();

	void writeFragments(const unsigned int);



};

#endif


/// Documentation for used types and functions:
/*! @typedef
 @abstract Structure for core alignment information.
 @field  tid     chromosome ID, defined by bam_hdr_t
 @field  pos     0-based leftmost coordinate
 @field  bin     bin calculated by bam_reg2bin()
 @field  qual    mapping quality
 @field  l_qname length of the query name
 @field  flag    bitwise flag
 @field  l_extranul length of extra NULs between qname & cigar (for alignment)
 @field  n_cigar number of CIGAR operations
 @field  l_qseq  length of the query sequence (read)
 @field  mtid    chromosome ID of next read in template, defined by bam_hdr_t
 @field  mpos    0-based leftmost coordinate of next read in template
 */
// typedef struct {
//     int32_t tid;
//     int32_t pos;
//     uint16_t bin;
//     uint8_t qual;
//     uint8_t l_qname;
//     uint16_t flag;
//     uint8_t unused1;
//     uint8_t l_extranul;
//     uint32_t n_cigar;
//     int32_t l_qseq;
//     int32_t mtid;
//     int32_t mpos;
//     int32_t isize;
// } bam1_core_t;

/*! @typedef
 @abstract Structure for one alignment.
 @field  core       core information about the alignment
 @field  l_data     current length of bam1_t::data
 @field  m_data     maximum length of bam1_t::data
 @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux

 @discussion Notes:

 1. qname is terminated by one to four NULs, so that the following
 cigar data is 32-bit aligned; core.l_qname includes these trailing NULs,
 while core.l_extranul counts the excess NULs (so 0 <= l_extranul <= 3).
 2. l_qseq is calculated from the total length of an alignment block
 on reading or from CIGAR.
 3. cigar data is encoded 4 bytes per CIGAR operation.
 4. seq is nybble-encoded according to bam_nt16_table.
 */
// typedef struct {
//     bam1_core_t core;
//     int l_data;
//     uint32_t m_data; 
//     uint8_t *data;
// #ifndef BAM_NO_ID
//     uint64_t id;
// #endif
// } bam1_t;
//
// typedef BGZF *bamFile;
//
// typedef hts_idx_t bam_index_t;
// bam_index_t *bam_index_load(const chart *fn);
//
// typedef int (*bam_fetch_f)(const bam1_t *b, void *data)


// void print_core(const bam1_core_t *core) {
// 	Rcpp::Rcout << "Core:\n\ttid " << core->tid << "\n";
// 	Rcpp::Rcout << "\tpos " << core->pos << "\n";
// 	Rcpp::Rcout << "\tbin " << core->bin << "\n";
// 	Rcpp::Rcout << "\tqual " << core->qual << "\n";
// 	Rcpp::Rcout << "\tl_qname " << core->l_qname << "\n";
// 	Rcpp::Rcout << "\tflag " << core->flag << "\n";
// 	Rcpp::Rcout << "\tunused1 " << core->unused1 << "\n";
// 	Rcpp::Rcout << "\tl_extranul " << core->l_extranul << "\n";
// 	Rcpp::Rcout << "\tn_cigar" << core->n_cigar << "\n";
// 	Rcpp::Rcout << "\tl_qseq " << core->l_qseq << "\n";
// 	Rcpp::Rcout << "\tl_mtid " << core->mtid << "\n";
// 	Rcpp::Rcout << "\tmpos " << core->mpos << "\n";
// 	Rcpp::Rcout << "\tisize" << core->isize << "\n";
// }
// int
// 	SimpleFragFunc(const bam1_t *b, void *data) {
// 		FragmentThread *frag = (FragmentThread *)data;

// 		std::ofstream myfile;
// 		myfile.open(frag->outname, std::ofstream::app);
// 		myfile << "THIS IS AN OUTPUT FILE\n";
// 		myfile.close();

// 		return 1;
// 	}