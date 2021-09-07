#include <string>
#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <regex>
#include "bam.h"
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "FragmentThread.hpp"

FragmentThread::FragmentThread(
	std::string _outname,
	std::string _contig,
	int _tid,
	unsigned int _end,
	std::string _bam,
	bam_header_t *_bam_header,
	uint8_t _min_mapq,
	std::string _cellbarcode,
	std::string _readname_barcode,
	Rcpp::CharacterVector _cells,
	unsigned int _max_distance,
	unsigned int _min_distance,
	unsigned int _chunksize
) {
	this->outname = _outname;
	this->contig = _contig;
	this->tid = _tid;
	this->end = _end;
	this->bam = _bam;
	this->bam_header = _bam_header;
	this->min_mapq = _min_mapq;
	this->cellbarcode = _cellbarcode;
	this->readname_barcode = _readname_barcode;
	this->cells = _cells;
	this->max_distance = _max_distance;
	this->min_distance = _min_distance;
	this->chunksize = _chunksize;
	
	this->fragment_count = 0;
}


// fetchCall gets called for every segment in the bam file,
// so multiple times for one FragmentThread
// the parent FragmentThread class is passed in through data?
// bam1_t is struct containing the info for this specific segment
int
	FragmentThread::fetchCall(const bam1_t *b, void *data) {
		FragmentThread *frag = (FragmentThread *)data;

		frag->updateFragmentDict(b);

		if (frag->updateFragmentCount()) {
			// Find the complete fragments, remove from the fragment dict,
			// collapse them and then write to file
			frag->writeFragments(b);
		}

		return 1; // safe return value
	}

//' Update dictionary of ATAC fragments
//'
//' @description Takes a new aligned seqment and adds information to the dictionary
//' Modifies the original dictionary
void
	FragmentThread::updateFragmentDict(const bam1_t *seqment) {
		// If we have a regex to match against
		// using readname_barcode as a regex
		std::string cell_barcode;
		if (this->readname_barcode.length() > 0) {
			std::string qname(bam1_qname(seqment));
			std::smatch res;
			std::regex readname_regex (readname_barcode);

			// do the regex searching, which results in an smatch object
			// holding the matches information
			std::regex_search(qname, res, readname_regex);

			cell_barcode = *(res.begin()); // begin() is first match, with other matches after
		
		} else {
			// get the tag data associated with the specified tag
			// default for cellbarcode is "CB". Must be 2 characters
			// get raw uint8_t data and then convert to string
			uint8_t *raw_data = bam_aux_get(seqment, this->cellbarcode.substr(0, 2).c_str());
			std::string barcode_data = bam_aux2string(raw_data); // FragmentThread.hpp

			if (barcode_data.length() != 0) {
				cell_barcode = barcode_data;
			}
		}

		// if we have cells (to retain) and found a barcode
		// we can make sure that we actually should process this read
		// at all. If not, just return
		if (this->cells.length() != 0 && cell_barcode.length() != 0) {
			bool contains = false;
			for (auto it = cells.begin(); it != cells.end(); it++) {
				if (std::strcmp(cell_barcode.c_str(), *it) == 0) {
					contains = true;
					break;
				}
			}
			if (!contains) {
				return;
			}
		}

		uint8_t mapq = bam_mapping_qual(seqment); // FragmentThread.hpp
		// recording a fragment requires a minimum mapping quality
		if (mapq >= this->min_mapq) {
			char *qname = bam1_qname(seqment); // bam.h
			int32_t rstart = bam_alignment_start(seqment); // FragmentThread.hpp
			int32_t rend = bam_endpos(seqment); // sam.h
			bool is_reverse = bam1_strand(seqment); // bam.h

			if (rstart == -1 || rend == -1) {
				return;
			}

			// Correct start and end for 9bp Tn5 shift
			if (is_reverse) {
				rend = rend - 5;
			} else {
				rstart = rstart + 4;
			}

			this->addToFragments(
				qname,
				this->contig,
				rstart,
				rend,
				cell_barcode,
				is_reverse
			);
		}
	}

//' Add new fragment information to the fragment dictionary
//' Checks to see if a fragment with this qname is already in the map
//' If so, performs quality checks and removes or updates the frag
//' Otherwise, creates and adds a new FragmentStruct
//' @param qname read name
//' @param chromosome chromsome name
//' @param start alignment start postiion
//' @param rend alignment end position
//' @param cell_barcode cell barcode sequence
//' @param is_reverse read is aligned to reverse strand
void
	FragmentThread::addToFragments(
		std::string qname,
		std::string chromosome,
		int32_t rstart,
		int32_t rend,
		std::string cell_barcode,
		bool is_reverse) {

		// check if qname is already in fragment dict
		if (fragment_dict[qname].chromosome.length()) {
			if (is_reverse) {
				int32_t current_coord = fragment_dict[qname].start;
				if (current_coord == -1) {
					// read aligned to the wrong strand
					fragment_dict.erase(qname);
				} else if (((rend - current_coord) > this->max_distance) ||
					((rend - current_coord) < this->min_distance)) {
					// too far away, don't include
					fragment_dict.erase(qname);
				} else {
					if (cell_barcode.empty() && fragment_dict[qname].cell_barcode.empty()) {
						// both fragment ends are present but no cell barcode
						fragment_dict.erase(qname);
					} else {
						if (fragment_dict[qname].cell_barcode.empty()) {
							fragment_dict[qname].cell_barcode = cell_barcode;
						}
						fragment_dict[qname].end = rend;
						fragment_dict[qname].complete = true;
					}
				}
			} else {
				int32_t current_coord = fragment_dict[qname].end;
				if (current_coord == -1) {
					fragment_dict.erase(qname);
				} else if (((current_coord - rstart) > this->max_distance) ||
					((current_coord - rstart) < this->min_distance)) {
					fragment_dict.erase(qname);
				} else {
					if (cell_barcode.empty() && fragment_dict[qname].cell_barcode.empty()) {
						fragment_dict.erase(qname);
					} else {
						if (fragment_dict[qname].cell_barcode.empty()) {
							fragment_dict[qname].cell_barcode = cell_barcode;
						}
						fragment_dict[qname].start = rend;
						fragment_dict[qname].complete = true;
					}
				}
			}
		} else {
			// make a new fragment
			fragment_dict[qname] = {
				chromosome,						// chromosome
				(!is_reverse ? rstart : -1),	// start
				(is_reverse ? rend : -1),		// end
				cell_barcode,					// cell_barcode
				false							// complete
			};
		}
	}



//' FragmentThread() is equivalent to old getFragments
//'
//' Execution thread to iterate over paired reads in BAM file and extract
//' ATAC fragment coordinates
//' Multithreaded. Needs every argument to be thread safe???
//' Only shared resource is inbam,
//' All required parameters are in class constructor and private variables
void
	FragmentThread::operator() () {
		Rcpp::Rcout << "This is a fragment thread with contig : " << this->contig << "\n";

		bamFile bam = bam_open(this->bam.c_str(), "r"); // bam.h
		bam_index_t *index = bam_index_load(this->bam.c_str()); // bam.h
		bam_fetch(bam, index, this->tid, 0, this->end, this, &FragmentThread::fetchCall);
		
		
		//this->writeFragments(); wtf do we do here?

		bam_close(bam); // bam.h
	}



//////////////////////////// utility functions for writing //////////////////
void
FragmentThread::writeFragments(const bam1_t *b) {
	int32_t current_position = bam_alignment_start(b);

	std::map<std::string, FragmentStruct> *complete = FragmentThread::findCompleteFragments(
		this->fragment_dict,
		this->max_distance,
		current_position
	)
	// 	std::map<std::string, FragmentStruct> *collapsed = FragmentThread::collapseFragments(complete);
	
	FragmentThread::writeFragmentsToFile(collapsed, this->outmame);

	// gc
	delete complete;
	// delete collapsed???????
}


void 
FragmentThread::writeFragmentsToFile(
		const std::map<std::string, FragmentStruct> *fragments,
		std::string filepath
) {
	std::ofstream fragfile;
	fragfile.open(filepath, std::ofstream::app);

	for (auto frag = fragments->begin(); frag != fragments->end(); frag++) {
		fragfile << frag->second.chromosome << "\t"
			 << frag->second.start << "\t"
			 << frag->second.end << "\t"
			 << frag->second.cell_barcode << "\n";
		
	}

	fragfile.close();
}

// Find complete fragments that are > max_dist bp away from the current BAM file position
// @param fragments the dictionary containing ATAC fragment information
// @param max_distance maximum allowed distance between fragment start and end postions
// @param current_position the current_position being looked at in the position-sorted BAM file
// @param max_collapse_dist maximum allowed distance for fragments from the same clel
// barcode that share one Tn5 integration site to be collapsed into a single fragment
// 
// @description Creates a new std::map of completed fragments and removes the completed ones from the original
std::map<std::string, FragmentStruct> *findCompleteFragments(
	const std::map<std::string, FragmentStruct>& fragments,
	unsigned int max_distance,
	int32_t current_position
) {
	unsigned int max_collapse_dist = 20;

	std::map<std::string, FragmentStruct> *completed = new std::map<std::string, FragmentStruct>;
	unsigned int d = max_distance + max_collapse_dist;
	for (auto frag = fragments.begin(); frag != fragments.end(); frag++) {
		if (frag->second.complete) { 
			// if we have a complete fragment, add it to completed and delete from old
			if ((frag->second.end + d) < current_position) {
				(*completed)[frag->first] = frag->second;
				fragments.erase(frag->first);
			}
		} else {
			if (frag->second.start < 0) {
				if ((frag->second.end + d) < current_position) {
					fragments.erase(frag->first);
				}
			} else if (frag->second.end < 0) {
				if ((frag->second.start + d) < current_position) {
					fragments.erase(frag->first);
				}
			} else {
				fragments.erase(frag->first);
			}
		}
	}

	return completed;
}

static std::map<std::string, FragmentStruct> *collapseFragments(std::map<std::string, FragmentStruct> *fragments) {
	return fragments;
}

bool
FragmentThread::updateFragmentCount() {
	if (++(this->fragment_count) > this->chunksize) {
		this->fragment_count = 0;
		return true;
	}
	return false;
}


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