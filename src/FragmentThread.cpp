#include <string>
#include <Rcpp.h>
#include "bam.h"
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "FragmentThread.hpp"

FragmentThread::FragmentThread(
	std::string _outname,
	std::string _contig,
	unsigned int _end,
	std::string _bam,
	unsigned int _min_mapq,
	std::string _cellbarcode,
	std::string _readname_barcode,
	Rcpp::CharacterVector _cells,
	unsigned int _max_distance,
	unsigned int _min_distance,
	unsigned int _chunksize
) {
	this->outname = _outname;
	this->contig = _contig;
	this->end = _end;
	this->bam = _bam;
	this->min_mapq = _min_mapq;
	this->cellbarcode = _cellbarcode;
	this->readname_barcode = _readname_barcode;
	this->cells = _cells;
	this->max_distance = _max_distance;
	this->min_distance = _min_distance;
	this->chunksize = _chunksize;
	this->fragment_count = 0;
}

void print_core(const bam1_core_t *core) {
    Rcpp::Rcout << "Core:\n\ttid " << core->tid << "\n";
    Rcpp::Rcout << "\tpos " << core->pos << "\n";
    Rcpp::Rcout << "\tbin " << core->bin << "\n";
    Rcpp::Rcout << "\tqual " << core->qual << "\n";
    Rcpp::Rcout << "\tl_qname " << core->l_qname << "\n";
    Rcpp::Rcout << "\tflag " << core->flag << "\n";
    Rcpp::Rcout << "\tunused1 " << core->unused1 << "\n";
    Rcpp::Rcout << "\tl_extranul " << core->l_extranul << "\n";
    Rcpp::Rcout << "\tn_cigar" << core->n_cigar << "\n";
    Rcpp::Rcout << "\tl_qseq " << core->l_qseq << "\n";
    Rcpp::Rcout << "\tl_mtid " << core->mtid << "\n";
    Rcpp::Rcout << "\tmpos " << core->mpos << "\n";
    Rcpp::Rcout << "\tisize" << core->isize << "\n";
}

// fetchCall gets called for every segment in the bam file,
// so multiple times for one FragmentThread
// the parent FragmentThread class is passed in through data?
// bam1_t is struct containing the info for this specific segment
int
FragmentThread::fetchCall(const bam1_t *b, void *data) {
	Rcpp::Rcout << "Inside fetch function\n";

	//FragmentThread frag = *((FragmentThread *) data); // type cast and dereference
	
	//print_core(&(b->core));
	//Rcpp::Rcout << "\t" << frag.contig << "\n";

	//frag.updateFragmentDict();

	return 1; // safe return value
}

//' Update dictionary of ATAC fragments
//' 
//' @description Takes a new aligned seqment and adds information to the dictionary
//' Modifies the original dictionary
//' @param seqment bam1_t - samtools region
//' @param cell_barcode tag used for cell barcode. Default is CB
//' @param readname_barocde regex for matching cell barcode in read name
//' @param cells list of cells to retain. If empty, retain all cells found
void
FragmentThread::updateFragmentDict() {
	Rcpp::Rcout << "\tfragment dict\n";
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
	int rstart,
	int rend,
	std::string cell_barcode,
	bool is_reverse) {
	
	// check if qname is already in fragment dict
	if (fragment_dict[qname].chromosome.length()) {
		if (is_reverse) {
			int current_coord = fragment_dict[qname].start;
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
			int current_coord = fragment_dict[qname].end;
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

	Rcpp::Rcout << "Successfully opened the BAM file\n";

	bam_index_t *index = bam_index_load(this->bam.c_str()); // bam.h
	bam_header_t *header = bam_header_read(bam); // bam.h

	// get the tid from the header file
	// aparently this isn't thread safe. What are we to do?
	int tid = bam_get_tid(header, this->contig.c_str()); // bam.h

	// actually to the bam fetching
	// this is equivalent to
	// for i in inputBam.fetch(this->contig, 0, this->end):
	// 		FragmentThread::fetchCall(i, this)
	bam_fetch(bam, index, tid, 0, this->end, nullptr, &FragmentThread::fetchCall);

	bam_close(bam); // bam.h
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
