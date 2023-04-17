#include "FragmentThread.h"

#include <Rcpp.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <sstream>
#include <regex>
#include <map>
#include <forward_list>
#include <numeric>
#include <mutex>
#include <memory>


#include "bam.h"
#include <htslib/bgzf.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "ThreadOutputFile.h"
#include "FragmentUtils.h"

FragmentThread::FragmentThread(
	// std::shared_ptr<ThreadOutputFile> _fragfile,
	std::string _fragfile,
	std::string _contig,
	int _tid,
	unsigned int _end,
	std::string _bam,
	bam_header_t *_bam_header,
	unsigned int _min_mapq,
	std::string _cellbarcode,
	std::string _readname_barcode,
	std::vector<std::string> _cells,
	unsigned int _max_distance,
	unsigned int _min_distance,
	unsigned int _chunksize
) {
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

	this->fragfile.setFile(_fragfile);

	// this->debug.setFile("/stornext/Home/data/allstaff/v/voogd.o/scPipeATACTiming/outputs/d3Haoyu/tmpDebug0");
}

FragmentThread::FragmentThread(const FragmentThread &old) {
	this->contig = old.contig;
	this->tid = old.tid;
	this->end = old.end;
	this->bam = old.bam;
	this->bam_header = old.bam_header;
	this->min_mapq = old.min_mapq;
	this->cellbarcode = old.cellbarcode;
	this->readname_barcode = old.readname_barcode;
	this->cells = old.cells;
	this->max_distance = old.max_distance;
	this->min_distance = old.min_distance;
	this->chunksize = old.chunksize;
	
	this->fragment_count = old.fragment_count;

	this->fragment_dict = FragmentMap(old.fragment_dict);

	// this->fragfile(old.fragfile.getPath());
	this->fragfile.setFile(old.fragfile.getPath());

	// this->debug.setFile(old.debug.getPath());
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

		cell_barcode = *(res.begin()); // begin() is first match
	
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
	if (this->cells.size() != 0 && cell_barcode.length() != 0) {
		bool contains = false;
		for (auto it = cells.begin(); it != cells.end(); it++) {
			if (std::strcmp(cell_barcode.c_str(), (*it).c_str()) == 0) {
				contains = true;
				break;
			}
		}
		if (!contains) {
			return;
		}
	}

	unsigned int mapq = (unsigned int) bam_mapping_qual(seqment); // FragmentThread.hpp
	// recording a fragment requires a minimum mapping quality
	if (mapq < this->min_mapq) {
		return;
	}

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
	if (fragment_dict.count(qname)) {
		if (is_reverse) {
			int32_t current_coord = fragment_dict[qname].start;
			if (current_coord == -1) {
				// read aligned to the wrong strand
				fragment_dict.erase(qname);
			} else if (((rend - current_coord) > (int32_t)this->max_distance) ||
				((rend - current_coord) < (int32_t)this->min_distance)) {
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
			// augment an existing fragment, using a non reversed fragment 
			// (using rstart to augment existing fragment start)
			int32_t current_coord = fragment_dict[qname].end;
			if (current_coord == -1) {
				// if we only have a fragment that contains a start coord, delete (will be replaced with a new frag)
				fragment_dict.erase(qname);
			} else if (((current_coord - rstart) > (int32_t)this->max_distance) ||
				((current_coord - rstart) < (int32_t)this->min_distance)) {
				// if we can't use this start coord, as it's too far beyond the end of the current fragment,
				// delete the old one.
				fragment_dict.erase(qname);
			} else {
				// if we can use the start coord:
				if (cell_barcode.empty() && fragment_dict[qname].cell_barcode.empty()) {
					// delete the old fragment if there's no barcode here or there
					fragment_dict.erase(qname);
				} else {
					// there's at least one barcode
					if (fragment_dict[qname].cell_barcode.empty()) {
						// if we can use the start coord, and there's no barcode
						fragment_dict[qname].cell_barcode = cell_barcode;
					}
					fragment_dict[qname].start = rstart;
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
void
FragmentThread::operator() () {

	bamFile bam = bam_open(this->bam.c_str(), "r"); // bam.h
	bam_index_t *index = bam_index_load(this->bam.c_str()); // bam.h

	bam_fetch(bam, index, this->tid, 0, this->end, this, &FragmentThread::fetchCall);

	// for the final writeFragments call, pass in inf 
	this->completeCollapseAndWriteFragments(4294967295); // 4294967295 is max unsigned int

	bam_close(bam); // bam.h
}


void
FragmentThread::completeCollapseAndWriteFragments(unsigned int current_position) {
	FragmentMap complete = this->findCompleteFragments(current_position);

	std::vector<FragmentStruct> collapsed = FragmentThread::collapseFragments(complete);

	this->fragfile.write(collapsed);
}


/// Increment the fragment_count member on this FragmentThread
/// If the count is greater than the chunksize, reset the count and return true to indicate
/// that the fragment dictionary needs to be written to the file
bool
FragmentThread::updateFragmentCount() {
	if (++(this->fragment_count) > this->chunksize) {
		this->fragment_count = 0;
		return true;
	}
	return false;
}



////////////////////////////////////////////////////////////////////////////////////////////////
/////       Below contains all static utility functions for the FragmentThread class   /////////
////////////////////////////////////////////////////////////////////////////////////////////////

// fetchCall gets called for every segment in the bam file,
// so multiple times for one FragmentThread
// bam1_t is struct containing the info for this specific segment
int
FragmentThread::fetchCall(const bam1_t *b, void *data) {
	FragmentThread *frag = (FragmentThread *)data;

	// update the fragment map with this segments data 
	// retrived from the bam1_t aligned segment
	frag->updateFragmentDict(b);

	if (frag->updateFragmentCount()) {
		// Find the complete fragments, remove from the fragment dict,
		// collapse them and then write to file
		frag->completeCollapseAndWriteFragments((unsigned int) bam_alignment_start(b));
	}

	return 1; // safe return value
}


std::vector<FragmentStruct>
FragmentThread::collapseFragments(FragmentMap &fragments) {
	// counts is a map where each key is the FragmentStruct information joined into a string
	// separated by '|'
	std::map<std::string, int> counts = CounterMapFragment(fragments, FragToString);

	// early break if we have an empty map
	if (counts.size() == 0) return std::vector<FragmentStruct> {};

	// map of fragment information associated with indices
	std::map<std::string, int> *frag_id_lookup = 
		id_lookup(fragments, 
					[](FragmentStruct frag)->std::string { 
						std::stringstream ss;
						ss << frag.chromosome << "|" << frag.start << "|" << frag.end;
						return ss.str();
					}
	);
	// map of barcode information associated with indices
	std::map<std::string, int> *bc_id_lookup = 
		id_lookup(fragments,
					[](FragmentStruct frag)->std::string {
						return frag.cell_barcode;
					}
	); 

	counts = collapseOverlapFragments(counts, true);
	counts = collapseOverlapFragments(counts, false);

	std::map<int, int> row_sum;
	std::map<int, std::pair <int, int> > row_max;
	int max_row = 0;
	for (auto item : counts) {
		// rowstr contains a string of "chromosome|start|end"
		FragmentStruct temp_s1 = StringToFrag(item.first);
		std::string rowstr = FragToString(temp_s1, true, true, true, false);
		// bcstr contains "cell_barcode"
		std::string bcstr = FragToString(temp_s1, false, false, false, true);

		int this_row = frag_id_lookup->at(rowstr);

		// save the size of the created sparse matrix
		if (this_row > max_row) {
			max_row = this_row;
		}

		row_sum[this_row] += item.second;

		int this_col = bc_id_lookup->at(bcstr);
		if (row_max[this_row].second < item.second) {
			row_max[this_row] = std::pair<int, int> {this_col, item.second};
		}

	}
	
	std::map<int, std::string> frag_inverse = invertMap(frag_id_lookup);
	std::map<int, std::string> bc_inverse = invertMap(bc_id_lookup);

	std::map<int, std::string> collapsed_frags;
	for (auto it : row_max) {
		collapsed_frags[it.first] = frag_inverse[it.first];
	}

	std::map<int, std::string> collapsed_barcodes;
	for (auto it : row_max) {
		if (bc_inverse[it.second.first].size() == 0) continue;
		collapsed_barcodes[it.first] = bc_inverse[it.second.first];
	}

	std::vector<FragmentStruct> collapsed;
	for (auto it : row_sum) {
		if (it.second < 1) {
			continue;
		}

		std::stringstream ss;
		ss << collapsed_frags[it.first] << "|" << collapsed_barcodes[it.first];
		FragmentStruct frag = StringToFrag(ss.str());
		frag.sum = row_sum[it.first];
		collapsed.push_back(frag);
	}
	
	delete frag_id_lookup;
	delete bc_id_lookup;

	return collapsed; 
}


// Find Complete fragments that are >max_dist bp away from the current BAM file position
// @param max_ist the maximum allowed distance between fragment start and end positions
// @param current_position the current position being looked at in the position-sorted BAM file
// @description moves completed fragments to a new dictinary, and deletes completed fragments from
// this instance fragment_dict.
FragmentMap
FragmentThread::findCompleteFragments(unsigned int current_position){
	FragmentMap completed;
	std::vector<std::string> to_erase;
	unsigned int d = this->max_distance + 20;
	for (auto item : this->fragment_dict) {
		// item.first is key
		// item.second is FragmentStruct
		if (item.second.complete) {
			// if this fragment is complete
			if (item.second.end + d < current_position) {
				completed[item.first] = item.second;
				to_erase.push_back(item.first);
			}
		} else {
			// remove incomplete fragments that are 
			// too far away to ever be complete
			if (item.second.start == -1) {
				if (item.second.end + d < current_position) {
					to_erase.push_back(item.first);
				}
			} else if (item.second.end == -1) {
				if (item.second.start + d < current_position) {
					to_erase.push_back(item.first);
				}
			} else {
				to_erase.push_back(item.first);
			}
		}
	}

	for (auto er : to_erase) {
		this->fragment_dict.erase(er);
	}


	return completed;
}

/// Collapse fragment counts that share a common start
/// or end coordinate and are from the same cell
// @param counts map of fragmentstructs as strings with a count value for each
// @param fragments FragmentMap of all fragments to build position map of fragments
// @param start bool flag if we are checking fragments at same start position, or same end position
// setting start=true includes the FragmentStruct.start field
std::map<std::string, int>
FragmentThread::collapseOverlapFragments(std::map<std::string, int> &counts, bool start) {
	std::map<std::string, int> out_counts;

	// startFrags is a map of keys (where key is "chromosome+start|end+cell_barcode")
	// and where value is a list of FragmentStructs at that position
	std::map<std::string, std::vector<std::string>> startFrags = 
		FragmentThread::createPositionLookup(counts, start);
	
	for (auto frag : startFrags) {
		// frag.first is string key
		// frag.second is forward list of strings of full fragments (fullfrags)
		int windex = -1; // keep track of the index with the highest count
		int i = 0;
		std::vector<int> countvec;
		for (auto x : frag.second) {
			int this_value = counts[x];
			countvec.push_back(this_value);

			if (windex == -1 || this_value > countvec[windex]) {
				windex = i;
			}

			i++;
		}

		std::string winner = frag.second[windex];

		out_counts[winner] = std::accumulate(countvec.begin(), countvec.end(), 0);
	}

	for (auto frag : counts) {
		if (startFrags.count(FragToString(StringToFrag(frag.first),true, start, !start, true)) == 0) {
			out_counts[frag.first] = frag.second;
		}
	}

	return out_counts;
}

/// create a map where key is subsetted fragment coords (no start or no stop)
/// value is the list of full fragments that share the coordinates in the key
/// only entries where >1 full fragments share the same coordinate are retained.
/// @param fragments is the FragmentMap containing all fragments
/// @param start indicates if start should be retained. If False, end is retained
std::map<std::string, std::vector<std::string>> 
FragmentThread::createPositionLookup(std::vector<std::string> &frags, bool start) {

	std::vector<FragmentStruct> fragsplit;
	for (auto it : frags) {
		fragsplit.push_back(StringToFrag(it));
	}
	std::vector<std::string> posfrags;
	for (auto it : fragsplit) {
		posfrags.push_back(FragToString(it, true, start, !start, true));
	}

	std::map<std::string, int> counts = 
		FragmentThread::CounterMapString(posfrags);

	std::map<std::string, std::vector<std::string>> starts;

	// iterate over all the fragments
	// for each string version which has >1 counts
	// add the full fragment to the forward list at that location.
	for (int i = 0; i < (int)frags.size(); i++) {
		if (counts[posfrags[i]] > 1) {
			starts[posfrags[i]].push_back(frags[i]);
		}
	}

	return starts;
}

// Method overloadind for createPositionLookup to
// enable passing a 'count' map instead of fragmentMap
std::map<std::string, std::vector<std::string>>
FragmentThread::createPositionLookup(std::map<std::string, int> &frags, bool start) {
	std::vector<std::string> vec_frags;
	for (auto it : frags) {
		vec_frags.push_back(it.first);
	}

	return FragmentThread::createPositionLookup(vec_frags, start);
}

/// General purpose function for generating a Counter map such that
/// each value is the number of times the given key appeared in the given map
/// Keys are generated using the provided function f.
std::map<std::string, int>
FragmentThread::CounterMapFragment(FragmentMap &fragments, std::function<std::string(FragmentStruct&, bool, bool, bool, bool)> f) {
	std::map<std::string, int> counts;
	for (auto frag : fragments) {
		counts[f(frag.second, true, true, true, true)] += 1;
	}
	return counts;
}

// Overloaded function for giving a vector of strings to save the manual conversion of fragment to string
std::map<std::string, int>
FragmentThread::CounterMapString(std::vector<std::string> &frags) {
	std::map<std::string, int> counts;
	for (auto frag : frags) {
		counts[frag]++;
	}
	return counts;
}
