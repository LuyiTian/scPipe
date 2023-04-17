#ifndef SC_ATAC_CREATE_FRAGMENTS
#define SC_ATAC_CREATE_FRAGMENTS

#include <string>
#include <vector>

#include <Rcpp.h>
#include <R.h>

void cpp_sc_atac_create_fragments(
		std::string inbam,
		std::string output,
		std::vector<std::string> contigs,
	    std::vector<int> ends,
		unsigned int min_mapq,
		unsigned int nproc,
		std::string cellbarcode,
		std::string chromosomes,
		std::string readname_barcode,
		std::vector<std::string> cells,
		unsigned int max_distance,
		unsigned int min_distance,
		unsigned int chunksize
		);

#endif // SC_ATAC_CREATE_FRAGMENTS