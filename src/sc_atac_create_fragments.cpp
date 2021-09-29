#include <Rcpp.h>
#include <R.h>
#include <thread>
#include <string>
#include "bam.h"
#include "FragmentThread.hpp"

using namespace Rcpp;


//' @name sc_atac_create_fragments_cpp
//' @title Generating the popular fragments for scATAC-Seq data using sinto
//' @description Takes in a tagged and sorted BAM file and outputs the associated fragments in a .bed file
//'
//' @param inbam The tagged, sorted and duplicate-free input BAM file
//' @param output_folder The path of the output folder
//' @param contigs character vector of chromosome names from get_chromosomes()
//' @param ends integer vector of reference lengths acquired from get_chromosomes()
//' @param fragment_path : str
//'    Path for output fragment file
//' @param min_mapq : int
//'    Minimum MAPQ to retain fragment
//' @param nproc : int, optional
//'    Number of processors to use. Default is 1.
//' @param cellbarcode : str
//'   Tag used for cell barcode. Default is CB (used by cellranger)
//' @param chromosomes : str, optional
//'    Regular expression used to match chromosome names to include in the
//'    output file. Default is "(?i)^chr" (starts with "chr", case-insensitive).
//'    If None, use all chromosomes in the BAM file.
//' @param readname_barcode : str, optional
//'    Regular expression used to match cell barocde stored in read name.
//'    If None (default), use read tags instead. Use "[^:]*" to match all characters
//'    before the first colon (":").
//' @param cells : str
//'    File containing list of cell barcodes to retain. If None (default), use all cell barcodes
//'    found in the BAM file.
//' @param max_distance : int, optional
//'    Maximum distance between integration sites for the fragment to be retained.
//'    Allows filtering of implausible fragments that likely result from incorrect
//'    mapping positions. Default is 5000 bp.
//' @param min_distance : int, optional
//'    Minimum distance between integration sites for the fragment to be retained.
//'    Allows filtering implausible fragments that likely result from incorrect
//'    mapping positions. Default is 10 bp.
//' @param chunksize : int
//'    Number of BAM entries to read through before collapsing and writing
//'    fragments to disk. Higher chunksize will use more memory but will be
//'    faster.
//' @import Rhtslib
//' @import Rcpp
//' @useDynLib scPipeFragments, .registration=TRUE
//' @export
// [[Rcpp::export]]
void sc_atac_create_fragments_cpp(
		std::string inbam,
		std::string output,
		CharacterVector contigs,
		IntegerVector ends,
		unsigned int min_mapq,
		unsigned int nproc,
		std::string cellbarcode,
		std::string chromosomes,
		String readname_barcode,
		StringVector cells,
		unsigned int max_distance,
		unsigned int min_distance,
		unsigned int chunksize
		) {

	Rcout << "Inside CPP function\n";

	//pool of threads of length nproc to use for fragments
	std::vector<std::string> outnames;
	std::vector<std::thread> pool;

	//temp TODO: Figure output how to use worker pool
	nproc = contigs.length();


	Rcpp::Rcout << "Output file name is: " << output << "\n";
	int i = 0;
	for (int i = 0; i < nproc; i++) {
		// setup fragment thread objects
		std::string outname = output + "/tempFragmentFile" + std::to_string(i);

		outnames.push_back(output);
		std::string contig = String(contigs[i]).get_cstring();

		// create the bam_header_t and find the tid for this contig
		bamFile bam = bam_open(inbam.c_str(), "r"); // bam.h
		bam_header_t *header = bam_header_read(bam); // bam.h
		int tid = bam_get_tid(header, contig.c_str()); // bam.h
		bam_close(bam);

		FragmentThread frag (
			output,
			String(contigs[i]).get_cstring(),
			tid,
			ends[i],
			inbam,
			header,
			min_mapq,
			cellbarcode,
			readname_barcode.get_cstring(),
			cells,
			max_distance,
			min_distance,
			chunksize
		);

		pool.push_back(std::thread(frag));
	}


	for (auto &th : pool) {
		th.join();
	}

	// now we can access the files

	Rcout << "Finished threading\n";
}
