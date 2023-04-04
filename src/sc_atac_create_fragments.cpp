#include "sc_atac_create_fragments.h"

#include <Rcpp.h>
#include <R.h>
#include <thread>
#include <string>
#include <fstream>
#include <cstdio>
#include <memory>

#include "bam.h"
#include "FragmentThread.h"
#include "ThreadOutputFile.h"
#include "utils.h"


void merge_files(std::vector<std::string> files, std::string outname) {
    if (fileExists(outname)) {
        remove(outname.c_str());
    }
	std::ofstream outfile;
	outfile.open(outname, std::ios::out | std::ios::app);
	for (auto i : files) {
		std::ifstream infile(i);
		std::string line;
		if (infile.is_open()) {
			while (std::getline(infile, line)) {
				outfile << line << "\n";
			}
			infile.close();
		}
		if (remove(i.c_str())) {
			// unsuccessful file removal
			Rcpp::Rcout << "File " << i << " unsuccessfully deleted\n";
		}
	}
	outfile.close();
}
//' @name cpp_sc_atac_create_fragments
//' @title Generating the popular fragments for scATAC-Seq data
//' @description Takes in a tagged and sorted BAM file and outputs the associated fragments in a .bed file
//'
//' @param inbam The tagged, sorted and duplicate-free input BAM file
//' @param output The path of the output folder
//' @param contigs character vector of chromosome names from get_chromosomes()
//' @param ends integer vector of reference lengths acquired from get_chromosomes()
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
//'
//' @returns returns NULL
//'
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
		) {
			
	std::string output_filename (output + "/fragments.bed");
	std::vector<std::string> outputs;

	std::vector<std::thread> pool;
	nproc = contigs.size(); // how can we accomodate a different number of processes if user specifices?
	Rcpp::Rcout << "Output folder name is: " << output << "\n";

	for (int i = 0; i < (int)nproc; i++) {
		std::string contig = contigs[i];

		// create the bam_header_t and find the tid for this contig
		bamFile bam = bam_open(inbam.c_str(), "r"); // bam.h
		bam_header_t *header = bam_header_read(bam); // bam.h
		int tid = bam_get_tid(header, contig.c_str()); // bam.h
		bam_close(bam);

		std::string threadOut = output + "/tmp" + std::to_string(i);

    	if (fileExists(threadOut)) remove(threadOut.c_str());

		FragmentThread frag (
			threadOut,
			contig,
			tid,
			ends[i],
			inbam,
			header,
			min_mapq,
			cellbarcode,
			readname_barcode,
			cells,
			max_distance,
			min_distance,
			chunksize
		);

		outputs.push_back(threadOut);
		pool.push_back(std::thread(frag));
	}

	for (auto &th : pool) {
		th.join();
	}
	
	// need to join the outputs
	merge_files(outputs, output_filename);
	Rcpp::Rcout << "Output BED file: " << output_filename << "\n";
}
