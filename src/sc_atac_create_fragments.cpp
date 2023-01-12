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

// void merge_files(std::vector<std::string> files, std::string outname) {
// 	std::ofstream outfile;
// 	outfile.open(outname, std::ios::out | std::ios::app);
// 	for (auto i : files) {
// 		std::ifstream infile(i);
// 		std::string line;
// 		if (infile.is_open()) {
// 			while (std::getline(infile, line)) {
// 				outfile << line << "\n";
// 			}
// 			infile.close();
// 		}
// 		if (remove(i.c_str())) {
// 			// unsuccessful file removal
// 			Rcout << "File " << i << " unsuccessfully deleted\n";
// 		}
// 	}
// 	outfile.close();
// }

//' @name sc_atac_create_fragments_cpp
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
//' @param readname_barcodeN : str, optional
//'    Regular expression used to match cell barocde stored in read name.
//'    If None (default), use read tags instead. Use "[^:]*" to match all characters
//'    before the first colon (":").
//' @param cellsN : str
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
//' @useDynLib scPipe, .registration = TRUE
// [[Rcpp::export]]
void sc_atac_create_fragments_cpp(
		std::string inbam,
		std::string output,
		Rcpp::CharacterVector contigs,
		Rcpp::IntegerVector ends,
		unsigned int min_mapq,
		unsigned int nproc,
		std::string cellbarcode,
		std::string chromosomes,
		Rcpp::Nullable<Rcpp::String> readname_barcodeN,
		Rcpp::Nullable<Rcpp::StringVector> cellsN,
		unsigned int max_distance,
		unsigned int min_distance,
		unsigned int chunksize
		) {

	Rcpp::StringVector cells = cellsN.isNotNull() ? Rcpp::StringVector(cellsN) : Rcpp::StringVector(0);

	Rcpp::String readname_barcode = readname_barcodeN.isNotNull() ? Rcpp::String(readname_barcodeN) : Rcpp::String();

	std::string output_filename (output + "/fragments.bed");
	std::shared_ptr<ThreadOutputFile> threadOutputFile = std::make_shared<ThreadOutputFile>(output_filename);

	std::vector<std::thread> pool;
	nproc = contigs.length();
	Rcpp::Rcout << "Output folder name is: " << output << "\n";
	for (int i = 0; i < (int)nproc; i++) {
		std::string contig = Rcpp::String(contigs[i]).get_cstring();

		// create the bam_header_t and find the tid for this contig
		bamFile bam = bam_open(inbam.c_str(), "r"); // bam.h
		bam_header_t *header = bam_header_read(bam); // bam.h
		int tid = bam_get_tid(header, contig.c_str()); // bam.h
		bam_close(bam);

		FragmentThread frag (
			threadOutputFile,
			contig,
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

	Rcpp::Rcout << "Output BED file: " << output_filename << "\n";
}
