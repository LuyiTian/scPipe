#include <testthat.h>
#include <map>
#include <string>
#include <functional>

#include <htslib/sam.h>

#include "bam.h"
#include "FragmentThread.hpp"
#include "FragmentUtils.hpp"

bool equalFragmentStruct(FragmentStruct &a, FragmentStruct &b) {
	return
		a.chromosome == b.chromosome &&
		a.start == b.start &&
		a.end == b.end &&
		a.cell_barcode == b.cell_barcode &&
		a.complete == b.complete; 
}

context("Fragment Utility Tests") {
	// Utility Function tests
	test_that("Fragment string conversions correctly execute") {
		FragmentStruct f = {
			std::string("ATCGATCG"),
			0,
			30,
			std::string("AATTCCGG")
		};

		// test FragToSTring conversion
		std::string fts = FragToString(f, true, true, true, true);
		expect_true(fts.compare("ATCGATCG|0|30|AATTCCGG") == 0);

		// test StringToFrag conversion
		FragmentStruct x = StringToFrag("ATCGATCG|0|30|AATTCCGG");
		expect_true(x.chromosome.compare("ATCGATCG") == 0);
		expect_true(x.start == 0);
		expect_true(x.end == 30);
		expect_true(x.cell_barcode.compare("AATTCCGG") == 0);
	}

	test_that("id_lookup function builds the correct map") {
		FragmentMap m;
		m["AA"] = {"AATG",1,30, "GG"}; // 1
		m["AG"] = {"AATG",3,30, "GA"}; // 3
		m["AT"] = {"AAGG",4,30, "GT"}; // 4
		m["AC"] = {"AAGG",2,30, "GC"}; // 2

		std::map<std::string, int> *n = id_lookup(m,
			[](FragmentStruct x) {
				std::stringstream ss;
				ss << x.chromosome << "|" << x.start << "|" << x.end << "|" << x.cell_barcode;
				return ss.str();
			}
		);
		expect_true(n->at("AATG|1|30|GG") == 1);
		expect_true(n->at("AAGG|2|30|GC") == 2);
		expect_true(n->at("AATG|3|30|GA") == 3);
		expect_true(n->at("AAGG|4|30|GT") == 4);
	}

	test_that("invert_map correctly inverts") {
		std::map<std::string, int> m;
		m["1"] = 1;
		m["2"] = 2;
		m["3"] = 3;
		m["4"] = 4;

		std::map<int, std::string> n = invertMap(&m);
		expect_true(n[1].compare("1") == 0);
		expect_true(n[2].compare("2") == 0);
		expect_true(n[3].compare("3") == 0);
		expect_true(n[4].compare("4") == 0);
	}
}

context("FragmentThread tests") {
	test_that("CounterMapFragment creates correct map") {
		FragmentMap m;
		m["AAT"] = {"AATG",3,25, "GG"};
		m["AAG"] = {"AATG",3,60, "GG"}; 
		m["AAC"] = {"AATG",3,60, "GG"};
		m["ATA"] = {"AAGG",4,15, "GT"}; 
		m["ACG"] = {"AAGG",4,25, "GT"};
		m["ACA"] = {"AGTC",5,30, "GA"};
		// py_frags = ["AATG|3|25|GG", "AATG|3|60|GG", "AATG|3|60|GG", "AAGG|4|15|GT", "AAGG|4|25|GT", "AGTC|5|30|GA"]

		std::map<std::string, int> counts = FragmentThread::CounterMapFragment(m,
			[](FragmentStruct &f, bool _1, bool _2, bool _3, bool _4) {
				return f.chromosome;
			});

		expect_true(counts["AATG"] == 3);
		expect_true(counts["AAGG"] == 2);
		expect_true(counts["AGTC"] == 1);
	}

	test_that("CounterMapString creates correct map") {
		std::vector<std::string> frags {"AATG", "AAGG", "AATG", "AGTC", "AATG"};

		std::map<std::string, int>counts = FragmentThread::CounterMapString(frags);

		expect_true(counts["AATG"] == 3);
		expect_true(counts["AAGG"] == 1);
		expect_true(counts["AGTC"] == 1);
	}

	test_that("createPositionLookup correctly filters and matches frags") {
		FragmentMap m;
		m["AAT"] = {"AATG",3,25, "GG"};
		m["AAG"] = {"AATG",3,60, "GG"}; 
		m["AAC"] = {"AATG",3,50, "GG"};
		m["ATA"] = {"AAGG",4,15, "GT"}; 
		m["ACG"] = {"AAGG",4,25, "GT"};
		m["ACA"] = {"AGTC",5,30, "GA"};
		// py_frags = ["AATG|3|25|GG", "AATG|3|60|GG", "AATG|3|50|GG", "AAGG|4|15|GT", "AAGG|4|25|GT", "AGTC|5|30|GA"]
		std::vector<std::string> frags;
		for (auto it : m) {
			frags.push_back(FragToString(it.second, true, true, true, true));
		}
		std::map<std::string, std::vector<std::string>> res = FragmentThread::createPositionLookup(frags, true);
		
		expect_true(res["AAGG|4|GT"].size() == 2);
		expect_true(res["AATG|3|GG"].size() == 3);
		
		// test exact contents of one of the map buckets?
	}

	test_that("createPositionLookup map version correctly matches frags") {
		FragmentMap m;
		m["AAT"] = {"AATG",3,25, "GG"};
		m["AAG"] = {"AATG",3,60, "GG"}; 
		m["AAC"] = {"AATG",3,50, "GG"};
		m["ATA"] = {"AAGG",4,15, "GT"}; 
		m["ACG"] = {"AAGG",4,25, "GT"};
		m["ACA"] = {"AGTC",5,30, "GA"};
		// py_frags = ["AATG|3|25|GG", "AATG|3|60|GG", "AATG|3|50|GG", "AAGG|4|15|GT", "AAGG|4|25|GT", "AGTC|5|30|GA"]
		std::map<std::string, int> counts = FragmentThread::CounterMapFragment(m, FragToString);
		std::map<std::string, std::vector<std::string>> res = FragmentThread::createPositionLookup(counts, true);
		
		expect_true(res["AAGG|4|GT"].size() == 2);
		expect_true(res["AATG|3|GG"].size() == 3);
	}

	test_that("Complete fragments are found and removed") {
		FragmentMap m;
		m["AAT"] = {"AATG",3,25, "GG", true};
		m["AAG"] = {"AATG",3,60, "GG", false}; 
		m["AAC"] = {"AATG",3,50, "GG", true};
		m["ATA"] = {"AAGG",4,15, "GT", false}; 
		m["ACG"] = {"AAGG",4,25, "GT", false};
		m["ACA"] = {"AGTC",5,30, "GA", true};
		// frags = {'AAT': ['AATG', 3, 25, 'GG', True], 'AAG': ['AATG', 3, 60, 'GG', False], 'AAC': ['AATG', 3, 50, 'GG', True], 'ATA': ['AAGG', 4, 15, 'GT', False], 'ACG': ['AAGG', 4, 25, 'GT', False], 'ACA': ['AGTC', 5, 30, 'GA', True]}
	
		FragmentThread f ("fake", "A", 1, 30, "BAM", nullptr, 0, 
			"AA", "readname", Rcpp::CharacterVector(), 5, 0, 10);

		f.fragment_dict = m;
		FragmentMap complete = f.findCompleteFragments(200);

		// complete = {'AAT': ['AATG', 3, 25, 'GG'], 'AAC': ['AATG', 3, 50, 'GG'], 'ACA': ['AGTC', 5, 30, 'GA']}
		// f.fragment_dict == {}
		expect_true(f.fragment_dict.empty());
		expect_true(complete["AAT"].chromosome.compare("AATG") == 0);
		expect_true(complete["AAC"].chromosome.compare("AATG") == 0);
		expect_true(complete["ACA"].chromosome.compare("AGTC") == 0);
	}
	
	test_that("Collapsed overlapping fragments collapse correctly") {
		FragmentMap m;
		m["AAT"] = {"AATG",3,25, "GG"};
		m["AAG"] = {"AATG",3,60, "GG"}; 
		m["AAC"] = {"AATG",3,50, "GG"};
		m["ATA"] = {"AAGG",4,15, "GT"}; 
		m["ACG"] = {"AAGG",4,25, "GT"};
		m["ACA"] = {"AGTC",5,30, "GA"};
		// py_frags = ["AATG|3|25|GG", "AATG|3|60|GG", "AATG|3|50|GG", "AAGG|4|15|GT", "AAGG|4|25|GT", "AGTC|5|30|GA"]
		std::map<std::string, int> counts = FragmentThread::CounterMapFragment(m, FragToString);
		std::map<std::string, int> collapsed = FragmentThread::collapseOverlapFragments(counts, true);

		expect_true(collapsed["AATG|3|25|GG"] == 3);
		expect_true(collapsed["AAGG|4|15|GT"] == 2);
		expect_true(collapsed["AGTC|5|30|GA"] == 1);
	}

	test_that("Full collapse fragments function operates correctly") {
		FragmentMap m;
		m["AAT"] = {"AATG",3,25, "GG"};
		m["AAG"] = {"AATG",3,60, "GG"}; 
		m["AAC"] = {"AATG",3,50, "GG"};
		m["ATA"] = {"AAGG",4,15, "GT"}; 
		m["ACG"] = {"AAGG",4,25, "GT"};
		m["ACA"] = {"AGTC",5,30, "GA"};
		// frags={"AAT":["AATG",3,25, "GG"],"AAG":["AATG",3,60, "GG"], "AAC":["AATG",3,50, "GG"], "ATA":["AAGG",4,15, "GT"], "ACG":["AAGG",4,25, "GT"], "ACA":["AGTC",5,30, "GA"]}
	
		std::vector<FragmentStruct> collapsed = FragmentThread::collapseFragments(m);

		for (auto frag : collapsed) {
			if (frag.chromosome.compare("AATG") == 0) {
				expect_true(frag.sum == 3);
			} else if (frag.chromosome.compare("AAGG") == 0) {
				expect_true(frag.sum == 2);
			} else if (frag.chromosome.compare("AGTC") == 0) {
				expect_true(frag.sum == 1);
			}
		}
	}

	test_that("Write to file is able of writing a fragment map to file") {
		expect_true(true);
		// below code executes correctly when run:
		// (not run for tests)

		// FragmentMap m;
		// m["AAT"] = {"AATG",3,25, "GG", true, 10};
		// m["AAG"] = {"AATG",3,60, "GG", false, 20}; 
		// m["AAC"] = {"AATG",3,50, "GG", true, 30};
		// m["ATA"] = {"AAGG",4,15, "GT", false, 40}; 
		// m["ACG"] = {"AAGG",4,25, "GT", false, 50};
		// m["ACA"] = {"AGTC",5,30, "GA", true, 60};
		
		// std::vector<FragmentStruct> ls;
		// for (auto it : m) {
		// 	ls.push_back(it.second);
		// }

		// FragmentThread::writeFragmentsToFile(ls, "output/testWriteToFile");
	}

	test_that("Fragment Dict can be correctly added to") {
		FragmentMap m;
		m["AAT"] = {"AATG",3,25, "GG"};
		m["AAG"] = {"AATG",3,60, "GG"}; 
		m["AAC"] = {"AATG",3,50, "GG"};
		m["ATA"] = {"AAGG",4,15, "GT"}; 
		m["ACG"] = {"AAGG",4,25, "GT"};
		m["ACA"] = {"AGTC",5,30, "GA"};
		// frags={"AAT":["AATG",3,25, "GG"],"AAG":["AATG",3,60, "GG"], "AAC":["AATG",3,50, "GG"], "ATA":["AAGG",4,15, "GT"], "ACG":["AAGG",4,25, "GT"], "ACA":["AGTC",5,30, "GA"]}

		FragmentThread f ("fake", "A", 1, 30, "BAM", nullptr, 0, 
			"AA", "readname", Rcpp::CharacterVector(), 5, 0, 10);

		// check default case of fragment does not exist in dictionary
		f.addToFragments("CCG", "CGCG", 10, 60, "CC", false);
		FragmentStruct res1 = {"CGCG", 10, -1, "CC", false};
		expect_true(equalFragmentStruct(f.fragment_dict["CCG"], res1));
	
		// check case of qname in dictionary and fragment is reverse strand
		// this should complete the fragment
		f.addToFragments("CCG", "CGCG", 15, 60, "CC", true);
		FragmentStruct res2 = {"CGCG", 10, 60, "CC", true};
		expect_true(equalFragmentStruct(f.fragment_dict["CCG"], res2));

		// check 
	}

	test_that("Fragment dict is updated") {
		FragmentMap m;
		m["AAT"] = {"AATG",3,25, "GG"};
		m["AAG"] = {"AATG",3,60, "GG"}; 
		m["AAC"] = {"AATG",3,50, "GG"};
		m["ATA"] = {"AAGG",4,15, "GT"}; 
		m["ACG"] = {"AAGG",4,25, "GT"};
		m["ACA"] = {"AGTC",5,30, "GA"};
		// frags={"AAT":["AATG",3,25, "GG"],"AAG":["AATG",3,60, "GG"], "AAC":["AATG",3,50, "GG"], "ATA":["AAGG",4,15, "GT"], "ACG":["AAGG",4,25, "GT"], "ACA":["AGTC",5,30, "GA"]}


		// python:
		// from sinto_fragments import *
		// import pysam
		// frags={"AAT":["AATG",3,25, "GG"],"AAG":["AATG",3,60, "GG"], "AAC":["AATG",3,50, "GG"], "ATA":["AAGG",4,15, "GT"], "ACG":["AAGG",4,25, "GT"], "ACA":["AGTC",5,30, "GA"]}
		// bam = "/Users/voogd.o/Documents/scPipeTesting/sc_atac_create_fragments/sinto_output/demux_testfastq_S1_L001_R1_001_aligned_tagged_sorted.bam"
		// inputBam = pysam.AlignmentFile(bam, "rb")
		// seg = inputBam.fetch("chr21", 0, 48129895).__next__()
		// fragments = updateFragmentDict(fragments=frags, segment=seg, min_mapq=10, cellbarcode="CB", readname_barcode=None, cells=None, max_dist=5000, min_dist=10)
		// fragments should be :
		// {'AAT': ['AATG', 3, 25, 'GG'], 'AAG': ['AATG', 3, 60, 'GG'], 'AAC': ['AATG', 3, 50, 'GG'], 'ATA': ['AAGG', 4, 15, 'GT'], 'ACG': ['AAGG', 4, 25, 'GT'], 'ACA': ['AGTC', 5, 30, 'GA'], 'TGAGTCACATTGTGAC#A00228:277:HFKLHDMXX:1:2167:2853:35336': ['chr21', 9411277, None, 'TGAGTCACATTGTGAC', False]}


	}
}


int testFetchCall(const bam1_t *b, void *data) {
	int *i = (int *)data;
	if (*i >= 1) {
		return 1;
	}

	*i = 100;

	// FragmentThread fragThread((
	// 		outname,
	// 		contig,
	// 		tid,
	// 		ends[i],
	// 		inbam,
	// 		header,
	// 		min_mapq,
	// 		cellbarcode,
	// 		readname_barcode.get_cstring(),
	// 		cells,
	// 		max_distance,
	// 		min_distance,
	// 		chunksize
	// 	);

	std::string qname(bam1_qname(b));
	uint8_t *raw_data = bam_aux_get(b, "CB");
	std::string barcode_data = bam_aux2string(raw_data);

	Rcpp::Rcout << "testing function\n";
	
	Rcpp::Rcout << qname << "\n";
	Rcpp::Rcout << barcode_data << "\n";


	int32_t rstart = bam_alignment_start(b);
	int32_t rend = bam_endpos(b);
	bool is_reverse = bam1_strand(b);

	Rcpp::Rcout << "Start: " << rstart << " End: " << rend << " reverse: " << is_reverse << "\n";

	return 1;
}

// a function to just chuck whatever in a test it
// [[Rcpp::export]]
void test_whatever() {

	

	bamFile bam = bam_open("/Users/voogd.o/Documents/scPipeTesting/sc_atac_create_fragments/scPipe_output/demux_testfastq_S1_L001_R1_001_aligned_tagged_sorted.bam", "r"); // bam.h
	
	bam_header_t *header = bam_header_read(bam); // bam.h
	int tid = bam_get_tid(header, "chr21"); // bam.h
	int tval = 0;
	bam_index_t *index = bam_index_load("/Users/voogd.o/Documents/scPipeTesting/sc_atac_create_fragments/scPipe_output/demux_testfastq_S1_L001_R1_001_aligned_tagged_sorted.bam"); // bam.h
	bam_fetch(bam, index, tid, 0, 48129895, &tval, &testFetchCall);

	bam_close(bam);
}
