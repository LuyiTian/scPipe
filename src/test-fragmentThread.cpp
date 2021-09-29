#include <testthat.h>
#include "FragmentThread.hpp"
#include "FragmentUtils.hpp"

context("Fragment Tests") {
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
}