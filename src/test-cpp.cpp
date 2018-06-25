/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.

#include "Gene.h"
#include "Interval.h"

#include "cellbarcode.h"
#include "parsecount.h"
#include "utils.h"


// ALWAYS INCLUDE TESTTHAT LAST
// Macros defined in this file can break other headers
#include <testthat.h>

context("Testing utilities") {
  test_that("Path constructor works") {
    std::string a = "aa/bb/cc";
    std::string b = "aa/bb/cc/";
    std::string c = "mm.csv";
    expect_true(join_path(a,c) == join_path(b,c));
  }

  test_that("Hamming distance works") {
    std::string a = "ATCGTAAC";
    std::string b = "ATGCTAAC";
    expect_true(hamming_distance(a,b) == 2);
  }
}

context("Testing classes") {
  Interval a(1, 5, 0);
  Interval b(2, 8, 0);
  Interval c(6, 10, 0);

  test_that("Intervals comparisons work") {
    expect_true(a == b);
    expect_true(b == c);
    expect_true(a < c);
    expect_false(a == c);
  }

  test_that("Gene comparisons work") {
    Gene ga("Gene1", 1, 5, 0);
    Gene gb("Gene2", 2, 8, 0);
    Gene gc("Gene3", 6, 10, 0);

    expect_true(ga == gb);
    expect_true(gb == gc);
    expect_true(ga < gc);
    expect_false(ga == gc);
  }

  test_that("Exon sorting and flattening works") {
    Gene g1;
    Gene g2;

    g1.add_exon(c);
    g1.add_exon(b);
    g1.add_exon(a);
    expect_true(g1.st == 1);
    expect_true(g1.en == 10);

    g1.sort_exon();
    expect_true(g1.exon_vec[0].st == a.st);
    expect_true(g1.exon_vec[0].en == a.en);
    expect_true(g1.exon_vec[1].st == b.st);
    expect_true(g1.exon_vec[1].en == b.en);
    expect_true(g1.exon_vec[2].st == c.st);
    expect_true(g1.exon_vec[2].en == c.en);

    g1.flatten_exon();

    expect_true(g1.exon_vec.size() == 1);
    expect_true(g1.exon_vec[0].st == 1);
    expect_true(g1.exon_vec[0].en == 10);
  }
}

/* 
 * Need to work out how to access internal data in C++ 
 * Can't actually use ../../inst/extdata/barcode_anno.csv because testing once
 * package is built occurs elsewhere, and / doesn't work on Windows anyway.
 */
// context("Test with toy inputs") {
//   Barcode bar;
//   std::string fn = "../../inst/extdata/barcode_anno.csv";
//   bar.read_anno(fn);

//   test_that("Barcode reading works") {
//     expect_true(bar.barcode_dict.count("ATCTGCCC") > 0);
//     expect_true(bar.barcode_dict["ACGATCGA"] == "CELL003");
//     expect_true(bar.cellid_list.size() == 10);
//     expect_true(bar.barcode_list.size() == 10);
//   }

//   test_that("Closest barcode match works") {
//     std::string bc1 = "ATGATAAA";
//     std::string bc2 = "AAAAAAAA";
//     expect_true(bar.get_closest_match(bc1, 2) == "ATGATCAT");
//     expect_true(bar.get_closest_match(bc2, 2).empty());
//   }
// }

context("Test UMI deduplication") {
  std::vector<std::pair<std::string,int>> v1;

  v1.push_back(std::make_pair("ATGCTAAC", 100));
  v1.push_back(std::make_pair("GTAGTAGC", 100));
  v1.push_back(std::make_pair("ATGCTAAC", 110));
  v1.push_back(std::make_pair("ATCTGCCC", 150));

  std::vector<std::pair<std::string,int>> v2;
  v2.push_back(std::make_pair("ATGCTAAC", 100));
  v2.push_back(std::make_pair("ATGCTAAT", 100));
  v2.push_back(std::make_pair("ATGCTAAC", 110));
  v2.push_back(std::make_pair("ATCTGCCC", 150));

  std::unordered_map<std::string, std::vector<std::pair<std::string,int>>> gene_read;
  gene_read["GENE01"] = v1;
  gene_read["GENE02"] = v1;
  gene_read["GENE03"] = v2;
  std::vector<int> UMI_dup_count(MAX_UMI_DUP + 1);
  UMI_dedup_stat s = {};
  std::unordered_map<std::string, int> tmp_res;
  tmp_res = UMI_dedup(gene_read, UMI_dup_count, s, 1, true);
  
  test_that("Genes with the same UMI are deduplicated") {
    expect_true(tmp_res["GENE01"]== 3);
    expect_true(tmp_res["GENE03"] == 2);
    expect_true(UMI_dup_count[2] == 1);
    expect_true(s.corrected_UMI == 4);
  }
}
