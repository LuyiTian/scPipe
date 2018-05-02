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
#include <testthat.h>
#include "utils.h"
#include "Interval.h"
#include "Gene.h"
#include "cellbarcode.h"

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
  test_that("Intervals comparisons work") {
    Interval a(1, 5, 0);
    Interval b(2, 8, 0);
    Interval c(6, 10, 0);

    expect_true(a == b);
    expect_true(b == c);
    expect_true(a < c);
    expect_false(a == c);
  }

  test_that("Gene comparisons work") {
    Gene a("Gene1", 1, 5, 0);
    Gene b("Gene2", 2, 8, 0);
    Gene c("Gene3", 6, 10, 0);

    expect_true(a == b);
    expect_true(b == c);
    expect_true(a < c);
    expect_false(a == c);
  }

  test_that("Exon sorting and flattening works") {
    Gene g1;
    Gene g2;
    Interval a(1, 5, 0);
    Interval b(2, 8, 0);
    Interval c(6, 10, 0);

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

context("Test with toy inputs") {
  test_that("Barcode reading works") {
    Barcode bar;
    std::string fn = "./inst/extdata/barcode_anno.csv";
    bar.read_anno(fn);

    expect_true(bar.barcode_dict.count("ATCTGCCC") > 0);
    expect_true(bar.barcode_dict["ACGATCGA"] == "CELL003");
    expect_true(bar.cellid_list.size() == 3);
    expect_true(bar.barcode_list.size() == 4);
  }
}
