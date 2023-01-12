#include <testthat.h>
#include <iostream>

#include "Trie.h"

context("Trie Data Structure") {

    test_that("Trie builds correctly & finds easy matches") {
        Trie t;

        t.Add_String("AGTC", 10, 10);

        expect_true(t.Locate_Seq_At_Pos("GGAGTCGGG", 2, 4) == 10);
        expect_true(t.Locate_Seq_At_Pos("ATTCTTGGTGT", 0, 4) == -1);
    }

    test_that("Trie can find match without guidance") {
        Trie t;

        t.Add_String("ATTC", 10, 10);
        t.Add_String("AGGT", 20, 20);

        int fp;
        expect_true(t.Locate_Seq_Subsection("AAAAATTCGGGGGGCGCGCGC", 0, 10, &fp) == 10);
        expect_true(fp == 4);
    }

	test_that("Trie can find mismatch sequences") {
		Trie t;

        t.Add_String("ATTC", 10, 10);
        t.Add_String("AGGT", 20, 20);
		t.Add_String("AGCC", 30, 30);

		// std::vector<MismatchResult> res = t.Locate_Seq_Mismatches("ATCC", 0, 4);

		// std::vector<MismatchResult> real {
		// 	{10, 2}, {30, 1}
		// };
		// expect_true(res.size() == real.size());
		// for (int i = 0; i < (int)res.size(); i++) {
		// 	expect_true(real[i].sequenceIndex == res[i].sequenceIndex);
		// 	expect_true(real[i].mismatchPosition == res[i].mismatchPosition);
		// }
	}

    test_that("Trie clears all data") {
        Trie t;

        t.Add_String("ATTC", 10, 10);
        t.Add_String("AGGT", 10, 10);
        t.Add_String("GTCG", 10, 10);

        t.Clear_Trie();
        // Trie contains no data
        expect_true(t.Locate_Seq_At_Pos("ATTCG", 0, 4) == -1);

        t.Add_String("AAAA", 10, 10);
        // Trie contains only a single barcode
        expect_true(t.Locate_Seq_At_Pos("AAAA", 0, 4) == 10);
    }
}
