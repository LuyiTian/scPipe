#include <testthat.h>

#include "Trie.h"

context("Trie Data Structure") {

    test_that("Trie builds correctly & finds easy matches") {
        Trie t;

        t.Add_String("ABCD", 10, 10);

        //expect_true(t.Locate_Seq_At_Pos("ABCDEEEEE", 0, 4) == 10);
        expect_true(t.Locate_Seq_At_Pos("ABEEEEE", 0, 4) == -1);
    }

    test_that("Trie can find match without guidance") {
        Trie t;

        t.Add_String("ATTC", 10, 10);
        t.Add_String("AGGT", 20, 20);

        int fp;
        //expect_true(t.Locate_Seq_Subsection("EATTC", 0, 6, &fp) == 10);
        //expect_true(fp == 1);
    }
}
