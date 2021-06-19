#include <testthat.h>
#include "ResizeArray.h"

context("Resize Array Tests") {

    // The format for specifying tests is similar to that of
    // testthat's R functions. Use 'test_that()' to define a
    // unit test, and use 'expect_true()' and 'expect_false()'
    // to test the desired conditions.
    test_that("Correct values in array") {
        ResizeArray ra(3);

        expect_true(ra[0] == 0);
        expect_true(ra[1] == 0);
        expect_true(ra[2] == 0);
    }

    test_that("Increment operation works") {
        ResizeArray ra(10);

        for (int i = 0; i < 10; i++) {
            for (int z = 0; z < 100; z++) {
                ra.Increment(i);
            }
        }

        expect_true(ra[1] == 100);
        expect_true(ra[5] == 100);
    }

    test_that("Resizing is successful") {
        ResizeArray ra(2);

        expect_true(ra.length() == 2);

        ra.Increment(3);

        expect_true(ra.length() == 4);

        ra.Increment(5);

        expect_true(ra.length() == 8);
    }

    test_that("Max function finds correct maximum") {
        ResizeArray ra(10);

        for (int i = 0; i < 100; i++) {
            ra.Increment(5);
        }

        for (int i = 0; i < 40; i++) {
            ra.Increment(8);
        }

        long value;
        int pos = ra.Max(&value);

        expect_true(pos == 5);
        expect_true(value == 100);
    }

}