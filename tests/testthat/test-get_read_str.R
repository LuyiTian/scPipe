context("Read structure helper")
test_that("Stripping and lower casing is correct", {
    expect_equal(
        get_read_str("celSeq"),
        list(bs1 = -1, bl1 = 0, bs2 = 6, bl2 = 8,  us = 0,  ul = 6)
    )

    expect_equal(
        get_read_str("cel-Seq"),
        list(bs1 = -1, bl1 = 0, bs2 = 6, bl2 = 8,  us = 0,  ul = 6)
    )

    expect_equal(
        get_read_str("CelSeq"),
        list(bs1 = -1, bl1 = 0, bs2 = 6, bl2 = 8,  us = 0,  ul = 6)
    )
})

test_that("Unknown protocols are properly caught", {
    expect_error(get_read_str("FakeSeq"))
    expect_error(get_read_str("errorseq"))
})
