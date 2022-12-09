context("Check Barcode Reads before Trimming")
test_that("Correct input positions are accepted", {
    small_r1 <- system.file("extdata", "simu_R1.fastq.gz", package="scPipe")
    small_barcode <- system.file("extdata", "smallbarcode1col.csv", package="scPipe")

    expect_true(check_barcode_start_position(small_r1, small_barcode, "smallbarcode1col.csv", 0, 4, 40, .5))
})

test_that("Bad input position causes rejection", {
    small_r1 <- system.file("extdata", "simu_R1.fastq.gz", package="scPipe")
    small_barcode <- system.file("extdata", "smallbarcode1col.csv", package="scPipe")

    expect_false(check_barcode_start_position(small_r1, small_barcode, "smallbarcode1col.csv", 10, 4, 40, .5))
})

test_that("Bad match rate causes rejection", {
    small_r1 <- system.file("extdata", "simu_R1.fastq.gz", package="scPipe")
    small_barcode <- system.file("extdata", "smallbarcode1col.csv", package="scPipe")

    expect_false(check_barcode_start_position(small_r1, small_barcode, "smallbarcode1col.csv", 0, 4, 100, .5))
})
