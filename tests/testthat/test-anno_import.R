context("Annotation import")
test_that("ENSEMBL gtf annotation can be imported", {
  anno_path <- system.file("extdata", "ensembl_y_tiny.gtf.gz", package = "scPipe")
  
  expected <- data.frame(
    GeneID = c("ENSG00000251841", "ENSG00000184895", "ENSG00000237659", "ENSG00000232195"),
    Chr = factor(c("Y", "Y", "Y", "Y")),
    Start = as.integer(c(2784749, 2786855, 2789827, 2827982)),
    End = as.integer(c(2784853, 2787699, 2790328, 2828218)),
    Strand = factor(c("+", "-", "+", "+"), levels = c("+", "-", "*")),
    stringsAsFactors = FALSE
  )
  
  expect_identical(anno_import(anno_path), expected)
})

test_that("Common annotation sources can be imported", {
    expect_silent(anno_import(system.file("extdata", "ens_tiny_anno.gff3.gz", package = "scPipe")))
    expect_silent(anno_import(system.file("extdata", "ens_tiny_anno.gtf.gz", package = "scPipe")))
    expect_equal(
        anno_import(system.file("extdata", "ens_tiny_anno.gff3.gz", package = "scPipe")),
        anno_import(system.file("extdata", "ens_tiny_anno.gtf.gz", package = "scPipe"))
    )
    expect_silent(anno_import(system.file("extdata", "gen_tiny_anno.gff3.gz", package = "scPipe")))
    expect_silent(anno_import(system.file("extdata", "gen_tiny_anno.gtf.gz", package = "scPipe")))
    expect_equal(
        anno_import(system.file("extdata", "gen_tiny_anno.gff3.gz", package = "scPipe")),
        anno_import(system.file("extdata", "gen_tiny_anno.gtf.gz", package = "scPipe"))
    )
    expect_silent(anno_import(system.file("extdata", "ref_tiny_anno.gff3.gz", package = "scPipe")))
})

test_that("Errors are properly reported", {
  anno <- rtracklayer::import(system.file("extdata", "ensembl_y_tiny.gtf.gz", package = "scPipe"))
  
  anno_missing_gene_id <- anno
  GenomicRanges::mcols(anno_missing_gene_id) <- subset(GenomicRanges::mcols(anno_missing_gene_id), select = -gene_id) # delete gene_id column
  expect_error(anno_to_saf(anno_missing_gene_id), "'gene_id' column missing from GRanges metadata and could not be inferred")
})

test_that("SAF validation works", {
  expect_error(validate_saf(1), "annotation must object of class data.frame")
  
  anno <- anno_import(system.file("extdata", "ens_tiny_anno.gff3.gz", package = "scPipe"))
  expect_error(validate_saf(anno[, 5:1]), "columns of SAF data.frame must be: GeneID, Chr, Start, End, Strand")
  
  anno_na <- anno
  anno_na[4, 1] <- NA
  expect_error(validate_saf(anno_na), "SAF data.frame must not contain any NA")
})
