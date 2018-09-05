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

    expect_silent(anno_import(system.file("extdata", "ens_tiny_anno.gff3.gz", package = "scPipe")))
    expect_silent(anno_import(system.file("extdata", "ens_tiny_anno.gtf.gz", package = "scPipe")))
    expect_silent(anno_import(system.file("extdata", "gen_tiny_anno.gff3.gz", package = "scPipe")))
    expect_silent(anno_import(system.file("extdata", "gen_tiny_anno.gtf.gz", package = "scPipe")))
    expect_silent(anno_import(system.file("extdata", "ref_tiny_anno.gff3.gz", package = "scPipe")))
})

test_that("Errors are properly reported", {
    anno <- rtracklayer::import(system.file("extdata", "ensembl_y_tiny.gtf.gz", package = "scPipe"))

    anno_missing_gene_id <- anno
    mcols(anno_missing_gene_id) <- subset(mcols(anno_missing_gene_id), select = -gene_id)
    expect_error(anno_to_saf(anno_missing_gene_id), "'gene_id' column missing from GRanges metadata and could not be inferred")
})
