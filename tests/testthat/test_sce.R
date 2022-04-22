context("scPipe testing")

test_that("new SingleCellExperiment does not contain QC information", {
  data("sc_sample_data")
  data("sc_sample_qc")
  sce = SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
  expect_that(sce, is_a("SingleCellExperiment"))
  expect_warning(val <- QC_metrics(sce),
                 "`scPipe` not in `metadata`. Cannot identify quality control columns")
  expect_null(val)

  expect_warning(val <- demultiplex_info(sce),
                 "`scPipe` not in `metadata`. Cannot find columns for cell barcode demultiplex results")
  expect_null(val)

  expect_warning(val <- UMI_dup_info(sce),
                 "`scPipe` not in `metadata`. Cannot find columns for cell barcode demultiplex results")
  expect_null(val)

})

test_that("new SingleCellExperiment should work properly", {
  data("sc_sample_data")
  data("sc_sample_qc")
  sce = SingleCellExperiment(assays = list(counts = as.matrix(sc_sample_data)))
  QC_metrics(sce) = sc_sample_qc
  demultiplex_info(sce) = cell_barcode_matching
  UMI_dup_info(sce) = UMI_duplication
  expect_equal(dim(sce), c(1000, 383))
  expect_equal(dim(QC_metrics(sce)), c(383, 12))
  expect_equal(dim(UMI_dup_info(sce)), c(1001, 2))
  expect_equal(dim(demultiplex_info(sce)), c(6, 2))
})
