library(SingleCellExperiment)

library(scPipe)

data("sc_sample_data")
data("sc_sample_qc")

sce = SingleCellExperiment(assays = list(counts =as.matrix(sc_sample_data)))


#sce@int_metadata[["Biomart"]] = list(organism="mmusculus_gene_ensembl",
#                                     gene_id_type = "external_gene_name")


sce_qc = calculate_QC_metrics(sce)

sce_qc_o = detect_outlier(sce_qc)

head(QC_metrics(sce_qc_o))


plotQC_pair(sce_qc_o)

