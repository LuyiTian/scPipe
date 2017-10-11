#' @name sc_sample_data
#' @title a small sample scRNA-seq counts dataset to 
#' demonstrate capabilities of scPipe
#' @description This data set contains counts for high variable genes for 
#' 100 cells. The cells have different cell types. The data contains
#' raw read counts. The cells are chosen randomly from 384 cells and
#' they did not go through quality controls. The rows names are 
#' Ensembl gene ids and the columns are cell names, which is the wall
#' position in the 384 plates.
#' @return NULL, but makes a matrix of count data
#' @docType data
#' @usage sc_sample_data
#' @format a matrix instance, one row per gene.
#' @source Christin Biben (WEHI). She FACS sorted cells from several immune
#' cell types including B cells, granulocyte and some early progenitors. 
#' @author Luyi Tian
NULL