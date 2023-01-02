
###########################################################
# sub-functions for report plot generation for scATAC-Seq data
###########################################################

#' @name sc_atac_plot_fragments_per_cell
#' @title A histogram of the log-number of fragments per cell
#'
#' @param sce The SingleExperimentObject produced by the sc_atac_create_sce function at the end of the pipeline
#'
#' @return returns NULL
#' @export
#'
sc_atac_plot_fragments_per_cell <- function(sce) {
    cell_stats <- as.data.frame(QC_metrics(sce))
    cell_stats$log_counts_per_cell <- log(cell_stats$counts_per_cell+1)
    
    ggplot2::ggplot(cell_stats, ggplot2::aes(x=log_counts_per_cell, y = ..count..)) +
        ggplot2::geom_histogram(color = "#E6AB02", fill = "#E6AB02", bins = 10) +
        ggplot2::stat_density(geom = "line", color = "#E6AB02") +
        ggplot2::ggtitle("Counts per cell") +
        ggplot2::xlab("log_counts_per_cell") + 
        ggplot2::ylab("count") 
}

#' @name sc_atac_plot_fragments_per_feature
#' @title A histogram of the log-number of fragments per feature
#'
#' @param sce The SingleExperimentObject produced by the sc_atac_create_sce function at the end of the pipeline
#'
#' @return returns NULL
#' @export
sc_atac_plot_fragments_per_feature <- function(sce) {
    feature_stats <- as.data.frame(feature_info(sce))
    feature_stats$log_counts_per_feature <- log(feature_stats$counts_per_feature+1)
    
    ggplot2::ggplot(feature_stats, ggplot2::aes(x=log_counts_per_feature, y = ..count..)) +
        ggplot2::geom_histogram(color = "#E6AB02", fill = "#E6AB02", bins = 10) +
        ggplot2::stat_density(geom = "line", color = "#E6AB02") +
        ggplot2::ggtitle("Counts per feature") +
        ggplot2::xlab("log_counts_per_feature") + 
        ggplot2::ylab("count") 
}


#' @name sc_atac_plot_features_per_cell
#' @title A histogram of the log-number of features per cell
#'
#' @param sce The SingleExperimentObject produced by the sc_atac_create_sce function at the end of the pipeline
#'
#' @return returns NULL
#' @export
#'
sc_atac_plot_features_per_cell <- function(sce) {
    cell_stats <- as.data.frame(QC_metrics(sce))
    cell_stats$log_features_per_cell <- log(cell_stats$features_per_cell+1)
    
    ggplot2::ggplot(cell_stats, ggplot2::aes(x=log_features_per_cell, y = ..count..)) +
        ggplot2::geom_histogram(color = "#7570B3", fill = "#7570B3", bins = 10) +
        ggplot2::stat_density(geom = "line", color = "#7570B3") +
        ggplot2::ggtitle("Features per cell") +
        ggplot2::xlab("log_features_per_cell") + 
        ggplot2::ylab("count") 
}

#' @name sc_atac_plot_features_per_cell_ordered
#' @title Plot showing the number of features per cell in ascending order
#'
#' @param sce The SingleExperimentObject produced by the sc_atac_create_sce function at the end of the pipeline
#'
#' @return returns NULL
#' @export
#'
sc_atac_plot_features_per_cell_ordered <- function(sce) {
    cell_stats <- QC_metrics(sce)
    plot(sort(cell_stats$features_per_cell), 
        xlab= 'cell', 
        log= 'y', 
        ylab = "features", 
        main= 'features per cell (ordered)',
        col = "#1B9E77")
}

#' @name sc_atac_plot_cells_per_feature
#' @title A histogram of the log-number of cells per feature
#'
#' @param sce The SingleExperimentObject produced by the sc_atac_create_sce function at the end of the pipeline
#'
#' @return returns NULL
#' @export
#'
sc_atac_plot_cells_per_feature <- function(sce) {
    feature_stats <- as.data.frame(feature_info(sce))
    feature_stats$log_cells_per_feature <- log(feature_stats$cells_per_feature+1)
    
   ggplot2::ggplot(feature_stats, ggplot2::aes(x=log_cells_per_feature, y = ..count..)) +
        ggplot2::geom_histogram(color = "#7570B3", fill = "#7570B3", bins = 10) +
        ggplot2::stat_density(geom = "line", color = "#7570B3") +
        ggplot2::ggtitle("Cells per feature") +
        ggplot2::xlab("log_cells_per_feature") + 
        ggplot2::ylab("count") 
    
}

#' @name sc_atac_plot_fragments_features_per_cell
#' @title A scatter plot of the log-number of fragments and log-number of features per cell
#'
#' @param sce The SingleExperimentObject produced by the sc_atac_create_sce function at the end of the pipeline
#'
#' @return returns NULL
#' @export
#'
sc_atac_plot_fragments_features_per_cell <- function(sce) {
    cell_stats <- as.data.frame(QC_metrics(sce))
    cell_stats$log_counts_per_cell <- log(cell_stats$counts_per_cell+1)
    cell_stats$log_features_per_cell <- log(cell_stats$features_per_cell+1)
    
    ggplot2::ggplot(cell_stats, ggplot2::aes(x=log_counts_per_cell, y=log_features_per_cell)) +
        ggplot2::geom_point(color = "#E6AB02") +
        ggplot2::xlab("log_counts_per_cell") + 
        ggplot2::ylab("log_features_per_cell") +
        ggplot2::geom_smooth(formula = y ~ x, method='lm', color = "#E6AB02", fill = "#E6AB02")
}

#' @name sc_atac_plot_fragments_cells_per_feature
#' @title A scatter plot of the log-number of fragments and log-number of cells per feature
#'
#' @param sce The SingleExperimentObject produced by the sc_atac_create_sce function at the end of the pipeline
#'
#' @return returns NULL
#' @export
#'
sc_atac_plot_fragments_cells_per_feature <- function(sce) {
    feature_stats <- as.data.frame(feature_info(sce))
    feature_stats$log_counts_per_feature <- log(feature_stats$counts_per_feature+1)
    feature_stats$log_cells_per_feature <- log(feature_stats$cells_per_feature+1)
    
    ggplot2::ggplot(feature_stats, ggplot2::aes(x=log_counts_per_feature, y=log_cells_per_feature)) +
        ggplot2::geom_point(color = "#7570B3") +
        ggplot2::xlab("log_counts_per_feature") + 
        ggplot2::ylab("log_cells_per_feature") +
        ggplot2::geom_smooth(formula = y ~ x, method='lm', color = "#7570B3", fill = "#7570B3")
}