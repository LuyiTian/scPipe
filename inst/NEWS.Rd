\name{NEWS}
\title{News for Package \pkg{scPipe}}

\section{Changes in version 0.99.0 (2017-07-28)}{
  \itemize{
    \item Package prepared for Bioconductor submission.
  }
}

\section{Changes in version 0.99.20 (2017-09-22)}{
  \itemize{
    \item scPipe now supports SingleCellExperiment class and use it as the base class
    \item add two functions `plot_demultiplex` and `plot_UMI_dup`
    \item scPipe support the offical bam tags for cell barcode and UMI
  }
}

\section{Changes in version 1.0.1 (2017-11-28)}{
  \itemize{
    \item Bug Fix: Fixed handling of colData through QC_metrics (https://github.com/LuyiTian/scPipe/issues/34)
  }
}

\section{Changes in version 1.0.2 (2017-12-01)}{
  \itemize{
    \item fix errors in unittest
  }
}

\section{Changes in version 1.0.3 (2017-12-03)}{
  \itemize{
    \item fix a bug in `validObject`. the default value for gene id and organism is set to NA
  }
}

\section{Changes in version 1.0.4 (2017-12-04)}{
  \itemize{
    \item In `detect_outlier`, give more informative error message when some cells or QC metrics have zero values.
  }
}

\section{Changes in version 1.0.5 (2017-12-14)}{
  \itemize{
    \item Fixed bugs in slim report and trimbarcode error message
    \item Fix incomplete error message
    \item Documentation updates, new functions and bug fixes
  }
}

\section{Changes in version 1.0.6 (2017-12-18)}{
  \itemize{
    \item Now the id conversion can also be done by using the bioconductor annotation package, when biomart fails to connect.
  }
}

\section{Changes in version 1.0.8 (2018-2-5)}{
  \itemize{
    \item update the gene id annotation code
  }
}

\section{Changes in version 1.0.9 (2018-4-25)}{
  \itemize{
    \item fix the bug that might misalign the exon mapping reads
    \item support gff files from gencode or refseq
  }
}

\section{Changes in version 1.1.0-1.1.3 (2018-4-29)}{
  \itemize{
    \item fix the bugs that related to the reimplementation of the exon mapping algorithm
  }
}

\section{Changes in version 1.1.4 (2018-5-1)}{
  \itemize{
    \item improved the `sc_detect_bc` function. add white list which is a list of known 10x barcode
    also added a new argument `number_of_cells` to limit the maximum number of cell barcodes in the result.
  }
}