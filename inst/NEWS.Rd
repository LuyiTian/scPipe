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

\section{Changes in version 1.3.1 (y version bumped by bioconductor) (2018-5-3)}{
  \itemize{
    \item updated the `sc_exon_mappping` function so it can accept multiple bam files now, together with
    a list of cell id or cell barcode.
  }
}

\section{Changes in version 1.3.2 (2018-5-3)}{
  \itemize{
    \item Added gzipped output for `sc_trim_barcode()`, if output filename ends with `.gz` then gzipped output will be produced.
    \item Added `get_read_str()` function for getting common read structures.
  }
}

\section{Changes in version 1.3.4 (2018-5-14)}{
  \itemize{
    \item put `distance_to_end` back.
  }
}

\section{Changes in version 1.3.5 (2018-5-24)}{
  \itemize{
    \item bug fix
  }
}

\section{Changes in version 1.3.6 (2018-6-12)}{
  \itemize{
    \item support multiple bam file for the same sample pooling as input for exon mapping
  }
}

\section{Changes in version 1.3.7 (2018-6-23)}{
  \itemize{
    \item added multithreading to `parse_align`
    \item added nthreads argument to `parse_align` and all upstream wrappers
  }
}

\section{Changes in version 1.3.8 (2018-6-27)}{
  \itemize{
    \item fix a minor bug in cell barcode demultiplexing
    \item corrected `suffix` argument in `sc_detect_bc` to `prefix`
  }
}

\section{Changes in version 1.3.10 (2018-8-29)}{
  \itemize{
    \item added option to use GRanges objects as annotation
    \item added option to use SAF formatted data.frame objects as annotation
    \item added support for gzipped annotation files
  }
}

\section{Changes in version 1.3.11 (2018-9-20)}{
  \itemize{
    \item added `get_ercc_anno()` function to retrieve ERCC annotation as data.frame
    \item fixed `sc_detect_bc()` returning one more cell than requested
  }
}

\section{Changes in version 1.3.13 (2018-10-22)}{
  \itemize{
    \item fixed bug in `create_report()` due to `changes` in scater parameters
  }
}

\section{Changes in version 1.3.16 (2018-10-31)}{
  \itemize{
    \item add travis integration and codecov report.
    \item fixed broken examples
  }
}

\section{Changes in version 1.5.2 (2018-12-18)}{
  \itemize{
    \item add new function `sc_correct_bam_bc`, which will update the cell barcode tag 
    in bam file with corrected barcode, using existing cell barcode annotation.
  }
}

\section{Changes in version 1.9.1 (2019-12-23)}{
  \itemize{
    \item Changed instances of deprecated `plotQC()` to `plotHighestExprs()`
  }
}

\section{Changes in version 1.9.4 (2020-04-06)}{
  \itemize{
    \item Changed hashmap package calls to hash package as hashmap is no longer on CRAN.
  }
}