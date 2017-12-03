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

\section{Changes in version 1.0.2 (2017-12-01)}{
  \itemize{
    \item fix a bug in `validObject`. the default value for gene id and organism is set to NA
  }
}