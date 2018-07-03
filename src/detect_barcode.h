#include <zlib.h> // for reading compressed .fq file
#include <string>
#include <stdio.h>
#include <iostream>
#include <unordered_map>
#include <limits>
#include <Rcpp.h>
#include "config_hts.h"
#include "utils.h"

#ifndef INIT_KSEQ
#define INIT_KSEQ
KSEQ_INIT(gzFile, gzread)
#endif

#ifndef DETECTBARCODE_H
#define DETECTBARCODE_H

std::unordered_map<std::string, int> summarize_barcode(
    std::string filename,
    int bc_len,
    int max_reads,
    int max_mismatch,
    int min_count,
    std::string whitelist_fn
);

void write_barcode_summary(
    std::string outfn,
    std::string prefix,
    std::unordered_map<std::string, int> counter,
    int number_of_cells
);
#endif
