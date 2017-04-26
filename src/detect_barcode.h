#include <zlib.h> // for reading compressed .fq file
#include <string>
#include <stdio.h>
#include <iostream>
#include <limits>
#include "config_hts.h"
#include "utils.h"

#ifndef DETECTBARCODE_H
#define DETECTBARCODE_H
KSEQ_INIT(gzFile, gzread)

std::unordered_map<std::string, int> summarize_barcode(std::string fn, int bc_len, int max_reads, int max_mismatch, int min_count);
void write_barcode_summary(std::string outfn, std::string surfix, std::unordered_map<std::string, int> counter);
#endif