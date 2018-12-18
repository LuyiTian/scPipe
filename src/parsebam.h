#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <unordered_map>
#include <thread>
#include <algorithm>
#include <Rcpp.h>
#include "config_hts.h"
#include "utils.h"
#include "cellbarcode.h"
#include "htslib/thread_pool.h"

#ifndef PARSEBAM_H
#define PARSEBAM_H

// overall_count_stat: count statistics for all cells
// chr_aligned_stat: per chromosome count
// per cell statistics:
// cell_mapped: mapped to exon
// cell_mapped_intron: mapped to intron
// cell_mapped_ambiguous: ambiguous map to multiple exon
// cell_align_unmapped: aligned to genome but not mapped to any feature
// cell_unaligned: not aligned
//
class Bamdemultiplex
{
public:
    Barcode bar;
    std::string c_tag;
    std::string m_tag;
    std::string g_tag;
    std::string a_tag;
    std::string out_dir;
    std::string mt_tag;

    std::unordered_map<std::string, int> overall_count_stat;
    std::unordered_map<std::string, int> chr_aligned_stat;
    std::unordered_map<std::string, int> cell_mapped_exon;
    std::unordered_map<std::string, int> cell_mapped_intron;
    std::unordered_map<std::string, int> cell_mapped_ambiguous;
    std::unordered_map<std::string, int> cell_align_unmapped;
    std::unordered_map<std::string, int> cell_unaligned;
    std::unordered_map<std::string, int> cell_ERCC;
    std::unordered_map<std::string, int> cell_MT;

    Bamdemultiplex(
        std::string odir,
        Barcode b,
        std::string cellular_tag,
        std::string molecular_tag,
        std::string gene_tag,
        std::string map_tag,
        std::string MT_tag
    );
    int barcode_demultiplex(std::string bam_path, int max_mismatch, bool has_UMI, int nthreads);
    int clean_bam_barcode(std::string bam_path, std::string out_bam, int max_mismatch, int nthreads);
    void write_statistics(
        std::string overall_stat_f,
        std::string chr_stat_f,
        std::string cell_stat_f
    );
};

#endif
