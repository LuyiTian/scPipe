#ifndef TRIMBARCODE_H
#define TRIMBARCODE_H

#include <zlib.h> // for reading compressed .fq file
#include <string>
#include <stdio.h>
#include <iostream>
#include <Rcpp.h>
#include "config_hts.h"
#include "utils.h"

// #ifndef INIT_KSEQ
// #define INIT_KSEQ
// KSEQ_INIT(gzFile, gzread)

// // static void REMOVE_KSEQ_WARNINGS(void) {
// //     (void)&kseq_read;
// //     // (void)&KSEQ_INIT;
// //     (void)&kseq_init;
// //     (void)&kseq_destroy;
// //     return;
// // }
// #endif

#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2) )

static const char empty_header[] = "@HD\tVN:1.4\tSO:unknown\n";

enum MatchType { Exact, Partial, NoMatch };

// Read structure
struct read_s
{
    int id1_st;     // id1 start
    int id1_len;    // id1 length
    int id2_st;     // id2 start
    int id2_len;    // id2 length
    int umi_st;     // umi start
    int umi_len;    // umi length
};

// Filter settings
struct filter_s
{
    bool if_check_qual;
    bool if_remove_N;
    int min_qual;
    int num_below_min;
};

// Conversion functions
//void kseq_t_to_bam_t(kseq_t *seq, bam1_t *b, int trim_n);
void paired_fastq_to_bam(char *fq1_fn, char *fq2_fn, char *bam_out, const read_s read_structure, const filter_s filter_settings);
void paired_fastq_to_fastq(char *fq1_fn, char *fq2_fn, char *fq_out, const read_s read_structure, const filter_s filter_settings, const bool write_gz);
void single_fastq_to_fastq(char *fq1_fn, char *fq_out, const read_s read_structure, const filter_s filter_settings);

std::vector<int> sc_atac_paired_fastq_to_fastq(
        const char *fq1_fn,
        std::vector<std::string> fq2_fn_list,
        const char *fq3_fn,
		const char *valid_barcode_fn,
        const char *fq_out,
        const bool write_gz,
        const bool rmN,
        const bool rmlow,
        int min_qual,
        int num_below_min,
		bool no_reverse_complement);

std::vector<int> sc_atac_paired_fastq_to_csv(
        char *fq1_fn,
        char *fq3_fn,
        char *fq_out, 
        char *bc_fn, 
		char *valid_barcode_fn,
        int umi_start,
        int umi_length,
        char *umi_in,
        const bool write_gz,
        const bool rmN,
        const bool rmlow,
        int min_qual,
        int num_below_min,
        int id1_st,
        int id1_len,
        int id2_st,
        int id2_len);

#endif
