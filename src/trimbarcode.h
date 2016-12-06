#include <zlib.h> // for reading compressed .fq file
#include <string>
#include <stdio.h>
#include <iostream>
#include <htslib/sam.h>
#include <htslib/kseq.h>
#include "utils.h"

#ifndef TRIMBARCODE_H
#define TRIMBARCODE_H
KSEQ_INIT(gzFile, gzread)
#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2) )

static const char empty_header[] = "@HD\tVN:1.4\tSO:unknown\n";

// the read structure, gives information where we can find the barcode and UMI.
// if there is no UMI, set umi_st = -1
struct read_s
{
    int id1_st;
    int id1_len;
    int id2_st;
    int id2_len;
    int umi_st;
    int umi_len;
};

struct filter_s
{
    bool if_check_qual;
    bool if_remove_N;
    int min_qual;
    int num_below_min;
};


void kseq_t_to_bam_t(kseq_t *seq, bam1_t *b, int trim_n);
//void paired_fastq_to_bam(char *fq1_fn, char *fq2_fn, char *bam_out, read_s r, filter_s fl, std::string cellular_tag, std::string molecular_tag);
void paired_fastq_to_bam(char *fq1_fn, char *fq2_fn, char *bam_out, read_s r, filter_s fl);
#endif