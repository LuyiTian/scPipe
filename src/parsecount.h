// parsecount.h
#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <map>
#include <algorithm>
#include <Rcpp.h>
#include "utils.h"
#include "cellbarcode.h"

#ifndef PARSECOUNT_H
#define PARSECOUNT_H

const int MAX_UMI_DUP = 1000;

struct UMI_dedup_stat
{
    int filtered_gene;
    int corrected_UMI;
    double A_prop;
    double T_prop;
    double G_prop;
    double C_prop;
};

std::unordered_map<std::string, std::vector<umi_pos_pair>> read_count(std::string fn, char sep);

int UMI_correct1(std::unordered_map<umi_pos_pair, int>& UMI_count); // sequence
int UMI_correct2(std::unordered_map<umi_pos_pair, int>& UMI_count); // sequence + position
int UMI_correct3(std::unordered_map<umi_pos_pair, int>& UMI_count); // sequence (2 edit distance)

std::unordered_map<std::string, int> UMI_dedup(
    std::unordered_map<std::string, std::vector<umi_pos_pair>> gene_read,
    std::vector<int>& UMI_dup_count,
    UMI_dedup_stat& s,
    int UMI_correct,
    bool read_filter
);

void write_mat(std::string fn, std::unordered_map<std::string, std::vector<int>> gene_cnt_matrix, std::vector<std::string> cellid_list);

void write_stat(std::string cnt_fn, std::string stat_fn, std::vector<int> UMI_dup_count, std::unordered_map<std::string, UMI_dedup_stat> UMI_dedup_stat_dict);

void get_counting_matrix(Barcode bar, std::string in_dir, int UMI_correct, bool read_filter);
#endif
