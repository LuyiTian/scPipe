// general functions
#include <string>
#include <fstream>
#include <stdexcept>
#include <map>
#include <sstream>
#include <vector>
#include <iomanip>
#include <zlib.h>
#include <Rcpp.h>
#include "config_hts.h"
#include "Gene.h"
#include "global_config.h"


#ifndef UTILS_H
#define UTILS_H

typedef std::pair<std::string, int> umi_pos_pair;

// join two paths, add separator if p1 does not end with it.
std::string join_path(const std::string p1, const std::string p2);

// calculate hamming distance of two strings A and B
// two strings should have equal length
int hamming_distance(const std::string &A, const std::string &B);

// calculate edit distance of two strings A and B (allow indel)
int edit_distance(const std::string& A, const std::string& B);

// since htslib does not validate file status we 
// need to check ourself
// if file not exist throw an exception
void check_file_exists (std::string name);

// count times of occurrence in a string vector
std::map<umi_pos_pair, int> vector_counter(const std::vector<umi_pos_pair> &v);

// split a string by given delimiter
// from: http://stackoverflow.com/questions/236129/split-a-string-in-c
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);

std::string padding(int count, int zero_num);

// stops program when file cannot be opened
void file_error(const char *filename);
// check if a file exists
bool fileExists(const std::string& filename);

char* getFileName(const char* path, const char* seperator = "/");
char* createFileWithAppend(const char *fq_out, const char* appendR1, const char *fq1_fn);
void openFile(gzFile &o_stream_gz_R1, std::ofstream &o_stream_R1, char* fqoutR1, bool write_gz);
#endif

