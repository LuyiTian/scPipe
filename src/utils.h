// general functions
#include <string>
#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <sstream>
#include <vector>
#include <iomanip>
#include <Rcpp.h>
#include "Gene.h"
#include "global_config.h"

#ifndef UTILS_H
#define UTILS_H

// join two paths, add separator if p1 does not end with it.
std::string join_path(const std::string p1, const std::string p2);

// calculate hamming distance of two strings A and B
// two strings should have equal length
int hamming_distance(std::string A, std::string B);

// since htslib does not validate file status we 
// need to check ourself
// if file not exist throw an exception
void check_file_exists (std::string name);

// count times of occurrence in a string vector
std::unordered_map<std::string, int> vector_counter(std::vector<std::string> v);

// split a string by given delimiter
// from: http://stackoverflow.com/questions/236129/split-a-string-in-c
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);

std::string padding(int count, int zero_num);

#endif