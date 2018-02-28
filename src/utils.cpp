#include "utils.h"

using namespace Rcpp;

std::string join_path(const std::string p1, const std::string p2)
{
    auto sep = '/';
	return p1.back() == sep ? p1 + p2 : p1 + sep + p2;
}

int hamming_distance(std::string A, std::string B)
{
    int dist = 0;
    for (int i = 0; i < A.length(); ++i)
    {
        if (A[i] != B[i])
        {
            dist++;
        }
    }
    return dist;
}

void check_file_exists(std::string fn)
{
    std::ifstream f(fn.c_str());
    if (f.good()) {
        f.close();
    }
    else
    {
        f.close();
        throw std::invalid_argument("cannot open file: " + fn + "\n");
    }   
}

// tally the element in vector
std::unordered_map<std::string, int> vector_counter(std::vector<std::string> v)
{
    std::unordered_map<std::string, int> counter;
    for(auto const& val: v)
    {
        if (counter.find(val) != counter.end())
        {
            counter[val] ++;
        }
        else
        {
            counter[val] = 1;
        }
    }

    return counter;
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::string padding(int count, int zero_num)
{
  std::ostringstream ss;
  ss << std::setw(zero_num) << std::setfill('0') << count;
  return ss.str();
}

void file_error(char *filename) {
    std::stringstream err_msg;
    err_msg << "Can't open file: " << filename << "\n";
    Rcpp::stop(err_msg.str());
}
