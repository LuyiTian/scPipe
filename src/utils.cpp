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
    } else {
        f.close();
        throw std::invalid_argument("cannot open file: "+ fn + "\n");
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

// Constructors for Interval
Interval::Interval(int s, int e): st(s), en(e), snd(0) {}
Interval::Interval(int s, int e, int sd): st(s), en(e), snd(sd) {}

int Interval::overlap(int st1, int en1)
{
  if (en1 < st)
  {
    return -1;
  }
  else if (st1 > en)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

bool operator < (const Interval &L, const Interval &R)
{
  return L.en < R.st;
}

bool operator > (const Interval &L, const Interval &R)
{
  return L.st > R.en;
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