#include "utils.h"


std::string join_path(const std::string p1, const std::string p2) {

   char sep = '/';
   std::string tmp = p1;


  if (p1[p1.length()-1] != sep) { // need to add a path separator
     tmp += sep;                
     return(tmp + p2);
  }
  else
     return(p1 + p2);
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


void check_file_exists (std::string fn) {
    std::ifstream f(fn.c_str());
    if (f.good()) {
        f.close();
    } else {
        f.close();
        throw std::invalid_argument("cannot open file: "+ fn + "\n");
    }   
}


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

Interval::Interval(int s, int e): st(s), en(e), snd(0){}
Interval::Interval(int s, int e, int sd): st(s), en(e), snd(sd){}


int Interval::overlap(int st1, int en1)
{
  int s, e;
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

//bool Interval::operator<(Interval & R)
//{
//  return en < R.st;
//}
//
//
//bool Interval::operator>(Interval & R)
//{
//  return st > R.en;
//}




std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}