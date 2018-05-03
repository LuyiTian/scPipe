#include "utils.h"

using namespace Rcpp;

using std::ifstream;
using std::invalid_argument;
using std::ostringstream;
using std::setfill;
using std::setw;
using std::string;
using std::unordered_map;
using std::vector;

string join_path(const string p1, const string p2)
{
    auto sep = '/';
	return p1.back() == sep ? p1 + p2 : p1 + sep + p2;
}

int hamming_distance(string A, string B)
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

void check_file_exists(string fn)
{
    ifstream f(fn.c_str());
    if (f.good()) {
        f.close();
    }
    else
    {
        f.close();
        throw invalid_argument("cannot open file: " + fn + "\n");
    }   
}

// tally the element in vector
unordered_map<umi_pos_pair, int> vector_counter(vector<umi_pos_pair> v)
{
    unordered_map<umi_pos_pair, int> counter;
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

vector<string> &split(const string &s, char delim, vector<string> &elems)
{
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim)
{
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

string padding(int count, int zero_num)
{
    ostringstream ss;
    ss << setw(zero_num) << setfill('0') << count;
    return ss.str();
}

void file_error(char *filename) {
    stringstream err_msg;
    err_msg << "Can't open file: " << filename << "\n";
    Rcpp::stop(err_msg.str());
}
