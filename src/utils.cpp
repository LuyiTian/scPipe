#include "utils.h"

#include "config_hts.h"
#include "Gene.h"
#include "global_config.h"

#include <string>
#include <fstream>
#include <stdexcept>
#include <map>
#include <sstream>
#include <vector>
#include <iomanip>
#include <zlib.h>
#include <Rcpp.h>
#include <sys/stat.h>

using namespace Rcpp;

using std::ifstream;
using std::invalid_argument;
using std::ostringstream;
using std::stringstream;
using std::setfill;
using std::setw;
using std::string;
using std::map;
using std::vector;

string join_path(const string p1, const string p2)
{
    auto sep = '/';
	return p1.back() == sep ? p1 + p2 : p1 + sep + p2;
}

int hamming_distance(const string &A, const string &B)
{
    int dist = 0;
    for (unsigned int i = 0; i < A.length(); ++i)
    {
        dist += (A[i] != B[i]);
    }
    return dist;
}

// not a fast implementation.
int edit_distance(const string& A, const string& B)
{
    int NA = A.size();
    int NB = B.size();
    double x, y, z;
    
    vector<vector<int>> M(NA + 1, vector<int>(NB + 1));
    
    for (int i = 0; i <= NA; ++i)
        M[i][0] = i;
    
    for (int i = 0; i <= NB; ++i)
        M[0][i] = i;
    
    for (int a = 1; a <= NA; ++a)
    {
        for (int b = 1; b <= NB; ++b)
        {
            x = M[a - 1][b] + 1;
            y = M[a][b - 1] + 1;
            z = M[a - 1][b - 1] + (A[a - 1] == B[b - 1] ? 0 : 1);
            
            M[a][b] = std::min(std::min(x, y), z);
        }
    }
    
    
    return M[A.size()][B.size()];
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
map<umi_pos_pair, int> vector_counter(const vector<umi_pos_pair> &v)
{
    map<umi_pos_pair, int> counter;
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

void file_error(const char *filename) {
    stringstream err_msg;
    err_msg << "Can't open file: " << filename << "\n";
    Rcpp::stop(err_msg.str());
}

// check file exists
bool fileExists(const std::string& filename) {
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1) {
        return true;
    }
    return false;
}



char* getFileName(const char* cpath, const char* seperator)
{
	char *path = (char *)malloc(std::strlen(cpath) * sizeof(char *));
	std::strcpy(path, cpath);
    char *ssc;
    int l = 0;
    ssc = std::strstr(path, seperator);
    while (ssc) {
        l = strlen(ssc) + 1;
        path = &path[strlen(path)-l+2];
        ssc = strstr(path, seperator);
    }
    return path;
}

char* createFileWithAppend(const char *fq_out, const char* appendR1, const char *fq1_fn){
    char* fqoutR1 = (char*)malloc(strlen(fq_out) + strlen(appendR1) + strlen(getFileName(fq1_fn)) + 1);
    strcpy(fqoutR1, fq_out);
    strcat(fqoutR1,appendR1);
    strcat(fqoutR1,getFileName(fq1_fn));
    return fqoutR1;
}

void openFile(gzFile &o_stream_gz_R1,std::ofstream &o_stream_R1,char* fqoutR1, bool write_gz){
    
    if (write_gz) {
        o_stream_gz_R1 = gzopen(fqoutR1, "wb2"); // open gz file
        if (!o_stream_gz_R1) {
            file_error(fqoutR1);
        }
        
    } else {
        o_stream_R1.open(fqoutR1); // output file
        if (!o_stream_R1.is_open()) {
            file_error(fqoutR1);
        }
    }
}
