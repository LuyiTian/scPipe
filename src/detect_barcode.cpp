//detect_barcode.cpp

#include "detect_barcode.h"

using namespace Rcpp;

void merge_barcode(std::unordered_map<std::string, int> &counter, int max_mismatch, int min_count)
{
    for (auto bc1 = counter.begin(); bc1 != counter.end();) // use normal iterator in order to use `erase`
    {
        if (bc1->second < min_count)
        {
            bc1 = counter.erase(bc1);
        }
        else
        {
            bc1++;
        }
    }

    bool found;
    for (auto bc1 = counter.begin(); bc1 != counter.end();) // use normal iterator in order to use `erase`
    {
        found = false;
        for (auto const& bc2: counter) // use range based
        {
            if (hamming_distance(bc1->first, bc2.first) == 1)
            {
                if (bc1->second == max_mismatch || bc1->second < bc2.second*0.5)
                {
                    found = true;
                    // merge two barcodes
                    counter[bc2.first] += counter[bc1->first];
                    if (__DEBUG) {Rcpp::Rcout << "merge: " <<  bc1->first << "::" << bc2.first << "\t" << bc1->second << "::" << bc2.second << std::endl;}
                    break;
                }
            }
        }
        if (found)
        {
            // delete bc1
            bc1 = counter.erase(bc1);
        }
        else
        {
            bc1++;
        }
    }
}



std::unordered_map<std::string, int> summarize_barcode(std::string fn, int bc_len, int max_reads, int max_mismatch, int min_count)
{
    check_file_exists(fn);
    gzFile fq = gzopen(fn.c_str(), "r"); // input fastq
    std::unordered_map<std::string, int> counter;
    std::string tmp_bc;
    if (max_reads <= 0)
    {
        max_reads = std::numeric_limits<int>::max();
    }
    else if ((max_reads > 0) && (max_reads < 1000))
    {
        // std::cerr << "max_reads should be larger than 1000." << std::endl;
        Rcpp::stop("max_reads should be larger than 1000.");
    }
    if (bc_len < 4)
    {
        // std::cerr << "bc_len should be larger than 3." << std::endl;
        Rcpp::stop("bc_len should be larger than 3.");
    }

    kseq_t *seq;
    seq =  kseq_init(fq);
    int cnt = 0;
    int l1 = 0;
    while((cnt < max_reads) && ((l1 = kseq_read(seq))>= 0))
    {
        tmp_bc = std::string(seq->name.s).substr(0, bc_len);
        if (counter.find(tmp_bc) != counter.end())
        {
            counter[tmp_bc] ++;
        }
        else
        {
            counter[tmp_bc] = 1;
        }
        cnt++;
    }
    kseq_destroy(seq);
    gzclose(fq);

    merge_barcode(counter, max_mismatch, min_count);



    return counter;
}

void write_barcode_summary(std::string outfn, std::string surfix, std::unordered_map<std::string, int> counter)
{
    std::ofstream o_file(outfn);  // output file
    int cnt = 0;
    int dig = std::to_string(counter.size()).length()+1;
    for (auto const& bc: counter)
    {
        o_file << surfix << padding(cnt, dig) << "," << bc.first << "," << bc.second << std::endl;
        cnt++;
    }
}
    