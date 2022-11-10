//detect_barcode.cpp

#include "detect_barcode.h"

#ifndef INIT_KSEQ
#define INIT_KSEQ
KSEQ_INIT(gzFile, gzread)
#endif

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
                    if (__DEBUG) {Rcpp::Rcout << "merge: " <<  bc1->first << "::" << bc2.first << "\t" << bc1->second << "::" << bc2.second << "\n";}
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


void merge_barcode_using_whitelist(std::unordered_map<std::string, int> &counter, int max_mismatch, int min_count, std::string whitelist_fn)
{
    std::string line;
    std::unordered_map<std::string, int> whitelist;
    std::ifstream wl_file(whitelist_fn);
    while (std::getline(wl_file, line))
    {
       whitelist[line] = 0;
    }

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

    for (auto const& bc: counter)  // compare against whitelist
    {
        if(whitelist.find(bc.first) != whitelist.end())
        {
            whitelist[bc.first] += counter[bc.first];
        }
        else
        {
            for (auto const& wt: whitelist)
            {
                if(hamming_distance(wt.first,bc.first)<=max_mismatch)
                {
                    whitelist[wt.first] += counter[bc.first];
                    break;
                }
            }
        }
    }

    for (auto wt = whitelist.begin(); wt != whitelist.end();) // use normal iterator in order to use `erase`
    {
        if (wt->second == 0)
        {
            wt = whitelist.erase(wt);
        }
        else
        {
            wt++;
        }
    }
    counter = whitelist;
    wl_file.close();
}



std::unordered_map<std::string, int> summarize_barcode(std::string fn, int bc_len, int max_reads, int max_mismatch, int min_count, std::string whitelist_fn)
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
        // std::cerr << "max_reads should be larger than 1000." << "\n";
        Rcpp::stop("max_reads should be larger than 1000.");
    }
    if (bc_len < 4)
    {
        // std::cerr << "bc_len should be larger than 3." << "\n";
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

    if(whitelist_fn.length()>1)
    {
        merge_barcode_using_whitelist(counter, max_mismatch, min_count, whitelist_fn);
    }
    else
    {
        merge_barcode(counter, max_mismatch, min_count);
    }

    return counter;
}

void write_barcode_summary(std::string outfn, std::string prefix, std::unordered_map<std::string, int> counter, int number_of_cells)
{
    // open output file and write header
    std::ofstream o_file(outfn);
    o_file << "cell_name" << "," << "barcode_sequence" << "," << "count" << "\n";

    // if no cells requested, write file with empty rows
    if (number_of_cells <= 0) {
        return;
    }

    std::vector<std::pair<int, std::string>> items;
    for (auto const& bc: counter)
    {
        auto const &barcode = bc.first;
        auto const &bc_count = bc.second;
        items.push_back(std::make_pair(bc_count, barcode));
    }

    // sort barcodes by descending count, breaking ties lexically
    // the standard comparator for std::pair would be enough, as it compares first then second member
    std::sort(items.begin(), items.end()); 
    std::reverse(items.begin(), items.end()); // highest -> lowest

    // digits for padding cell number with 0s
    const int digits = std::to_string(counter.size()).length() + 1;
    int cell_ind = 1;
    for (auto const& bc: items)
    {
        o_file << prefix << padding(cell_ind, digits) << "," << bc.second << "," << bc.first << "\n";
        cell_ind++;
        if (cell_ind > number_of_cells)
        {
            break;
        }
    }
    o_file.close();
}
