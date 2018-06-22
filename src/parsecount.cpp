// parse count
#include "parsecount.h"

using namespace Rcpp;

using std::getline;
using std::ifstream;
using std::make_pair;
using std::map;
using std::ofstream;
using std::stoi;
using std::string;
using std::stringstream;
using std::unordered_map;
using std::vector;

unordered_map<string, vector<umi_pos_pair>> read_count(string fn, char sep)
{
    ifstream infile(fn);
    unordered_map<string, vector<umi_pos_pair>> gene_read;
    string line;
    getline(infile, line); // skip header
    
    while(getline(infile, line))
    {
        size_t comma1_pos = line.find(',');
        size_t comma2_pos = line.find(',', comma1_pos + 1);

        string gene_id = line.substr(0, comma1_pos);
        string UMI = line.substr(comma1_pos + 1, comma2_pos - comma1_pos - 1);
        int pos = stoi(line.substr(comma2_pos + 1));

        gene_read[gene_id].push_back(make_pair(UMI, pos));
    }
    infile.close();
    return gene_read;
}

int UMI_correct1(map<umi_pos_pair, int>& UMI_count)
{
    bool found = false;
    int corrected_UMI = 0;
    for (auto UMI1 = UMI_count.begin(); UMI1 != UMI_count.end();) // use normal iterator in order to use `erase`
    {
        found = false;
        for (auto const& UMI2: UMI_count) // use range based
        {
            if (hamming_distance(UMI1->first.first, UMI2.first.first) == 1) // sequencing errors
            {
                if (UMI1->second == 1 || UMI1->second < UMI2.second*0.1)
                {
                    found = true;
                    // merge two UMIs
                    UMI_count[UMI2.first] += UMI_count[UMI1->first];
                    if (__DEBUG){Rcout << "merge: " <<  "<"<< UMI1->first.first << ", " << UMI1->first.second << ">" << "::" << "<"<< UMI2.first.first << ", " << UMI2.first.second << ">" << "\t" << UMI1->second << "::" << UMI2.second << "\n";}
                    break;
                }
            } else if (UMI1->first.first == UMI2.first.first) // possiblely same sequence in different position
            {
                if ((UMI1->first.second != UMI2.first.second) && (UMI1->second == 1 || UMI1->second < UMI2.second*0.5))
                {
                    found = true;
                    // merge two UMIs
                    UMI_count[UMI2.first] += UMI_count[UMI1->first];
                    if (__DEBUG){Rcout << "merge: " <<  "<"<< UMI1->first.first << ", " << UMI1->first.second << ">" << "::" << "<"<< UMI2.first.first << ", " << UMI2.first.second << ">" << "\t" << UMI1->second << "::" << UMI2.second << "\n";}
                    break;
                }
            }
        }
        if (found)
        {
            // delete UMI1
            corrected_UMI++;
            UMI1 = UMI_count.erase(UMI1);
        }
        else
        {
            UMI1++;
        }
    }
    return corrected_UMI;
}

int UMI_correct2(map<umi_pos_pair, int>& UMI_count)
{
    bool found = false;
    int corrected_UMI = 0;

    // iterate backwards to get better erase() performance
    // erase must move all elements in the tail, starting from the tail is much
    // more efficient than starting from the beginning
    for (auto UMI1 = UMI_count.end(); UMI1 != UMI_count.begin();) // use normal iterator in order to use `erase`
    {
        found = false;
        for (auto const& UMI2: UMI_count) // use range based
        {   
            if (abs(UMI1->first.second - UMI2.first.second) < 2)
            {
                if (UMI1->second == 1 || 2 * UMI1->second < UMI2.second)
                {
                    if (hamming_distance(UMI1->first.first, UMI2.first.first) == 1) // sequencing errors
                    {
                        found = true;
                        // merge two UMIs
                        UMI_count[UMI2.first] += UMI_count[UMI1->first];
                        if (__DEBUG){Rcout << "merge: " <<  UMI1->first.first << "::" << UMI2.first.first << "\t" << UMI1->second << "::" << UMI2.second << "\n";}
                        break;
                    }
                    else if ((UMI1->first.second != UMI2.first.second) && (UMI1->first.first == UMI2.first.first)) // diff pos
                    {
                        found = true;
                        // merge two UMIs
                        UMI_count[UMI2.first] += UMI_count[UMI1->first];
                        if (__DEBUG){Rcout << "merge: " <<  UMI1->first.first << "::" << UMI2.first.first << "\t" << UMI1->second << "::" << UMI2.second << "\n";}
                        break;
                    }
                }
            }
        }
        if (found)
        {
            // delete UMI1
            corrected_UMI++;
            UMI1 = UMI_count.erase(UMI1);
        }
        else
        {
            UMI1--;
        }
    }
    return corrected_UMI;
}


unordered_map<string, int> UMI_dedup(
    unordered_map<string, vector<umi_pos_pair>> gene_read,
    vector<int>& UMI_dup_count,
    struct UMI_dedup_stat& dedup_stat,
    int UMI_correct,
    bool read_filter
)
{
    unordered_map<string, int> gene_counter;

    for(auto const& a_gene: gene_read)
    {
        if (read_filter && a_gene.second.size() == 1)
        {
            dedup_stat.filtered_gene++;
            continue;
        }

        map<umi_pos_pair, int> UMI_count = vector_counter(a_gene.second);
        if (UMI_correct == 1)
        {
            dedup_stat.corrected_UMI += UMI_correct1(UMI_count);
        }
        else if (UMI_correct == 2)
        {
            dedup_stat.corrected_UMI += UMI_correct2(UMI_count);
        }
        else
        {
            Rcpp::stop("not implemented\n");
        }

        for (auto const& UMI: UMI_count)
        {
            if (UMI.second > MAX_UMI_DUP)
            {
                UMI_dup_count[MAX_UMI_DUP] ++;
            }
            else
            {
                UMI_dup_count[UMI.second-1] ++;
            }

        }

        //TODO: add ATCG percentage
        gene_counter[a_gene.first] = UMI_count.size();
    }

    return gene_counter;
}

void write_mat(string fn, unordered_map<string, vector<int>> gene_cnt_matrix, vector<string> cellid_list)
{
    ofstream o_file(fn);
    //write header
    o_file << "gene_id";
    for (auto const& ce : cellid_list)
    {
        o_file << "," << ce;
    }
    o_file << "\n";

    for (auto const& ge : gene_cnt_matrix)
    {
        o_file << ge.first;
        for (auto const& n : ge.second)
        {
            o_file << "," << n;
        }
        o_file << "\n";
    }
    o_file.close();
}

void write_stat(string cnt_fn, string stat_fn, vector<int> UMI_dup_count, unordered_map<string, UMI_dedup_stat> UMI_dedup_stat_dict)
{
    ofstream cnt_file(cnt_fn);
    cnt_file << "duplication number,count" << "\n";
    for (int i=0; i<UMI_dup_count.size(); i++)
    {
        cnt_file << i+1 << "," << UMI_dup_count[i] << "\n";
    }
    cnt_file.close();

    ofstream stat_file(stat_fn);
    stat_file << "cell_id,number of filtered gene,number of corrected UMI,UMI A percentage,UMI T percentage,UMI G percentage,UMI C percentage" << "\n";
    for (auto const& n : UMI_dedup_stat_dict)
    {
        stat_file << n.first << "," << n.second.filtered_gene << "," << n.second.corrected_UMI << "," \
            << n.second.A_prop << "," << n.second.T_prop << "," << n.second.G_prop << "," << n.second.C_prop << "," << "\n";
    }
    stat_file.close();
}

void get_counting_matrix(Barcode bar, string in_dir, int UMI_correct, bool read_filter)
{
    char sep = ',';
    unordered_map<string, string> cnt_files = bar.get_count_file_path(join_path(in_dir, "count"));
    unordered_map<string, vector<int>> gene_cnt_matrix; // store gene count matrix
    vector<string> all_gene_list; // store all gene ids
    vector<int> UMI_dup_count(MAX_UMI_DUP+1, 0); // store UMI duplication statistics
    unordered_map<string, UMI_dedup_stat> UMI_dedup_stat_dict;
    int cell_number = bar.cellid_list.size();
    int ind = 0;
    for (auto const& ce : bar.cellid_list) // for each cell
    {
        UMI_dedup_stat_dict[ce] = {}; // init zero
        unordered_map<string, vector<umi_pos_pair>> gene_read = read_count(cnt_files[ce], sep);
        unordered_map<string, int> gene_cnt =  UMI_dedup(gene_read, UMI_dup_count, UMI_dedup_stat_dict[ce], UMI_correct, read_filter);

        for (auto const& ge : gene_cnt) // for each gene
        {
            if (gene_cnt_matrix.find(ge.first) == gene_cnt_matrix.end())
            {
                auto & vec = gene_cnt_matrix[ge.first];
                vec.resize(cell_number, 0); // init with all zeros
                vec[ind] = ge.second;
            }
            else
            {
                gene_cnt_matrix[ge.first][ind] = ge.second;
            }
        }
        ind++;
    }

    // write to file
    write_mat(join_path(in_dir, "gene_count.csv"), gene_cnt_matrix, bar.cellid_list);
    string stat_dir = join_path(in_dir, "stat");
    write_stat(join_path(stat_dir, "UMI_duplication_count.csv"), join_path(stat_dir, "UMI_dedup_stat.csv"), UMI_dup_count, UMI_dedup_stat_dict);

}

