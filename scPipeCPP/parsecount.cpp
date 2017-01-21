// parse count
#include "parsecount.h"


std::unordered_map<std::string, std::vector<std::string>> read_count(std::ifstream& infile, char sep)
{
    std::unordered_map<std::string, std::vector<std::string>> gene_read;
    std::string line;
    std::getline(infile, line); // skip header
    while(std::getline(infile, line))
    {
        std::stringstream linestream(line);
        std::string gene_id;
        std::string UMI;
        std::getline(linestream, gene_id, sep);
        std::getline(linestream, UMI, sep);

        gene_read[gene_id].push_back(UMI);

    }
    return gene_read;
}

int UMI_correct1(std::unordered_map<std::string, int>& UMI_count)
{
    bool found = false;
    int corrected_UMI = 0;
    for (auto UMI1 = UMI_count.begin(); UMI1 != UMI_count.end();) // use normal iterator in order to use `erase`
    {
        found = false;
        for (auto const& UMI2: UMI_count) // use range based
        {
            if (hamming_distance(UMI1->first, UMI2.first) == 1)
            {
                if (UMI1->second == 1 || UMI1->second < UMI2.second*0.5)
                {
                    found = true;
                    // merge two UMIs
                    UMI_count[UMI2.first] += UMI_count[UMI1->first];
                    // std::cout << "merge: " <<  UMI1->first << "::" << UMI2.first << "\t" << UMI1->second << "::" << UMI2.second << std::endl;
                    break;
                }
            }
        }
        if (found)
        {
            // delete UMI1
            corrected_UMI ++;
            UMI1 = UMI_count.erase(UMI1);
        }
        else
        {
            UMI1 ++;
        }
    }
    return corrected_UMI;
}


std::unordered_map<std::string, int> UMI_dedup(std::unordered_map<std::string, std::vector<std::string>> gene_read, std::vector<int>& UMI_dup_count, UMI_dedup_stat& s, int UMI_correct, bool read_filter)
{
    std::unordered_map<std::string, int> gene_counter;

    for(auto const& a_gene: gene_read)
    {
        if (read_filter && a_gene.second.size() == 1)
        {
            s.filtered_gene ++;
            continue;
        }

        std::unordered_map<std::string, int> UMI_count = vector_counter(a_gene.second);
        if (UMI_correct == 1)
        {
            s.corrected_UMI += UMI_correct1(UMI_count);
        }
        else
        {
            std::cout << "not implemented" << std::endl; 
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


void write_mat(std::string fn, std::unordered_map<std::string, std::vector<int>> gene_cnt_matrix, std::vector<std::string> cellid_list)
{
    std::ofstream o_file(fn);
    //write header
    o_file << "gene_id";
    for (auto const& ce : cellid_list)
    {
        o_file << "," << ce;
    }
    o_file << std::endl;

    for (auto const& ge : gene_cnt_matrix)
    {
        o_file << ge.first;
        for (auto const& n : ge.second)
        {
            o_file << "," << n;
        }
        o_file << std::endl;
    }
    o_file.close();
}


void write_stat(std::string cnt_fn, std::string stat_fn, std::vector<int> UMI_dup_count, std::unordered_map<std::string, UMI_dedup_stat> UMI_dedup_stat_dict)
{
    std::ofstream cnt_file(cnt_fn);
    cnt_file << "duplication number,count" << std::endl;
    for (int i=0; i<UMI_dup_count.size(); i++)
    {
        cnt_file << i+1 << "," << UMI_dup_count[i] << std::endl;
    }
    cnt_file.close();

    std::ofstream stat_file(stat_fn);
    stat_file << "cell_id,number of filtered gene,number of corrected UMI,UMI A percentage,UMI T percentage,UMI G percentage,UMI C percentage" << std::endl;
    for (auto const& n : UMI_dedup_stat_dict)
    {
        stat_file << n.first << "," << n.second.filtered_gene << "," << n.second.corrected_UMI << "," \
            << n.second.A_prop << "," << n.second.T_prop << "," << n.second.G_prop << "," << n.second.C_prop << "," << std::endl;
    }
    stat_file.close();
}


void get_counting_matrix(Barcode bar, std::string in_dir, int UMI_correct, bool read_filter)
{
    char sep = ',';
    std::unordered_map<std::string, std::ifstream> cnt_files = bar.get_count_file_r(join_path(in_dir,"count"));
    std::unordered_map<std::string, std::vector<int>> gene_cnt_matrix; // store gene counting matrix
    std::vector<std::string> all_gene_list; //store all gene ids
    std::vector<int> UMI_dup_count(MAX_UMI_DUP+1, 0); // store UMI duplication statistics
    std::unordered_map<std::string, UMI_dedup_stat> UMI_dedup_stat_dict;
    int cell_number = bar.cellid_list.size();
    int ind = 0;
    for (auto const& ce : bar.cellid_list) // for each cells
    {
        UMI_dedup_stat_dict[ce] = {}; // init with all zero
        std::unordered_map<std::string, std::vector<std::string>> gene_read = read_count(cnt_files[ce], sep);
        std::unordered_map<std::string, int> gene_cnt =  UMI_dedup(gene_read, UMI_dup_count, UMI_dedup_stat_dict[ce], UMI_correct, read_filter);

        for (auto const& ge : gene_cnt) // for each genes
        {
            if (gene_cnt_matrix.find(ge.first) == gene_cnt_matrix.end())
            {
                auto & vec = gene_cnt_matrix[ge.first];
                vec.resize(cell_number, 0); // init with all zero
                vec[ind] = ge.second; 
            }
            else
            {
                gene_cnt_matrix[ge.first][ind] = ge.second;
            }
        }
        ind ++;
    }

    // write to file
    write_mat(join_path(in_dir, "gene_count.csv"), gene_cnt_matrix, bar.cellid_list);
    std::string stat_dir = join_path(in_dir, "stat");
    write_stat(join_path(stat_dir, "UMI_duplication_count.csv"), join_path(stat_dir, "UMI_dedup_stat.csv"), UMI_dup_count, UMI_dedup_stat_dict);

}

