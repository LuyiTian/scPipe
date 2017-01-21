// demultiplexing
#include "cellbarcode.h"


// read annotation from files, the file should have two columns
// first column is cell id and second column is barcode.
// protocols with two barcodes are merged as one.
void Barcode::read_anno(std::string fn)
{
    std::ifstream infile(fn);
    std::string line;
    std::getline(infile, line); // skip header
    char sep;
    if (std::string::npos != line.find(","))
    {
        sep = ',';
    }
    else if (line.find("\t") != std::string::npos)
    {
        sep = '\t';
    }
    else
    {
        std::cout << "the annotation file should be comma or tab separated" << std::endl;
        exit(-1);
    }
    while(std::getline(infile, line))
    {
        std::stringstream linestream(line);
        std::string cell_id;
        std::string barcode;
        std::getline(linestream, cell_id, sep);
        std::getline(linestream, barcode, sep);

        barcode_dict[barcode] = cell_id;
        if (std::find(cellid_list.begin(), cellid_list.end(), cell_id) == cellid_list.end())
        {
            cellid_list.push_back(cell_id);
        }
        barcode_list.push_back(barcode);
    }
}


std::unordered_map<std::string, std::ofstream> Barcode::get_count_file_w(std::string out_dir)
{
    std::string csv_fmt = ".csv";
    std::unordered_map<std::string, std::ofstream> outfn_dict;
    for(const auto& n : cellid_list) 
    {
        outfn_dict[n].open(join_path(out_dir, n+csv_fmt));
        outfn_dict[n] << "gene_id,UMI,position" << std::endl;
    }
    return outfn_dict;
}

std::unordered_map<std::string, std::ifstream> Barcode::get_count_file_r(std::string in_dir)
{
    std::string csv_fmt = ".csv";
    std::unordered_map<std::string, std::ifstream> infn_dict;
    for(const auto& n : cellid_list) 
    {
        infn_dict[n].open(join_path(in_dir, n+csv_fmt));
    }
    return infn_dict;
}


std::string Barcode::get_closest_match(std::string bc_seq, int max_mismatch)
{
    if (barcode_dict.find(bc_seq) != barcode_dict.end())
    {
        return bc_seq;
    }
    int sml1st = 99999, sml2ed = 99999;
    int tmp;
    std::string cloest_match;



    for ( auto &barcode : barcode_list )
    {
        tmp = hamming_distance(barcode, bc_seq);
        if (tmp <= max_mismatch)
        {
            if (tmp < sml1st)
            {
                sml1st = tmp;
                cloest_match = barcode;
            }
            else if (sml1st < tmp <= sml2ed)
            {
                sml2ed = tmp;
            }
        }
    }
    if (sml1st < sml2ed)
    {
        return cloest_match;
    }
    else
    {
        return std::string();
    }

}



std::ostream& operator<< (std::ostream& out, const Barcode& obj)
{
    for( const auto& n : obj.barcode_dict ) 
    {
        out << "Barcode:[" << n.first << "] Cell Id:[" << n.second << "]\n";
    }
    return out;
}