// demultiplexing
#include "cellbarcode.h"

using std::string;
using std::unordered_map;
using namespace Rcpp;

// read annotation from files, the file should have two columns
// first column is cell id and second column is barcode.
// protocols with two barcodes are merged as one.
void Barcode::read_anno(string fn)
{
    check_file_exists(fn);
    std::ifstream infile(fn);
    string line;
    std::getline(infile, line); // skip header

    char sep;
    if (line.find(",") != string::npos)
    {
        sep = ',';
    }
    else if (line.find("\t") != string::npos)
    {
        sep = '\t';
    }
    else
    {
        // std::cout << "the annotation file should be comma or tab separated" << std::endl;
        Rcpp::stop("the annotation file should be comma or tab separated");
    }

    while(std::getline(infile, line))
    {
        std::stringstream linestream(line);
        string cell_id;
        string barcode;

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

unordered_map<string, string> Barcode::get_count_file_path(string out_dir)
{
    string csv_fmt = ".csv";
    unordered_map<string, string> out_fn_dict;
    for(const auto& n : cellid_list) 
    {
        out_fn_dict[n] = join_path(out_dir, n+csv_fmt);
    }
    return out_fn_dict;
}


string Barcode::get_closest_match(string bc_seq, int max_mismatch)
{
    if (barcode_dict.find(bc_seq) != barcode_dict.end())
    {
        return bc_seq;
    }
    int sml1st = 99999, sml2ed = 99999;
    int tmp;
    string closest_match;



    for ( auto &barcode : barcode_list )
    {
        tmp = hamming_distance(barcode, bc_seq);
        if (tmp <= max_mismatch)
        {
            if (tmp < sml1st)
            {
                sml1st = tmp;
                closest_match = barcode;
            }
            else if (sml1st < tmp <= sml2ed)
            {
                sml2ed = tmp;
            }
        }
    }
    if (sml1st < sml2ed)
    {
        return closest_match;
    }
    else
    {
        return string();
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