// transcriptmapping.cpp
#include "transcriptmapping.h"

using std::string;
using namespace Rcpp;

int GeneAnnotation::get_strand(char st)
{
    int strand = 0;
    if (st == '+')
    {
        strand = 1;
    }
    else if (st == '-')
    {
        strand = -1;
    }
    return strand;
}

string GeneAnnotation::get_parent(string tok)
{
    string parent = "";
    auto subtoken = split(tok, ';');
    for (auto attr : subtoken)
    {
        if (attr.substr(0,6) == "Parent")
        {
            parent = split(split(attr, '=')[1],':')[1];
        }
    }
    return parent;
}

string GeneAnnotation::get_ID(string tok)
{
    string ID = "";
    auto subtoken = split(tok, ';');
    for (auto attr : subtoken)
    {
        if (attr.substr(0,2) == "ID")
        {
            ID = split(split(attr, '=')[1],':')[1];
        }
    }
    return ID;
}

string GeneAnnotation::fix_name(string na)
{
    string new_na;
    if (na.compare(0,3,"chr") == 0)
    {
        return na;
    }
    else if (na.length() > 4) // just fix 1-22, X, Y, MT. ignore contig and ERCC
    {
        return na;
    }
    else
    {
        if (na == "MT")
        {
            new_na = "chrM";
        }
        else
        {
            new_na = "chr"+na;
        }
        return new_na;
    }
}

void GeneAnnotation::parse_gff3_annotation(string gff3_fn, bool fix_chrname)
{
    std::ifstream infile(gff3_fn);

    string line;
    string ID;
    string parent;
    std::unordered_map<string, string> tmp_trans_dict; // store transcript - gene mapping
    std::unordered_map<string, std::unordered_map<string, Gene>> tmp_gene_dict;
    int strand = 0;
    std::vector<string> token;

    // Rcpp::Rcout << "build annotation from gff3 file..." << "\n";
    while(std::getline(infile, line))
    {
        if (line[0] == '#')
        {
            continue; // skip header
        }

        token = split(line, '\t');
        string chr_name = token[0];
        parent = get_parent(token[8]);
        strand = get_strand(token[6][0]);

        string type = token[2];
        if (type == "exon")
        {
            if (tmp_trans_dict.end() != tmp_trans_dict.find(parent))
            {
                if (fix_chrname)
                {
                    chr_name = fix_name(chr_name);
                }
                auto &genes_list = tmp_gene_dict[chr_name];

                int interval_start = std::atoi(token[3].c_str());
                int interval_end = std::atoi(token[4].c_str());

                string target = tmp_trans_dict[parent];
                genes_list[target].add_exon(Interval(interval_start, interval_end, strand));
                genes_list[target].set_ID(target);
            }
            else
            {
                std::stringstream err_msg;
                err_msg << "cannot find grandparent for exon:" << "\n";
                err_msg << line << "\n";
                Rcpp::stop(err_msg.str());
            }

        }
        else if (!parent.empty())
        {
            ID = get_ID(token[8]);
            if (!ID.empty())
            {
                tmp_trans_dict[ID] = parent;
            }
        }

    }

    for (auto iter : tmp_gene_dict)
    {
        for (auto sub_iter : iter.second)
        {
            sub_iter.second.sort_exon();
            gene_dict[iter.first].push_back(sub_iter.second);
        }
        std::sort(gene_dict[iter.first].begin(), gene_dict[iter.first].end());
    }

}

void GeneAnnotation::parse_bed_annotation(string bed_fn, bool fix_chrname)
{
    std::ifstream infile(bed_fn);

    string line;
    std::unordered_map<string, std::unordered_map<string, Gene>> tmp_gene_dict;
    int strand = 0;
    std::vector<string> token;

    std::getline(infile, line); // skip the header
    while(std::getline(infile, line))
    {
        token = split(line, '\t');
        strand = get_strand(token[4][0]);
        if (fix_chrname)
        {
            tmp_gene_dict[fix_name(token[1])][token[0]].add_exon(Interval(std::atoi(token[2].c_str()), std::atoi(token[3].c_str()), strand));
            tmp_gene_dict[fix_name(token[1])][token[0]].set_ID(token[1]);
        }
        else
        {
            tmp_gene_dict[token[1]][token[0]].add_exon(Interval(std::atoi(token[2].c_str()), std::atoi(token[3].c_str()), strand));
            tmp_gene_dict[token[1]][token[0]].set_ID(token[1]);
        }

    }

    for (auto iter : tmp_gene_dict)
    {
        for (auto sub_iter : iter.second)
        {
            if (sub_iter.second.exon_vec.size()>1)
            {
                sub_iter.second.sort_exon();
            }
            gene_dict[iter.first].push_back(sub_iter.second);
        }
        if (gene_dict[iter.first].size()>1)
        {
            std::sort(gene_dict[iter.first].begin(), gene_dict[iter.first].end());
        }
    }
}


int GeneAnnotation::ngenes()
{
    int gene_number = 0;
    for (auto iter : gene_dict)
    {
        for (auto sub_iter : iter.second)
        {
            gene_number ++;
        }
    }

    return gene_number;
}


std::vector<string> GeneAnnotation::get_genelist()
{
    std::vector<string> gene_list;
    for (auto iter : gene_dict)
    {
        for (auto sub_iter : iter.second)
        {
            gene_list.push_back(sub_iter.gene_id);
        }
    }

    return gene_list;
}


std::ostream& operator<< (std::ostream& out, const GeneAnnotation& obj)
{
    out << "annotation statistics:" << "\n";
    for ( const auto& n : obj.gene_dict )
    {
        out << "\t" << "chromosome:[" << n.first << "] number of genes:[" << n.second.size() << "]\n";
    }
    for ( const auto& n : obj.gene_dict )
    {
        out << "first gene in chromosome " << n.first << " :" << "\n";
        out << n.second[0] << "\n";
        //break;
    }
    return out;
}


void Mapping::add_annotation(string gff3_fn, bool fix_chrname)
{
    if (gff3_fn.substr(gff3_fn.find_last_of(".") + 1) == "gff3")
    {
        Rcpp::Rcout << "adding gff3 annotation: " << gff3_fn << "\n";
        Anno.parse_gff3_annotation(gff3_fn, fix_chrname);
    }
    else
    {
        Anno.parse_bed_annotation(gff3_fn, fix_chrname);
        Rcpp::Rcout << "adding bed annotation: " << gff3_fn << "\n";
    }

}

int Mapping::map_exon(bam_hdr_t *header, bam1_t *b, string& gene_id, bool m_strand)
{
    int ret = 9999;
    int rev = bam_is_rev(b)?(-1):1;
    uint32_t* cig = bam_get_cigar(b);
    int tmp_pos = b->core.pos;
    int tmp_rest = 9999999; // distance to end pos
    int tmp_ret;
    string tmp_id;
    gene_id = "";

    for (int c=0; c<b->core.n_cigar; c++)
    {
        tmp_ret = 9999;
        // *   bit 1 set if the cigar operation consumes the query
        // *   bit 2 set if the cigar operation consumes the reference
        if (((bam_cigar_type(cig[c]) >> 0) & 1) && ((bam_cigar_type(cig[c]) >> 1) & 1))
        {
            Interval it = Interval(tmp_pos, tmp_pos+bam_cigar_oplen(cig[c]), rev);
            auto &gene_list = Anno.gene_dict[header->target_name[b->core.tid]];
            auto iter = std::equal_range(gene_list.begin(), gene_list.end(), it);
            if ((iter.second - iter.first) == 0)
            {
                tmp_ret = tmp_ret>3?3:tmp_ret;
            }
            else
            {
                tmp_id = "";
                for (auto i = iter.first; i < iter.second; ++i)
                {
                    if (i->in_exon(it, m_strand))
                    {
                        if (tmp_id != "")
                        {
                            if (tmp_id != i->gene_id)
                            {
                                tmp_ret = 1; // ambiguous mapping
                                break;
                            }
                            else
                            {
                                // update the distance to end pos
                                tmp_rest = tmp_rest<(i->distance_to_end(it))?tmp_rest:i->distance_to_end(it);
                            }
                        }
                        else
                        {
                            tmp_id = i->gene_id;
                            tmp_ret = 0;
                            tmp_rest = i->distance_to_end(it);
                        }
                    }
                    else if ((it > *i) || (it < *i))
                    {
                        tmp_ret = tmp_ret>=3?3:tmp_ret;
                    }
                    else
                    {
                        tmp_ret = tmp_ret>=2?2:tmp_ret;
                    }
                }
            }

            tmp_pos = tmp_pos+bam_cigar_oplen(cig[c]);
            if (ret == 0 && tmp_ret == 0)
            {
                if (gene_id != "" && gene_id != tmp_id)
                {
                    ret = 1; // still ambiguous
                    break;
                }
            }
            else if (tmp_ret == 0)
            {
                ret = 0;
                gene_id = tmp_id;
            }
            else
            {
                ret = ret<tmp_ret?ret:tmp_ret; // choose the smallest
            }
        }
        else if (!((bam_cigar_type(cig[c]) >> 0) & 1) && ((bam_cigar_type(cig[c]) >> 1) & 1))
        {
            tmp_pos = tmp_pos+bam_cigar_oplen(cig[c]);
        }
    }
    if (ret == 0)
    {
        return -tmp_rest;
    }
    else
    {
        // return codes:
        // 0: unique map to exon
        // 1: ambiguous map to multiple exon
        // 2: map to intron
        // 3: not mapped
        return ret;
    }
}

// anonymous namespace for functions only used in this file
namespace {
    void report_every_3_mins(std::atomic<unsigned long long> &cnt, std::atomic<bool> &running) {
        Timer timer;
        timer.start();

        do {
            // sleep thread for a total of 3 minutes (180 seconds)
            // wake up at shorter intervals to check if process has stopped running
            for (int i = 0; i < 36; i++)
            {
                std::this_thread::sleep_for(std::chrono::seconds(5));
                if (!running)
                {
                    break;
                }
            }

            std::cout
                << cnt << " reads processed" << ", "
                << cnt / timer.seconds_elapsed() / 1000 << "k reads/sec" << "\n";
        } while (running);
    }
}

void Mapping::parse_align(string fn, string fn_out, bool m_strand, string map_tag, string gene_tag, string cellular_tag, string molecular_tag, int bc_len, int UMI_len)
{
    int unaligned = 0;
    int ret;

    check_file_exists(fn); // htslib does not check if file exist so we do it manually

    // open files
    bam1_t *b = bam_init1();
    BGZF *fp = bgzf_open(fn.c_str(), "r"); // input file
    samFile *of = sam_open(fn_out.c_str(), "wb"); // output file

    bam_hdr_t *header = bam_hdr_read(fp);
    sam_hdr_write(of, header);

    string gene_id;
    int tmp_c[4] = {0,0,0,0};

    bool found_any = false;
    for (int i = 0; i < header->n_targets; ++i)
    {
        if (Anno.gene_dict.end() == Anno.gene_dict.find(header->target_name[i]))
        {
            Rcpp::Rcout << header->target_name[i] << " not found in exon annotation." << "\n";
        }
        else
        {
            found_any = true;
        }

    }
    if (!found_any)
    {
        std::stringstream err_msg;
        err_msg << "ERROR: The annotation and .bam file contains different chromosome." << "\n";
        Rcpp::stop(err_msg.str());
    }
    // for moving barcode and UMI from sequence name to bam tags
    const char * g_ptr = gene_tag.c_str();
    const char * c_ptr = cellular_tag.c_str();
    const char * m_ptr = molecular_tag.c_str();
    const char * a_ptr = map_tag.c_str();
    char buf[999] = ""; // assume the length of barcode or UMI is less than 999
    std::atomic<unsigned long long> cnt{0};
    std::atomic<bool> running{true};

    Rcout << "annotating exon features..." << "\n";
    Rcout << "updating progress every 3 minutes..." << "\n";
    // spawn thread to report progress every 3 minutes
    std::thread reporter_thread(
        [&cnt, &running]() {
            report_every_3_mins(cnt, running); 
        }
    );

    while (bam_read1(fp, b) >= 0)
    {
        if (__DEBUG)
        {
            if (cnt % 1000000 == 0)
            {
                Rcpp::Rcout << "number of read processed:" << cnt << "\n";
                Rcpp::Rcout << tmp_c[0] <<"\t"<< tmp_c[1] <<"\t"<<tmp_c[2] <<"\t"<<tmp_c[3] <<"\t" << "\n";
            }
        }
        cnt++;

        if ((b->core.flag&BAM_FUNMAP) > 0)
        {
            unaligned++;
            ret = 4;
        }
        else
        {
            //  chromosome not found in annotation:
            if (Anno.gene_dict.end() == Anno.gene_dict.find(header->target_name[b->core.tid]))
            {
                ret = 3;
            }
            else
            {
                ret = map_exon(header, b, gene_id, m_strand);
            }
            if (ret <= 0)
            {
                tmp_c[0]++;
                bam_aux_append(b, g_ptr, 'Z', gene_id.size()+1, (uint8_t*)gene_id.c_str());
            }
            else
            {
                tmp_c[ret]++;
            }
        }
        if (bc_len > 0)
        {
            memcpy(buf, bam_get_qname(b), bc_len * sizeof(char));
            buf[bc_len] = '\0';
            bam_aux_append(b, c_ptr, 'Z', bc_len+1, (uint8_t*)buf);
        }
        if (UMI_len > 0)
        {
            memcpy(buf, bam_get_qname(b)+bc_len+1, UMI_len * sizeof(char)); // `+1` to add separator
            buf[UMI_len] = '\0';
            bam_aux_append(b, m_ptr, 'Z', UMI_len+1, (uint8_t*)buf);
        }

        bam_aux_append(b, a_ptr, 'i', sizeof(uint32_t), (uint8_t*)&ret);

        int re = sam_write1(of, header, b);
        if (re < 0)
        {
            std::stringstream err_msg;
            err_msg << "fail to write the bam file: " << bam_get_qname(b) << "\n";
            err_msg << "return code: " << re << "\n";
            Rcpp::stop(err_msg.str());
        }
    }

    running = false;
    reporter_thread.join();

    Rcpp::Rcout << "\t" << "number of read processed:" << cnt << "\n";
    Rcpp::Rcout << "\t" << "unique map to exon:" << tmp_c[0] << "\n";
    Rcpp::Rcout << "\t" << "ambiguous map to multiple exon:" << tmp_c[1] << "\n";
    Rcpp::Rcout << "\t" << "map to intron:" << tmp_c[2] << "\n";
    Rcpp::Rcout << "\t" << "not mapped:" << tmp_c[3] << "\n";
    Rcpp::Rcout << "\t" << "unaligned:" << unaligned << "\n";
    sam_close(of);
    bgzf_close(fp);
}

