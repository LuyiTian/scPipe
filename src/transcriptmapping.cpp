// transcriptmapping.cpp
#include "transcriptmapping.h"

using std::atoi;
using std::atomic;
using std::endl;
using std::fixed;
using std::getline;
using std::ifstream;
using std::ostream;
using std::setprecision;
using std::sort;
using std::string;
using std::stringstream;
using std::thread;
using std::unordered_map;
using std::vector;

using namespace std::this_thread;
using namespace std::chrono;
using namespace Rcpp;

string GeneAnnotation::get_attribute(
    const vector<string> &all_attributes,
    const string &target_attribute
)
{
    for (const string &attr : all_attributes) {
        auto sep_loc = attr.find("=");
        // get key
        const string key = attr.substr(0, sep_loc);
        // get value
        const string val = attr.substr(sep_loc + 1);
        if (key == target_attribute) {
            return val;
        }
    }
    return "";
}

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

string GeneAnnotation::get_ID(const vector<string> &attributes)
{
    for (const auto &attr : attributes)
    {
        if (attr.substr(0, 2) == "ID")
        {
            // check for ENSEMBL notation
            if (anno_source == "ensembl")
            {
                return attr.substr(attr.rfind(':') + 1);
            }
            else
            {
                return attr.substr(attr.find('=') + 1);
            }
        }
    }
    return "";
}

const string GeneAnnotation::get_parent(const vector<string> &attributes)
{
    for (const auto &attr : attributes)
    {
        if (attr.substr(0, 6) == "Parent")
        {
            // check for ENSEMBL notation
            if (anno_source == "ensembl")
            {
                return attr.substr(attr.rfind(':') + 1);
            }
            else
            {
                return attr.substr(attr.find('=') + 1);
            }
        }
    }
    return "";
}

string GeneAnnotation::fix_name(string chr_name)
{
    string new_chr_name;
    if (chr_name.compare(0, 3, "chr") == 0)
    {
        return chr_name;
    }
    else if (chr_name.length() > 4) // just fix 1-22, X, Y, MT. ignore contig and ERCC
    {
        return chr_name;
    }
    else
    {
        if (chr_name == "MT")
        {
            new_chr_name = "chrM";
        }
        else
        {
            new_chr_name = "chr" + chr_name;
        }
        return new_chr_name;
    }
}

string GeneAnnotation::get_gene_id(const vector<string> &attributes)
{
    if (anno_source == "gencode")
    {
        return get_gencode_gene_id(attributes);
    }
    else if (anno_source == "refseq")
    {
        return get_refseq_gene_id(attributes);
    }
    return "";
}

string GeneAnnotation::get_gencode_gene_id(const vector<string> &attributes)
{
    return get_attribute(attributes, "gene_id");
}

string GeneAnnotation::get_refseq_gene_id(const vector<string> &attributes)
{
    string dbxref = get_attribute(attributes, "Dbxref");

    // GeneID may be missing
    if (dbxref.find("GeneID") == string::npos)
    {
        return "";
    }

    auto start = dbxref.find("GeneID") + 7; // start after "GeneID:"
    auto end = dbxref.find(",", start);
    auto id_length = end - start;

    return dbxref.substr(start, id_length);
}

void GeneAnnotation::parse_anno_entry(
    const bool &fix_chrname,
    const string &line,
    unordered_map<string, unordered_map<string, Gene>> &chr_to_genes_dict,
    unordered_map<string, string> &transcript_to_gene_dict
)
{
    const vector<string> fields = split(line, '\t');
    const vector<string> attributes = split(fields[ATTRIBUTES], ';');

    string chr_name = fields[SEQID];
    const string parent = get_parent(attributes);
    const string type = fields[TYPE];
    const string ID = get_ID(attributes);
    const int strand = get_strand(fields[STRAND][0]);
    const int interval_start = atoi(fields[START].c_str());
    const int interval_end = atoi(fields[END].c_str());

    if (fix_chrname)
    {
        chr_name = fix_name(chr_name);
    }

    // DEBUG USE
    // Rcout << "Parsing: " << line << "\n";
    // Rcout << "Type: " << type << " "
    //       << "ID: " << ID << " "
    //       << "Parent: " << parent << "\n\n";
    // DEBUG USE

    string target_gene;
    if (anno_source == "ensembl")
    {
        if (is_exon(fields, attributes))
        {
            if (parent_is_known_transcript(transcript_to_gene_dict, parent))
            {
                target_gene = transcript_to_gene_dict[parent];
            }
            else
            {
                stringstream err_msg;
                err_msg << "cannot find grandparent for exon:" << "\n";
                err_msg << line << "\n";
                stop(err_msg.str());
            }
        }
        else if (is_transcript(fields, attributes))
        {
            if (!ID.empty() && !parent.empty())
            {
                transcript_to_gene_dict[ID] = parent;
            }
            return;
        }
        else if (is_gene(fields, attributes)) {
            recorded_genes.insert(ID);
            return;
        }
    }
    else if (anno_source == "gencode" || anno_source == "refseq")
    {
        if (type == "exon")
        {
            target_gene = get_gene_id(attributes);
        }
    }

    if (!target_gene.empty())
    {
        auto &current_chr = chr_to_genes_dict[chr_name];
        current_chr[target_gene].add_exon(Interval(interval_start, interval_end, strand));
        current_chr[target_gene].set_ID(target_gene);
    }

    return;
}

string GeneAnnotation::guess_anno_source(string gff3_fn)
{
    ifstream infile(gff3_fn);
    string line;

    while (getline(infile, line))
    {
        if (line.find("GENCODE") != string::npos) {
            Rcout << "guessing annotation source: GENCODE" << "\n";
            return "gencode";
        }
        else if (line.find("1\tEnsembl") != string::npos)
        {
            Rcout << "guessing annotation source: ENSEMBL" << "\n";
            return "ensembl";
        }
        else if (line.find("RefSeq\tregion") != string::npos)
        {
            Rcout << "guessing annotation source: RefSeq" << "\n";
            return "refseq";
        }
    }

    Rcout << "Annotation source not recognised, defaulting to ENSEMBL. Current supported sources: ENSEMBL, GENCODE and RefSeq\n";
    return "ensembl";
}

const bool GeneAnnotation::parent_is_gene(const string &parent)
{
    return recorded_genes.find(parent) != recorded_genes.end();
}

const bool GeneAnnotation::parent_is_known_transcript(const unordered_map<string, string> &transcript_to_gene_dict, const string &parent)
{
    return transcript_to_gene_dict.find(parent) != transcript_to_gene_dict.end();
}

const bool GeneAnnotation::is_gene(const vector<string> &fields, const vector<string> &attributes)
{
    string type = fields[TYPE];
    if (type.find("gene") != string::npos)
    {
        return true;
    }

    string id = get_attribute(attributes, "ID");
    if (id.find("gene:") != string::npos)
    {
        return true;
    }

    return false;
}

const bool GeneAnnotation::is_exon(const vector<string> &fields, const vector<string> &attributes)
{
    return fields[TYPE] == "exon";
}

const bool GeneAnnotation::is_transcript(const vector<string> &fields, const vector<string> &attributes)
{
    // assume feature is transcript is it has a gene as parent
    return parent_is_gene(get_parent(attributes));
}

void GeneAnnotation::parse_gff3_annotation(string gff3_fn, bool fix_chrname)
{
    ifstream infile(gff3_fn);

    string line;
    unordered_map<string, unordered_map<string, Gene>> chr_to_genes_dict;
    unordered_map<string, string> transcript_to_gene_dict; // store transcript - gene mapping

    // assigned to class member
    anno_source = guess_anno_source(gff3_fn);

    size_t _interrupt_ind = 0;
    // create transcript-gene mapping
    while (getline(infile, line))
    {
        if (++_interrupt_ind % 256 == 0) checkUserInterrupt();
        // skip header lines
        if (line[0] == '#')
        {
            continue;
        }

        parse_anno_entry(fix_chrname, line, chr_to_genes_dict, transcript_to_gene_dict);
    }

    // push genes into annotation class member
    for (auto &chr : chr_to_genes_dict)
    {
        const auto &chr_name = chr.first;

        // merge overlapping exons in each gene
        for (auto &gene : chr.second)
        {
            gene.second.sort_exon();
            gene.second.flatten_exon();
            gene_dict[chr_name].push_back(gene.second);
        }

        auto &current_genes = gene_dict[chr_name];
        // sort genes based on starting position
        sort(current_genes.begin(), current_genes.end(),
            [] (const Gene &g1, const Gene &g2) { return g1.st < g2.st; }
        );

        // create bins of genes
        bins_dict[chr_name].make_bins(current_genes);
    }
}

void GeneAnnotation::parse_bed_annotation(string bed_fn, bool fix_chrname)
{
    ifstream infile(bed_fn);

    string line;
    unordered_map<string, unordered_map<string, Gene>> tmp_gene_dict;
    int strand = 0;
    vector<string> token;

    getline(infile, line); // skip the header
    while(getline(infile, line))
    {
        token = split(line, '\t');
        strand = get_strand(token[4][0]);
        if (fix_chrname)
        {
            tmp_gene_dict[fix_name(token[1])][token[0]].add_exon(Interval(atoi(token[2].c_str()), atoi(token[3].c_str()), strand));
            tmp_gene_dict[fix_name(token[1])][token[0]].set_ID(token[1]);
        }
        else
        {
            tmp_gene_dict[token[1]][token[0]].add_exon(Interval(atoi(token[2].c_str()), atoi(token[3].c_str()), strand));
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
                sub_iter.second.flatten_exon();
            }
            gene_dict[iter.first].push_back(sub_iter.second);
        }
        if (gene_dict[iter.first].size()>1)
        {
            sort(gene_dict[iter.first].begin(), gene_dict[iter.first].end(),
                [] (const Gene &g1, const Gene &g2) { return g1.st < g2.st; }
            );
        }

        // create bins of genes
        bins_dict[iter.first].make_bins(gene_dict[iter.first]);
    }
}

void GeneAnnotation::parse_saf_dataframe(DataFrame anno_df, bool fix_chrname)
{
    CharacterVector gene_ids = anno_df["GeneID"];
    CharacterVector chrs = anno_df["Chr"];
    NumericVector starts = anno_df["Start"];
    NumericVector ends = anno_df["End"];
    CharacterVector strands = anno_df["Strand"];

    int n_entries = gene_ids.size();

    using Chr = std::string;
    using GeneID = std::string;
    unordered_map<Chr, unordered_map<GeneID, Gene>> tmp_gene_dict;
    for (int i = 0; i < n_entries; i++)
    {
        std::string const &gene_id = as<std::string>(gene_ids[i]);
        std::string const &chr = as<std::string>(chrs[i]);
        int start = starts[i];
        int end = ends[i];
        int const &strand = strands[i] == "+" ? 1 :
                            strands[i] == "-" ? -1 : 0;

        if (fix_chrname)
        {
            tmp_gene_dict[fix_name(chr)][gene_id].add_exon(Interval(start, end, strand));
            tmp_gene_dict[fix_name(chr)][gene_id].set_ID(chr);
        }
        else
        {
            tmp_gene_dict[chr][gene_id].add_exon(Interval(start, end, strand));
            tmp_gene_dict[chr][gene_id].set_ID(gene_id);
        }

    }

    for (auto iter : tmp_gene_dict)
    {
        for (auto sub_iter : iter.second)
        {
            if (sub_iter.second.exon_vec.size()>1)
            {
                sub_iter.second.sort_exon();
                sub_iter.second.flatten_exon();
            }
            gene_dict[iter.first].push_back(sub_iter.second);
        }
        if (gene_dict[iter.first].size()>1)
        {
            sort(gene_dict[iter.first].begin(), gene_dict[iter.first].end(),
                [] (const Gene &g1, const Gene &g2) { return g1.st < g2.st; }
            );
        }

        // create bins of genes
        bins_dict[iter.first].make_bins(gene_dict[iter.first]);
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


vector<string> GeneAnnotation::get_genelist()
{
    vector<string> gene_list;
    for (auto iter : gene_dict)
    {
        for (auto sub_iter : iter.second)
        {
            gene_list.push_back(sub_iter.gene_id);
        }
    }

    return gene_list;
}


ostream& operator<< (ostream& out, const GeneAnnotation& obj)
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
    if (gff3_fn.substr(gff3_fn.find_last_of(".")) == ".gff3" ||
        gff3_fn.substr(gff3_fn.find_last_of(".")) == ".gff")
    {
        Rcout << "adding gff3 annotation: " << gff3_fn << "\n";
        Anno.parse_gff3_annotation(gff3_fn, fix_chrname);
    }
    else
    {
        Anno.parse_bed_annotation(gff3_fn, fix_chrname);
        Rcout << "adding bed annotation: " << gff3_fn << "\n";
    }
}

void Mapping::add_annotation(DataFrame anno, bool fix_chrname)
{
    Anno.parse_saf_dataframe(anno, fix_chrname);
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
        string chr_name{header->target_name[b->core.tid]};
        // *   bit 1 set if the cigar operation consumes the query
        // *   bit 2 set if the cigar operation consumes the reference
        const bool consumes_qry = (bam_cigar_type(cig[c]) >> 0) & 1;
        const bool consumes_ref = (bam_cigar_type(cig[c]) >> 1) & 1;
        if (consumes_qry && consumes_ref)
        {
            Interval it = Interval(tmp_pos, tmp_pos+bam_cigar_oplen(cig[c]), rev);
            auto &bins_list = Anno.bins_dict[chr_name];
            const vector<GeneBin*> &matched_gene_bins = bins_list.get_bins(it);

            vector<Gene> matched_genes;

            for (auto &gene_list_ptr : matched_gene_bins) {
                for (auto &gene : gene_list_ptr->genes) {
                    if (gene == it) {
                        matched_genes.push_back(gene);
                    }
                }
            }

            if (matched_genes.size() == 0)
            {
                // no matching gene
                tmp_ret = (tmp_ret > 3) ? 3 : tmp_ret;
            }
            else
            {
                tmp_id = "";
                for (auto &gene : matched_genes)
                {
                    if (gene.in_exon(it, m_strand))
                    {
                        if (tmp_id != "")
                        {
                            if (tmp_id != gene.gene_id)
                            {
                                tmp_ret = 1; // ambiguous mapping
                                break;
                            }
                            else
                            {
                                // update the distance to end pos
                                tmp_rest = tmp_rest<(gene.distance_to_end(it))?tmp_rest:gene.distance_to_end(it);
                            }
                        }
                        else
                        {
                            tmp_id = gene.gene_id;
                            tmp_ret = 0;
                            tmp_rest = gene.distance_to_end(it);
                        }
                    }
                    else if ((it > gene) || (it < gene))
                    {
                        tmp_ret = (tmp_ret >= 3) ? 3 : tmp_ret;
                    }
                    else
                    {
                        tmp_ret = (tmp_ret >= 2) ? 2 : tmp_ret;
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
        else if (!consumes_qry && consumes_ref)
        {
            tmp_pos += bam_cigar_oplen(cig[c]);
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

namespace {
    void report_every_3_mins(
        atomic<unsigned long long> &cnt,
        atomic<bool> &running,
        atomic<bool> &report_message
    )
    {
        do {
            // sleep thread for a total of 3 minutes (180 seconds)
            // wake up at shorter intervals to check if process has stopped running
            for (int i = 0; i < 180; i++)
            {
                sleep_for(seconds(1));
                if (!running)
                {
                    break;
                }
            }

            report_message = true;
        } while (running);
    }
}

void Mapping::parse_align_warpper(vector<string> fn_vec, vector<string> cell_id_vec, string fn_out, bool m_strand, string map_tag, string gene_tag, string cellular_tag, string molecular_tag, int bc_len, int UMI_len, int nthreads)
{
  if (fn_vec.size()>1)
  {
    if((fn_vec.size()!=cell_id_vec.size()) & (bc_len==0))
    {
      stringstream err_msg;
      err_msg << "size of bam file and cell id vector should be the same: \n";
      err_msg << "\t number of bam files: " << fn_vec.size() << "\n";
      err_msg << "\t number of cell ids: " << cell_id_vec.size() << "\n";
      stop(err_msg.str());
    }
    if (bc_len==0)
    {
      parse_align(fn_vec[0], fn_out, m_strand, map_tag, gene_tag, cellular_tag, molecular_tag, bc_len, "wb", cell_id_vec[0], UMI_len, nthreads);
      for (int i=1;i<fn_vec.size();i++)
      {
        parse_align(fn_vec[i], fn_out, m_strand, map_tag, gene_tag, cellular_tag, molecular_tag, bc_len, "ab", cell_id_vec[i], UMI_len, nthreads);
      }
    }
    else
    {
      parse_align(fn_vec[0], fn_out, m_strand, map_tag, gene_tag, cellular_tag, molecular_tag, bc_len, "wb", "", UMI_len, nthreads);
      for (int i=1;i<fn_vec.size();i++)
      {
        parse_align(fn_vec[i], fn_out, m_strand, map_tag, gene_tag, cellular_tag, molecular_tag, bc_len, "ab", "", UMI_len, nthreads);
      }
    }
  }
  else
  {
    parse_align(fn_vec[0], fn_out, m_strand, map_tag, gene_tag, cellular_tag, molecular_tag, bc_len, "wb", "", UMI_len, nthreads);
  }
}

namespace {
    std::pair<int, int> get_bc_umi_lengths(string bam_fn) {
        BGZF *fp = bgzf_open(bam_fn.c_str(), "r"); // input file
        bam_hdr_t *bam_hdr = bam_hdr_read(fp);

        bam1_t *bam_record = bam_init1();

        if (bam_read1(fp, bam_record) >= 0) {
            string read_header = bam_get_qname(bam_record);
            int break_pos = read_header.find("#");
            // start from 1 to exclude @
            string first_section = read_header.substr(1, break_pos);
            bool valid_pattern = std::regex_search(first_section, std::regex("[ACTGN]+_[ACRGN]+"));
            if (!valid_pattern) {
                throw std::runtime_error("Read header does not contain valid barcode-umi data.");
            }

            // header has structure {BARCODE}_{UMI}
            int bc_len = first_section.find("_") + 1;
            int umi_len = first_section.length() - bc_len - 1;

            std::cout << "detected barcode length: " << bc_len << "\n";
            std::cout << "detected UMI length: " << umi_len << "\n";

            return std::make_pair(bc_len, umi_len);
        }
        else
        {
            throw std::runtime_error("BAM file reading failed.");
        }
    }
}

void Mapping::parse_align(string bam_fn, string fn_out, bool m_strand, string map_tag, string gene_tag, string cellular_tag, string molecular_tag, int bc_len, string write_mode, string cell_id, int UMI_len, int nthreads)
{
    int unaligned = 0;
    int ret;

    check_file_exists(bam_fn); // htslib does not check if file exist so we do it manually

    // FUTURE: guess BC and UMI lengths
    // int bc_len;
    // int UMI_len;
    // std::tie(bc_len, UMI_len) = get_bc_umi_lengths(bam_fn);

    const char * c_write_mode = write_mode.c_str();
    // open files
    bam1_t *b = bam_init1();
    BGZF *fp = bgzf_open(bam_fn.c_str(), "r"); // input file
    samFile *of = sam_open(fn_out.c_str(), c_write_mode); // output file

    // set up htslib threadpool for output
    int out_threads = std::max(nthreads - 1, 1);
    htsThreadPool p = {NULL, 0};
    p.pool = hts_tpool_init(out_threads);
    hts_set_opt(of, HTS_OPT_THREAD_POOL, &p);

    int hts_retcode;

    bam_hdr_t *header = bam_hdr_read(fp);
    hts_retcode = sam_hdr_write(of, header);

    int tmp_c[4] = {0,0,0,0};

    bool found_any = false;
    for (int i = 0; i < header->n_targets; ++i)
    {
        if (Anno.gene_dict.end() == Anno.gene_dict.find(header->target_name[i]))
        {
            Rcout << header->target_name[i] << " not found in exon annotation." << "\n";
        }
        else
        {
            found_any = true;
        }

    }
    if (!found_any)
    {
        stringstream err_msg;
        err_msg << "ERROR: The annotation and .bam file contains different chromosome." << "\n";
        stop(err_msg.str());
    }
    uint8_t* c_cell_id = const_cast<uint8_t*>(reinterpret_cast<const uint8_t*>(cell_id.c_str()));
    // for moving barcode and UMI from sequence name to bam tags
    const char * g_ptr = gene_tag.c_str();
    const char * c_ptr = cellular_tag.c_str();
    const char * m_ptr = molecular_tag.c_str();
    const char * a_ptr = map_tag.c_str();
    char buf[999] = ""; // assume the length of barcode or UMI is less than 999

    atomic<unsigned long long> cnt{0};
    atomic<bool> running{true};
    atomic<bool> report_message{false};

    Rcout << "updating progress every 3 minutes..." << "\n";
    // spawn thread to set report_message to true every 3 minutes
    thread reporter_thread(
        [&]() {
            report_every_3_mins(cnt, running, report_message);
        }
    );
    Timer timer;
    timer.start();

    while (bam_read1(fp, b) >= 0)
    {
        string gene_id;

        if (__DEBUG)
        {
            if (cnt % 1000000 == 0)
            {
                Rcout << "number of read processed:" << cnt << "\n";
                Rcout << tmp_c[0] <<"\t"<< tmp_c[1] <<"\t"<<tmp_c[2] <<"\t"<<tmp_c[3] <<"\t" << "\n";
            }
        }
        cnt++;
        if (cnt % 32768 == 0) checkUserInterrupt();

        // The Rcout would be conceptually cleaner if it lived inside the spawned thread
        // but only the master thread can interact with R without error so this code CANNOT
        // be run inside a child thread
        if (report_message)
        {
            Rcout
                << cnt << " reads processed" << ", "
                << cnt / timer.seconds_elapsed() / 1000 << "k reads/sec" << endl;
            report_message = false;
        }

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
                if (ret >= 0 && ret <= 3)
                    tmp_c[ret]++;
            }
        }
        if (bc_len > 0)
        {
            memcpy(buf, bam_get_qname(b), bc_len * sizeof(char));
            buf[bc_len] = '\0';
            bam_aux_append(b, c_ptr, 'Z', bc_len+1, (uint8_t*)buf);
        } else if (cell_id.size()>0)
        {
            bam_aux_append(b, c_ptr, 'Z', cell_id.size()+1, c_cell_id);
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
            stringstream err_msg;
            err_msg << "fail to write the bam file: " << bam_get_qname(b) << "\n";
            err_msg << "return code: " << re << "\n";
            stop(err_msg.str());
        }
    }

    running = false;
    reporter_thread.join();

    // final report of processing speed
    Rcout
        << cnt << " reads processed" << ", "
        << cnt / timer.seconds_elapsed() / 1000 << "k reads/sec" << endl;

    Rcout << "number of read processed: " << cnt << "\n";
    Rcout << "unique map to exon: " << tmp_c[0]
        << " (" << fixed << setprecision(2) << 100. * tmp_c[0]/cnt << "%)" << "\n";

    Rcout << "ambiguous map to multiple exon: " << tmp_c[1]
        << " ("  << fixed << setprecision(2) << 100. * tmp_c[1]/cnt << "%)" << "\n";

    Rcout << "map to intron: " << tmp_c[2]
        << " (" << fixed << setprecision(2) << 100. * tmp_c[2]/cnt << "%)" << "\n";

    Rcout << "not mapped: " << tmp_c[3]
        << " ("  << fixed << setprecision(2) << 100. * tmp_c[3]/cnt << "%)" << "\n";

    Rcout << "unaligned: " << unaligned
        << " (" << fixed << setprecision(2) << 100. * unaligned/cnt << "%)" << "\n";
    sam_close(of);
    bgzf_close(fp);
    if (p.pool) hts_tpool_destroy(p.pool);
}
