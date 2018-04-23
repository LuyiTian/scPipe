// transcriptmapping.h
// 

#include <algorithm>
#include <atomic>
#include <cstring>
#include <fstream>
#include <iostream>
#include <Rcpp.h>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include "config_hts.h"
#include "utils.h"
#include "Gene.h"
#include "Interval.h"
#include "Timer.h"

#ifndef TRANSCRIPTMAPPING_H
#define TRANSCRIPTMAPPING_H

// parse gff3 genome annotation
class GeneAnnotation
{
public:
    std::string anno_source = "";
    std::unordered_set<std::string> recorded_genes;

    std::unordered_map<std::string, std::vector<Gene>> gene_dict;

    //get number of genes
    int ngenes();

    //return all gene id as a std::vector
    std::vector<std::string> get_genelist();

    void parse_gff3_annotation(std::string gff3_fn, bool fix_chrname);

    // https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    // bed file specification:
    // 1. chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
    // 2. chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    // 3. chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
    // 4. name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
    // 5. score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray).
    // 6. strand - Defines the strand - either '+' or '-'.
    // 7. thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
    // 8. thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
    // 9. itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
    // 10. blockCount - The number of blocks (exons) in the BED line.
    // 11. blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
    // 12. blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
    void parse_bed_annotation(std::string bed_fn, bool fix_chrname);

    friend std::ostream& operator<< (std::ostream& out, const GeneAnnotation& obj); 
    
    private:
    // index variables for gff3 fields
    const int SEQID = 0;
    const int SOURCE = 1;
    const int TYPE = 2;
    const int START = 3;
    const int END = 4;
    const int SCORE = 5;
    const int STRAND = 6;
    const int PHASE = 7;
    const int ATTRIBUTES  = 8;

    std::string get_attribute(const std::vector<std::string> &all_attributes, const std::string &target_attribute) {
        for (const std::string &attr : all_attributes) {
            auto sep_loc = attr.find("=");
            std::string key = attr.substr(0, sep_loc);
            std::string val = attr.substr(sep_loc + 1);
            if (key == target_attribute) {
                return val;
            }
        }
        return "";
    }

    int get_strand(char st)
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

    const std::string get_parent(const std::vector<std::string> &attributes)
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

    std::string get_ID(const std::vector<std::string> &attributes)
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

    std::string fix_name(std::string chr_name)
    {
        std::string new_chr_name;
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

    void parse_anno_entry(const bool &fix_chrname, const std::string &line, std::unordered_map<std::string, std::unordered_map<std::string, Gene>> &chr_to_genes_dict, std::unordered_map<std::string, std::string> &transcript_to_gene_dict)
    {
        const std::vector<std::string> fields = split(line, '\t');
        const std::vector<std::string> attributes = split(fields[ATTRIBUTES], ';');

        std::string chr_name = fields[SEQID];
        const std::string parent = get_parent(attributes);
        const std::string type = fields[TYPE];
        const std::string ID = get_ID(attributes);
        const int strand = get_strand(fields[STRAND][0]);
        const int interval_start = std::atoi(fields[START].c_str());
        const int interval_end = std::atoi(fields[END].c_str());

        if (fix_chrname)
        {
            chr_name = fix_name(chr_name);
        }

        // DEBUG USE
        // Rcpp::Rcout << "Parsing: " << line << "\n";
        // Rcpp::Rcout << "Type: " << type << " "
        //       << "ID: " << ID << " "
        //       << "Parent: " << parent << "\n\n";
        // DEBUG USE

        std::string target_gene;
        if (anno_source == "ensembl")
        {
            if (is_gene(fields, attributes)) {
                recorded_genes.insert(ID);
                return;
            }
            else if (is_transcript(fields, attributes))
            {
                if (!ID.empty() && !parent.empty())
                {
                    transcript_to_gene_dict[ID] = parent;
                }
                return;
            }
            else if (is_exon(fields, attributes))
            {
                if (parent_is_known_transcript(transcript_to_gene_dict, parent))
                {
                    target_gene = transcript_to_gene_dict[parent];
                }
                else
                {
                    std::stringstream err_msg;
                    err_msg << "cannot find grandparent for exon:" << "\n";
                    err_msg << line << "\n";
                    Rcpp::stop(err_msg.str());
                }
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

    std::string get_gencode_gene_id(const std::vector<std::string> &attributes)
    {
        return get_attribute(attributes, "gene_id");
    }

    std::string get_refseq_gene_id(const std::vector<std::string> &attributes)
    {
        std::string dbxref = get_attribute(attributes, "Dbxref");

        // GeneID may be missing
        if (dbxref.find("GeneID") == std::string::npos)
        {
            return "";
        }
        
        auto start = dbxref.find("GeneID") + 7; // start after "GeneID:"
	    auto end = dbxref.find(",", start);
        auto id_length = end - start;

        return dbxref.substr(start, id_length);
    }

    std::string guess_anno_source(std::string gff3_fn) {
        std::ifstream infile(gff3_fn);
        std::string line;

        while (std::getline(infile, line))
        {
            if (line.find("GENCODE") != std::string::npos) {
                Rcpp::Rcout << "guessing annotation source: GENCODE" << "\n";
                return "gencode";
            }
            else if (line.find("1\tEnsembl") != std::string::npos)
            {
                Rcpp::Rcout << "guessing annotation source: ENSEMBL" << "\n";
                return "ensembl";
            }
            else if (line.find("RefSeq\tregion") != std::string::npos)
            {
                Rcpp::Rcout << "guessing annotation source: RefSeq" << "\n";
                return "refseq";
            }
        }

        Rcpp::stop("Annotation source not recognised. Current supported sources: ENSEMBL, GENCODE and RefSeq");
    }

    std::string get_gene_id(const std::vector<std::string> &attributes) {
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

    const bool parent_is_gene(const std::string &parent)
    {
        return recorded_genes.find(parent) != recorded_genes.end();
    }

    const bool parent_is_known_transcript(const std::unordered_map<std::string, std::string> &transcript_to_gene_dict, const std::string &parent)
    {
        return transcript_to_gene_dict.find(parent) != transcript_to_gene_dict.end();
    }

    const bool is_gene(const std::vector<std::string> &fields, const std::vector<std::string> &attributes)
    {
        std::string type = fields[TYPE];
        if (type == "gene")
        {
            return true;
        }

        std::string id = get_attribute(attributes, "ID");
        if (id.find("gene:") != std::string::npos)
        {
            return true;
        }

        return false;
    }

    const bool is_exon(const std::vector<std::string> &fields, const std::vector<std::string> &attributes) {
        return fields[TYPE] == "exon";
    }

    const bool is_transcript(const std::vector<std::string> &fields, const std::vector<std::string> &attributes) {
        // assume feature is transcript is it has a gene as parent
        return parent_is_gene(get_parent(attributes));
    }
};


class Mapping
{
public:
    GeneAnnotation Anno;
    void add_annotation(std::string gff3_fn, bool fix_chrname);
    // return:
    //  <=0 - unique map to exon, number indicate the distance to transcript end pos
    //  1 - ambiguous map to multiple exon
    //  2 - map to intron
    //  3 - unmapped
    //  4 - unaligned
    int map_exon(bam_hdr_t *header, bam1_t *b, std::string& gene_id, bool m_strand);


    // @param: m_strand, match based on strand or not
    // @param: fn_out, output bam file
    // @param: 
    void parse_align(std::string fn, std::string fn_out, bool m_strand, std::string map_tag, std::string gene_tag, std::string cellular_tag, std::string molecular_tag, int bc_len, int UMI_len);

    /* data */
};

#endif
