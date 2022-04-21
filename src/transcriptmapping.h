// transcriptmapping.h
//

#include <algorithm>
#include <atomic>
#include <cstring>
#include <fstream>
#include <iostream>
#include <Rcpp.h>
#include <regex>
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

class GeneBin {
using ull_int = unsigned long long;
public:
    // class data
    ull_int start;
    ull_int end = 0;
    std::vector<Gene> genes;

    // add gene to bin
    void add_gene(Gene gene)
    {
        if (genes.size() == 0) {
            start = gene.st;
        }
        genes.push_back(gene);
        // extend boundaries of bin to fully include bin
        if ((ull_int) gene.st < start)
        {
            start = gene.st;
        }
        if ((ull_int) gene.en > end)
        {
            end = gene.en;
        }
    }

    // check if interval overlaps bin
    const bool overlaps(const Interval &it) {
        return !(start > (ull_int) it.en) && !(end < (ull_int) it.st);
    }
};

class GeneBins {
public:
    // class data
    std::vector<GeneBin> gene_bins;

    // get vector of bins overlapping query interval
    std::vector<GeneBin*> get_bins(Interval it)
    {
        std::vector<GeneBin*> overlapped_bins;
        for (auto &bin : gene_bins) {
            if (bin.overlaps(it)) {
                overlapped_bins.push_back(&bin);
            }
        }
        return overlapped_bins;
    }

    // create bins from vector of genes
    void make_bins(std::vector<Gene> &genes) {
        unsigned int bin_index = 0;
        unsigned int count = 0;

        for (auto gene : genes) {
            if (bin_index + 1 > gene_bins.size()) {
                // add bin if required
                gene_bins.resize(gene_bins.size() + 1);
            }

            gene_bins[bin_index].add_gene(gene);
            count++;

            if (count == bin_size) {
                // bin is full, move to next bin
                count = 0;
                bin_index++;
            }
        }
    }
private:
    const unsigned int bin_size = 64;
};

enum AnnotationEnum {
    SEQID=0,
    SOURCE=1,
    TYPE=2,
    START=3,
    END=4,
    SCORE=5,
    STRAND=6,
    PHASE=7,
    ATTRIBUTES=8
};

// parse gff3 genome annotation
class GeneAnnotation
{
public:
    std::string anno_source = "";
    std::unordered_set<std::string> recorded_genes;

    std::unordered_map<std::string, std::vector<Gene>> gene_dict;
    std::unordered_map<std::string, GeneBins> bins_dict;

    //get number of genes
    int ngenes();

    //return all gene id as a std::vector
    std::vector<std::string> get_genelist();

    void parse_gff3_annotation(std::string gff3_fn, bool fix_chrname);
    void parse_saf_dataframe(Rcpp::DataFrame anno_df, bool fix_chrname);

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
    // // index variables for gff3 fields
    // const int SEQID      = 0;
    // const int SOURCE     = 1;
    // const int TYPE       = 2;
    // const int START      = 3;
    // const int END        = 4;
    // const int SCORE      = 5;
    // const int STRAND     = 6;
    // const int PHASE      = 7;
    // const int ATTRIBUTES = 8;

    // get attribute from gff3 standard columns
    std::string get_attribute(const std::vector<std::string> &all_attributes, const std::string &target_attribute);
    // convert strand from +- symbols to -1 or 1
    int get_strand(char st);

    const std::string get_parent(const std::vector<std::string> &attributes);
    std::string get_ID(const std::vector<std::string> &attributes);

    // add chr to molecule names if requested
    std::string fix_name(std::string chr_name);

    // parse entry of gff3 annotation and add its information to this object
    void parse_anno_entry(
        const bool &fix_chrname,
        const std::string &line,
        std::unordered_map<std::string, std::unordered_map<std::string, Gene>> &chr_to_genes_dict,
        std::unordered_map<std::string, std::string> &transcript_to_gene_dict
    );

    // generic gene_id getter for gff3 entries
    std::string get_gene_id(const std::vector<std::string> &attributes);

    // specific gene_id getter for gff3 entries
    std::string get_gencode_gene_id(const std::vector<std::string> &attributes);
    std::string get_refseq_gene_id(const std::vector<std::string> &attributes);

    // guess the source of annotation
    std::string guess_anno_source(std::string gff3_fn);

    const bool parent_is_gene(const std::string &parent);
    const bool parent_is_known_transcript(const std::unordered_map<std::string, std::string> &transcript_to_gene_dict, const std::string &parent);
    const bool is_gene(const std::vector<std::string> &fields, const std::vector<std::string> &attributes);
    const bool is_exon(const std::vector<std::string> &fields, const std::vector<std::string> &attributes);
    const bool is_transcript(const std::vector<std::string> &fields, const std::vector<std::string> &attributes);
};


class Mapping
{
public:
    GeneAnnotation Anno;
    void add_annotation(std::string gff3_fn, bool fix_chrname);
    void add_annotation(Rcpp::DataFrame anno, bool fix_chrname);
    // return:
    //  <=0 - unique map to exon, number indicate the distance to transcript end pos
    //  1 - ambiguous map to multiple exon
    //  2 - map to intron
    //  3 - unmapped
    //  4 - unaligned
    int map_exon(bam_hdr_t *header, bam1_t *b, std::string& gene_id, bool m_strand);

    void parse_align_warpper(std::vector<std::string> fn_vec, std::vector<std::string> cell_id_vec, std::string fn_out, bool m_strand, std::string map_tag, std::string gene_tag, std::string cellular_tag, std::string molecular_tag, int bc_len, int UMI_len, int nthreads);
    // @param: m_strand, match based on strand or not
    // @param: fn_out, output bam file
    // @param: write_mode, whether to write (wb) or attach (ab) to a bam file.
    // @param: cell_id, to provide the cell barcode, if not in the header of fastq file.
    void parse_align(std::string fn, std::string fn_out, bool m_strand, std::string map_tag, std::string gene_tag, std::string cellular_tag, std::string molecular_tag, int bc_len, std::string write_mode, std::string cell_id, int UMI_len, int nthreads);

    void sc_atac_parse_align_warpper(std::vector<std::string> fn_vec, std::string fn_out,  std::string cellular_tag, std::string molecular_tag, int nthreads);
    void sc_atac_parse_align(std::string fn, std::string fn_out, std::string cellular_tag, std::string molecular_tag, int nthreads);
    
    /* data */
};

#endif
