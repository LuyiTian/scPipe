// transcriptmapping.h
// 

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <cstring>
#include <Rcpp.h>
#include "config_hts.h"
#include "utils.h"
#include "Gene.h"
#include "Interval.h"

#ifndef TRANSCRIPTMAPPING_H
#define TRANSCRIPTMAPPING_H

// parse gff3 genome annotation
class GeneAnnotation
{
public:
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
