// fq simulator
#include <random>
#include <iostream>
#include <fstream>
//#include <math>  // floor
#include <chrono>  // high_resolution_clock::now()
#include <unordered_map>
#include "Eigen/Dense"
#include "config_hts.h"
#include "transcriptmapping.h"
#include "cellbarcode.h"
#include "utils.h"

#ifndef FQSIMULATOR_H
#define FQSIMULATOR_H

const std::vector<std::string> BP ({"A", "T", "G", "C"});

struct Fa_rec
{
    std::string name;
    std::string seq;
};

struct Fq_rec
{
    std::string name;
    std::string seq;
    std::string qual;
};


// simple fasta IO
class FaReader
{
public:
    std::ifstream fafile;
    std::string line;
    Fa_rec fa;
    FaReader(std::string fafn);
    bool readone();
};


// simple fastq IO
class FqWriter
{
public:
    std::ofstream fqfile;
    Fq_rec fq;
    FqWriter(std::string fqfn);
    void writeone();
};


// to simulate gene expression matrix
class CountSimulator
{
public:
    unsigned sim_seed;
    Eigen::MatrixXi cnt_mat;
    std::unordered_map<std::string, int> gene_ix;
    std::vector<std::string> ix_gene;
    std::unordered_map<std::string, int> cell_ix;
    std::vector<std::string> ix_cell;
    CountSimulator(unsigned seed);
    CountSimulator();
    void init_mat(std::vector<std::string> gene_v, std::vector<std::string> cell_v);

    int get_cnt(std::string gene_id, std::string cell_id);

    // generate a random gene counting matrix that follows gamma distribution.
    void gamma_count_matrix(double alpha, double beta);
};


// the parent class for all simulator
class FastqSimulator
{
public:
    // the annotation used to store exon information for each genes,
    // `Bar` include the cell barcode information, if Bar not provided, use random seq as cell barcode
    // `cell_cnt` is the number of cells
    // `sim_seed` is the seed used for all random process, to reproduce the data.
    // `cnt_mat` is the gene counting matrix
    GeneAnnotation Anno;
    Barcode Bar;
    CountSimulator Cnt_sim;
    int cell_cnt;
    unsigned sim_seed;
    std::mt19937 eng{std::random_device{}()};
    std::uniform_int_distribution<int> uni_int_dist{0, 3};

    FastqSimulator(std::string annofn, unsigned seed);
    FastqSimulator(std::string annofn);
    // get transcript ATCG sequence from fasta reference and gene annotation
    // skip overlap exons
    // only do on gene level. do not use transcirpt annotation
    // TODO: do on transcript level
    std::string get_transcript_seq(Gene ge, Fa_rec fa);  // kseq_t is the return value of htslib fasta reader

    // generate random barcode sequence of length `len`
    std::string gen_random_seq(int len);

    void gen_gene_expression(std::string method, std::vector<double> param_vec);

    virtual void makefq(std::string R1fn, std::string R2fn, std::string reffa) = 0;


};

// fastq reads that simulate the CEL-seq2 protocol
class Celseq2Simulator: public FastqSimulator
{
public:
    Celseq2Simulator(std::string annofn, std::string barfn);
    Celseq2Simulator(std::string annofn, std::string barfn, unsigned seed);
    void makefq(std::string R1fn, std::string R2fn, std::string reffa);
    void makefq(std::string R1fn, std::string R2fn, std::string reffa, int UMI_len, int r_len, int frag_mean, int dup_mean);
};

#endif