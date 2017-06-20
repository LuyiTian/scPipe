// fq simulator
#include "fq_simulator.h"

using namespace Rcpp;

FaReader::FaReader(std::string fafn)
{
    fafile.open(fafn);
}

bool FaReader::readone()
{
    if (!line.empty()){
        fa.name = line.substr(1);
    }
    if (!fa.seq.empty())
    {
        fa.seq.clear();
    }
    while (std::getline(fafile, line))
    {
        if(line.empty())
            continue;

        if (line[0] == '>') 
        {
            if(!fa.name.empty())
            {
                return true;
            }
            else
            {
                fa.name = line.substr(1);
            }
        }
        else
        {
            fa.seq += line;
        }
    }
    if (!fa.seq.empty())
    {
        return true;
    }
    else
    {
        return false;  // finished
    }

}


FqWriter::FqWriter(std::string fqfn)
{
    fqfile.open(fqfn);
}


void FqWriter::writeone()
{
    fqfile << "@" << fq.name << std::endl << fq.seq << std::endl << "+" << std::endl << fq.qual << std::endl;
}


void CountSimulator::init_mat(std::vector<std::string> gene_v, std::vector<std::string> cell_v){
    cnt_mat.resize(gene_v.size(), cell_v.size());  // assume the matrix is not too big
    cnt_mat.fill(0);

    ix_gene = gene_v;
    for (int ix = 0; ix < gene_v.size(); ix++)
    {
        gene_ix[gene_v[ix]] = ix;
    }

    ix_cell = cell_v;
    for (int ix = 0; ix < cell_v.size(); ix++)
    {
        cell_ix[cell_v[ix]] = ix;
    }
}

CountSimulator::CountSimulator(): CountSimulator::CountSimulator(std::chrono::high_resolution_clock::now().time_since_epoch().count()) {}

CountSimulator::CountSimulator(unsigned seed): sim_seed(seed){}


void CountSimulator::gamma_count_matrix(double alpha, double beta)
{
    auto rand_int = std::bind(std::gamma_distribution<double>(alpha, beta),
                           std::mt19937(sim_seed));

    for (int ri=0; ri< cnt_mat.rows(); ri++)
    {
        for (int ci=0; ci< cnt_mat.cols(); ci++)
        {
            cnt_mat(ri, ci) = std::floor(rand_int());
        }
    }
}


int CountSimulator::get_cnt(std::string gene_id, std::string cell_id)
{
    return cnt_mat(gene_ix[gene_id], cell_ix[cell_id]);
}


FastqSimulator::FastqSimulator(std::string annofn, unsigned seed): eng(seed)
{
    if(annofn.substr(annofn.find_last_of(".") + 1) == "gff3") 
    {
        Anno.parse_gff3_annotation(annofn, false);
    } 
    else 
    {
        Anno.parse_bed_annotation(annofn, false);
    }
    Cnt_sim.sim_seed = seed;
}


FastqSimulator::FastqSimulator(std::string annofn): FastqSimulator(annofn, std::chrono::high_resolution_clock::now().time_since_epoch().count()) {}

std::string FastqSimulator::get_transcript_seq(Gene ge, Fa_rec fa)
{
    std::string transcript = "";
    Interval tmp_cache(-1, -1, 0);
    for (auto i = ge.exon_vec.begin(); i != ge.exon_vec.end(); ++i)
    {
        if ((*i) > tmp_cache)
        {
            tmp_cache.st = i->st;
            tmp_cache.en = i->en;
            if (i->en > fa.seq.size())
            {
                std::stringstream err_msg;
                err_msg << "ERROR: the right coordinate of exon ("<< i->en <<") is larger than sequence length ("<< fa.seq.size() <<"). check your annotation.";
                Rcpp::stop(err_msg.str());
            }
            transcript += fa.seq.substr((i->st)-1, (i->en)-(i->st)+1);  // TODO the gff/bed file use one based inclusive coordinate.
        }
    }
    return transcript;
}


std::string FastqSimulator::gen_random_seq(int len)
{
    //auto dice_rand = std::bind(std::uniform_int_distribution<int>(0, BP.size()-1),
    //                       std::mt19937(sim_seed));
    std::string ran_seq = "";
    for (int i = 0; i < len; ++i)
    {
        ran_seq += BP[uni_int_dist(eng)];
    }

    return ran_seq;

}

void FastqSimulator::gen_gene_expression(std::string method, std::vector<double> param_vec)
{
    if (method == "gamma_random")
    {
        Cnt_sim.gamma_count_matrix(param_vec[0], param_vec[1]);
    }
}


Celseq2Simulator::Celseq2Simulator(std::string annofn, std::string barfn): Celseq2Simulator(annofn, barfn, std::chrono::high_resolution_clock::now().time_since_epoch().count()) {}

Celseq2Simulator::Celseq2Simulator(std::string annofn, std::string barfn, unsigned seed): FastqSimulator(annofn, seed)
{
    Bar.read_anno(barfn);
    cell_cnt = Bar.cellid_list.size();
    Cnt_sim.init_mat(Anno.get_genelist(), Bar.cellid_list);
}


void Celseq2Simulator::makefq(std::string R1fn, std::string R2fn, std::string reffa, int UMI_len, int r_len, int frag_mean, int dup_mean)
{
    int dup_count, molecular_count, frag_len;
    auto get_dup = std::bind(std::gamma_distribution<double>(1, dup_mean),
                           eng);  // assume the PCR duplicate number follows gamma distribution
    auto get_frag = std::bind(std::normal_distribution<double>(frag_mean, frag_mean/5),
                           eng);  // assume the fragment length follows normal distribution
     // used for getting the random UMI sequence from an ATGC vector

    std::string UMI_seq, transcript_seq;

    int gene_digit = std::to_string(Anno.get_genelist().size()).length()+1;
    int cell_digit = std::to_string(cell_cnt).length()+1;
    int molecule_digit = std::to_string(Cnt_sim.cnt_mat.maxCoeff()).length()+1;
    FqWriter R1(R1fn);
    FqWriter R2(R2fn);

    int global_count = 0;
    FaReader fareader = FaReader(reffa);
    while (fareader.readone())
    {
        if(__DEBUG)
        {
            std::cout << "read fasta sequence " << fareader.fa.name << " with lenth " << fareader.fa.seq.length() << std::endl;
        }
        if (Anno.gene_dict.end() == Anno.gene_dict.find(fareader.fa.name))
        {
            std::stringstream err_msg;
            err_msg << "cannot find chromosome (" << fareader.fa.name << ") in exon annotation." << "\n";
            Rcpp::stop(err_msg.str());
        }
        for (int gene_ix=0; gene_ix<Anno.gene_dict[fareader.fa.name].size(); ++gene_ix)  // for each genes in that chromosome
        {
            if(__DEBUG){std::cout<< "simulate gene " << gene_ix << std::endl;}
            transcript_seq = get_transcript_seq(Anno.gene_dict[fareader.fa.name][gene_ix], fareader.fa);  // gene transcript

            for (int cell_ix=0; cell_ix<cell_cnt; ++cell_ix)  // for each cells
            {
                if(__DEBUG){std::cout<< "\tsimulate cell " << cell_ix << std::endl;}
                molecular_count = Cnt_sim.get_cnt(Anno.gene_dict[fareader.fa.name][gene_ix].gene_id, Bar.cellid_list[cell_ix]);

                for (int m = 0; m < molecular_count; ++m)  // for each mRNA molecule
                {
                    UMI_seq = gen_random_seq(UMI_len);  // generate UMI sequence for this mRNA
                    dup_count = std::ceil(get_dup());
                    dup_count = (dup_count<2)?2:dup_count;  // should be larger than one

                    frag_len = std::ceil(get_frag());
                    frag_len = (frag_len<(r_len+1))?(r_len+1):frag_len;  // should be larger than read length
                    frag_len = (frag_len<(transcript_seq.length()))?frag_len:transcript_seq.length();  // should be smaller than transcript length

                    for (int d = 0; d < dup_count; ++d)  // for each PCR duplicates
                    {
                        // set name
                        // fastq name : SIMULATE_SEQ::cell_idx::gene_idx::molecule_idx::duplicate_idx
                        R1.fq.name = "SIMULATE_SEQ::" + padding(cell_ix, cell_digit) + "::" + padding(gene_ix, gene_digit) + \
                            "::" + padding(m, molecule_digit) + "::" + std::to_string(d);
                        R2.fq.name = R1.fq.name;
                        // set sequence
                        R1.fq.seq = transcript_seq.substr(transcript_seq.length()-frag_len, r_len);  // no ligation. read one is transcipt
                        R2.fq.seq = UMI_seq + Bar.barcode_list[cell_ix]+std::string(10, 'N');
                        // set quality
                        R1.fq.qual = std::string(R1.fq.seq.length(), 'A');  // assume all seq has high quality (quality score=32)
                        R2.fq.qual = std::string(R2.fq.seq.length(), 'A');

                        R1.writeone();
                        R2.writeone();
                        global_count++;
                    }
                }
                if(__DEBUG){std::cout<< "\t\tsimulate molecule. count: " << molecular_count << std::endl;}
            }
        }
    }
    std::cerr << "simulate finished" << std::endl;
    std::cerr << "\tcell number: " << cell_cnt << std::endl;
    std::cerr << "\tgene number: " << Anno.get_genelist().size() << std::endl;
    std::cerr << "\ttotal number of reads: " << global_count << std::endl;
}

void Celseq2Simulator::makefq(std::string R1fn, std::string R2fn, std::string reffa)
{
    makefq(R1fn, R2fn, reffa, 6, 75, 300, 10);
}




