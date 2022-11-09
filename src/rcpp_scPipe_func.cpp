#include <Rcpp.h>
#include "trimbarcode.h"
#include "parsecount.h"
#include "parsebam.h"
#include "cellbarcode.h"
#include "transcriptmapping.h"
#include "detect_barcode.h"
#include "Timer.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace Rcpp;


read_s get_read_structure(Rcpp::NumericVector bs1,
                          Rcpp::NumericVector bl1,
                          Rcpp::NumericVector bs2,
                          Rcpp::NumericVector bl2,
                          Rcpp::NumericVector us,
                          Rcpp::NumericVector ul)
{
  read_s s = {};
  s.id1_st = Rcpp::as<int>(bs1);  // id1 start
  s.id1_len = Rcpp::as<int>(bl1);    // id1 length
  s.id2_st = Rcpp::as<int>(bs2);     // id2 start
  s.id2_len = Rcpp::as<int>(bl2);    // id2 length
  s.umi_st = Rcpp::as<int>(us);     // umi start
  s.umi_len = Rcpp::as<int>(ul);    // umi length
  return s;
}

filter_s get_filter_structure(Rcpp::NumericVector rmlow,
                              Rcpp::NumericVector rmN,
                              Rcpp::NumericVector minq,
                              Rcpp::NumericVector numbq)
{
  int c_rmlow = Rcpp::as<int>(rmlow);
  int c_rmN = Rcpp::as<int>(rmN);
  filter_s f = {};
  f.if_check_qual = c_rmlow==1?true:false;
  f.if_remove_N = c_rmN==1?true:false;
  
  f.min_qual = Rcpp::as<int>(minq);
  f.num_below_min = Rcpp::as<int>(numbq);
  return f;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

void rcpp_sc_trim_barcode_paired(Rcpp::CharacterVector outfq,
                                 Rcpp::CharacterVector r1,
                                 Rcpp::CharacterVector r2,
                                 Rcpp::NumericVector bs1,
                                 Rcpp::NumericVector bl1,
                                 Rcpp::NumericVector bs2,
                                 Rcpp::NumericVector bl2,
                                 Rcpp::NumericVector us,
                                 Rcpp::NumericVector ul,
                                 Rcpp::NumericVector rmlow,
                                 Rcpp::NumericVector rmN,
                                 Rcpp::NumericVector minq,
                                 Rcpp::NumericVector numbq,
                                 Rcpp::LogicalVector write_gz) {
  
  std::string c_outfq = Rcpp::as<std::string>(outfq);
  std::string c_r1 = Rcpp::as<std::string>(r1);
  std::string c_r2 = Rcpp::as<std::string>(r2);
  read_s s = get_read_structure(bs1, bl1, bs2, bl2, us, ul);
  filter_s fl = get_filter_structure(rmlow, rmN, minq, numbq);
  bool c_write_gz = Rcpp::as<bool>(write_gz);
  
  Rcpp::Rcout << "trimming fastq file..." << "\n";
  
  Timer timer;
  timer.start();
  
  paired_fastq_to_fastq((char *)c_r1.c_str(), (char *)c_r2.c_str(), (char *)c_outfq.c_str(), s, fl, c_write_gz);
  
  Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

void rcpp_sc_exon_mapping(Rcpp::CharacterVector inbam,
                          Rcpp::CharacterVector outbam,
                          Rcpp::CharacterVector annofn,
                          Rcpp::CharacterVector am,
                          Rcpp::CharacterVector ge,
                          Rcpp::CharacterVector bc,
                          Rcpp::CharacterVector mb,
                          Rcpp::NumericVector bc_len,
                          Rcpp::CharacterVector bc_vector,
                          Rcpp::NumericVector UMI_len,
                          Rcpp::NumericVector stnd,
                          Rcpp::NumericVector fix_chr,
                          Rcpp::NumericVector nthreads)
{
  //std::string c_inbam = Rcpp::as<std::string>(inbam);
  std::string c_outbam = Rcpp::as<std::string>(outbam);
  
  std::string c_am = Rcpp::as<std::string>(am);
  std::string c_ge = Rcpp::as<std::string>(ge);
  std::string c_bc = Rcpp::as<std::string>(bc);
  std::string c_mb = Rcpp::as<std::string>(mb);
  
  int c_bc_len = Rcpp::as<int>(bc_len);
  int c_UMI_len = Rcpp::as<int>(UMI_len);
  bool c_stnd = Rcpp::as<int>(stnd)==1?true:false;
  bool c_fix_chr = Rcpp::as<int>(fix_chr)==1?true:false;
  std::vector<std::string> c_inbam_vec = Rcpp::as<std::vector<std::string>>(inbam);
  std::vector<std::string> c_bc_vec = Rcpp::as<std::vector<std::string>>(bc_vector);
  std::vector<std::string> token = Rcpp::as<std::vector<std::string>>(annofn);
  int c_nthreads = Rcpp::as<int>(nthreads);
  
  Mapping a = Mapping();
  Rcpp::Rcout << "adding annotation files..." << "\n";
  
  Timer timer;
  for (const auto& n : token)
  {
    timer.start();
    a.add_annotation(n, c_fix_chr);
    Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
  }
  
  Rcpp::Rcout << "annotating exon features..." << "\n";
  
  timer.start();
  a.parse_align_warpper(c_inbam_vec, c_bc_vec, c_outbam, c_stnd, c_am, c_ge, c_bc, c_mb, c_bc_len, c_UMI_len, c_nthreads);
  //a.parse_align(c_inbam, c_outbam, c_stnd, c_am, c_ge, c_bc, c_mb, c_bc_len, c_UMI_len);
  
  Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
}

// [[Rcpp::export]]
void rcpp_sc_exon_mapping_df_anno(
    Rcpp::CharacterVector inbam,
    Rcpp::CharacterVector outbam,
    Rcpp::DataFrame anno,
    Rcpp::CharacterVector am,
    Rcpp::CharacterVector ge,
    Rcpp::CharacterVector bc,
    Rcpp::CharacterVector mb,
    Rcpp::NumericVector bc_len,
    Rcpp::CharacterVector bc_vector,
    Rcpp::NumericVector UMI_len,
    Rcpp::NumericVector stnd,
    Rcpp::NumericVector fix_chr,
    Rcpp::NumericVector nthreads)
{
  //std::string c_inbam = Rcpp::as<std::string>(inbam);
  std::string c_outbam = Rcpp::as<std::string>(outbam);
  
  std::string c_am = Rcpp::as<std::string>(am);
  std::string c_ge = Rcpp::as<std::string>(ge);
  std::string c_bc = Rcpp::as<std::string>(bc);
  std::string c_mb = Rcpp::as<std::string>(mb);
  
  int c_bc_len = Rcpp::as<int>(bc_len);
  int c_UMI_len = Rcpp::as<int>(UMI_len);
  bool c_stnd = Rcpp::as<int>(stnd)==1?true:false;
  bool c_fix_chr = Rcpp::as<int>(fix_chr)==1?true:false;
  std::vector<std::string> c_inbam_vec = Rcpp::as<std::vector<std::string>>(inbam);
  std::vector<std::string> c_bc_vec = Rcpp::as<std::vector<std::string>>(bc_vector);
  int c_nthreads = Rcpp::as<int>(nthreads);
  
  Mapping a = Mapping();
  Rcpp::Rcout << "adding annotation files..." << "\n";
  
  Timer timer;
  timer.start();
  a.add_annotation(anno, c_fix_chr);
  Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
  
  Rcpp::Rcout << "annotating exon features..." << "\n";
  
  a.parse_align_warpper(c_inbam_vec, c_bc_vec, c_outbam, c_stnd, c_am, c_ge, c_bc, c_mb, c_bc_len, c_UMI_len, c_nthreads);
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

void rcpp_sc_demultiplex(Rcpp::CharacterVector inbam,
                         Rcpp::CharacterVector outdir,
                         Rcpp::CharacterVector bc_anno,
                         Rcpp::NumericVector max_mis,
                         Rcpp::CharacterVector am,
                         Rcpp::CharacterVector ge,
                         Rcpp::CharacterVector bc,
                         Rcpp::CharacterVector mb,
                         Rcpp::CharacterVector mito,
                         Rcpp::LogicalVector has_UMI,
                         Rcpp::NumericVector nthreads)
{
  std::string c_inbam = Rcpp::as<std::string>(inbam);
  std::string c_outdir = Rcpp::as<std::string>(outdir);
  std::string c_bc_anno = Rcpp::as<std::string>(bc_anno);
  std::string c_mito = Rcpp::as<std::string>(mito);
  
  std::string c_am = Rcpp::as<std::string>(am);
  std::string c_ge = Rcpp::as<std::string>(ge);
  std::string c_bc = Rcpp::as<std::string>(bc);
  std::string c_mb = Rcpp::as<std::string>(mb);
  
  int c_max_mis = Rcpp::as<int>(max_mis);
  bool c_has_UMI = Rcpp::as<bool>(has_UMI);
  int c_nthreads = Rcpp::as<int>(nthreads);
  
  Barcode bar;
  bar.read_anno(c_bc_anno);
  
  Rcpp::Rcout << "demultiplexing reads by barcode..." << "\n";
  
  Timer timer;
  timer.start();
  
  Bamdemultiplex bam_de = Bamdemultiplex(c_outdir, bar, c_bc, c_mb, c_ge, c_am, c_mito);
  
  bam_de.barcode_demultiplex(c_inbam, c_max_mis, c_has_UMI, c_nthreads);
  bam_de.write_statistics("overall_stat", "chr_stat", "cell_stat");
  Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

void rcpp_sc_clean_bam(Rcpp::CharacterVector inbam,
                       Rcpp::CharacterVector outbam,
                       Rcpp::CharacterVector bc_anno,
                       Rcpp::NumericVector max_mis,
                       Rcpp::CharacterVector am,
                       Rcpp::CharacterVector ge,
                       Rcpp::CharacterVector bc,
                       Rcpp::CharacterVector mb,
                       Rcpp::CharacterVector mito,
                       Rcpp::NumericVector nthreads)
{
  std::string c_inbam = Rcpp::as<std::string>(inbam);
  std::string c_outbam = Rcpp::as<std::string>(outbam);
  std::string c_bc_anno = Rcpp::as<std::string>(bc_anno);
  std::string c_mito = Rcpp::as<std::string>(mito);
  
  std::string c_am = Rcpp::as<std::string>(am);
  std::string c_ge = Rcpp::as<std::string>(ge);
  std::string c_bc = Rcpp::as<std::string>(bc);
  std::string c_mb = Rcpp::as<std::string>(mb);
  
  int c_max_mis = Rcpp::as<int>(max_mis);
  int c_nthreads = Rcpp::as<int>(nthreads);
  
  Barcode bar;
  bar.read_anno(c_bc_anno);
  
  Rcpp::Rcout << "update bam file with corrected barcode..." << "\n";
  
  Timer timer;
  timer.start();
  
  Bamdemultiplex bam_de = Bamdemultiplex("", bar, c_bc, c_mb, c_ge, c_am, c_mito);
  bam_de.clean_bam_barcode(c_inbam, c_outbam, c_max_mis, c_nthreads);
  Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

void rcpp_sc_gene_counting(Rcpp::CharacterVector outdir,
                           Rcpp::CharacterVector bc_anno,
                           Rcpp::NumericVector UMI_cor,
                           Rcpp::NumericVector gene_fl)
{
  std::string c_outdir = Rcpp::as<std::string>(outdir);
  std::string c_bc_anno = Rcpp::as<std::string>(bc_anno);
  int c_UMI_cor = Rcpp::as<int>(UMI_cor);
  bool c_gene_fl = Rcpp::as<int>(gene_fl)==1?true:false;
  
  Barcode bar;
  
  Rcpp::Rcout << "summarising gene counts..." << "\n";
  
  Timer timer;
  timer.start();
  
  bar.read_anno(c_bc_anno);
  get_counting_matrix(bar, c_outdir, c_UMI_cor, c_gene_fl);
  
  Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

void rcpp_sc_detect_bc(Rcpp::CharacterVector infq,
                       Rcpp::CharacterVector outcsv,
                       Rcpp::CharacterVector prefix,
                       Rcpp::NumericVector bc_len,
                       Rcpp::NumericVector max_reads,
                       Rcpp::NumericVector number_of_cells,
                       Rcpp::NumericVector min_count,
                       Rcpp::NumericVector max_mismatch,
                       Rcpp::CharacterVector white_list)
{
  std::string c_infq = Rcpp::as<std::string>(infq);
  std::string c_outcsv = Rcpp::as<std::string>(outcsv);
  std::string c_prefix = Rcpp::as<std::string>(prefix);
  std::string c_white_list = Rcpp::as<std::string>(white_list);
  int c_bc_len = Rcpp::as<int>(bc_len);
  int c_max_reads = Rcpp::as<int>(max_reads);
  int c_min_count = Rcpp::as<int>(min_count);
  int c_max_mismatch = Rcpp::as<int>(max_mismatch);
  int c_number_of_cells = Rcpp::as<int>(number_of_cells);
  
  std::unordered_map<std::string, int> counter = summarize_barcode(c_infq, c_bc_len, c_max_reads, c_max_mismatch, c_min_count, c_white_list);
  write_barcode_summary(c_outcsv, c_prefix, counter, c_number_of_cells);
}

















// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

std::vector<int> rcpp_sc_atac_trim_barcode(
    Rcpp::CharacterVector outfq,
    Rcpp::CharacterVector r1,
    Rcpp::CharacterVector r3,
    Rcpp::StringVector barcode_file,
	Rcpp::CharacterVector valid_barcode_file,
    Rcpp::NumericVector umi_start,
    Rcpp::NumericVector umi_len,
    Rcpp::CharacterVector umi_in,
    Rcpp::LogicalVector write_gz,
    Rcpp::LogicalVector rmN,
    Rcpp::LogicalVector rmlow,
    Rcpp::IntegerVector min_qual,
    Rcpp::IntegerVector num_below_min,
    Rcpp::IntegerVector id1_st,
    Rcpp::IntegerVector id1_len,
    Rcpp::IntegerVector id2_st,
    Rcpp::IntegerVector id2_len) {
  
  std::string c_outfq = Rcpp::as<std::string>(outfq);
  std::string c_r1 = Rcpp::as<std::string>(r1);
  std::string c_r3 = Rcpp::as<std::string>(r3);
  std::string c_barcode_file = Rcpp::as<std::string>(barcode_file(0));
  std::string c_valid_barcode_file = Rcpp::as<std::string>(valid_barcode_file);
  //std::string c_barcode_file = Rcpp::as<std::string>(barcode_file);
  std::string c_umi_in = Rcpp::as<std::string>(umi_in);
  
  
  //int c_start = Rcpp::as<int>(start);
  //int c_len = Rcpp::as<int>(len);
  int c_umi_start = Rcpp::as<int>(umi_start);
  int c_umi_len = Rcpp::as<int>(umi_len);
  bool c_write_gz = Rcpp::as<bool>(write_gz);
  bool c_rmN = Rcpp::as<bool>(rmN);
  bool c_rmlow = Rcpp::as<bool>(rmlow);
  int c_min_qual = Rcpp::as<int>(min_qual);
  int c_num_below_min = Rcpp::as<int>(num_below_min);
  int c_id1_st = Rcpp::as<int>(id1_st);
  int c_id1_len = Rcpp::as<int>(id1_len);
  int c_id2_st = Rcpp::as<int>(id2_st);
  int c_id2_len = Rcpp::as<int>(id2_len);
  
  Timer timer;
  timer.start();
  
  std::vector<int> out_vec = sc_atac_paired_fastq_to_csv(
    (char *)c_r1.c_str(), 
    (char *)c_r3.c_str(), 
    (char *)c_outfq.c_str(), 
    (char *)c_barcode_file.c_str(), 
	(char *)c_valid_barcode_file.c_str(),
    c_umi_start, 
    c_umi_len, 
    (char*)c_umi_in.c_str(), 
    c_write_gz, 
    c_rmN,
    c_rmlow,
    c_min_qual,
    c_num_below_min,
    c_id1_st,
    c_id1_len,
    c_id2_st,
    c_id2_len
  );
  
  Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
  
  return(out_vec);
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
std::vector<int> rcpp_sc_atac_trim_barcode_paired(
		Rcpp::CharacterVector outfq,
		Rcpp::CharacterVector r1,
		Rcpp::StringVector r2_list,
		Rcpp::CharacterVector r3,
		Rcpp::CharacterVector valid_barcode_file,
		Rcpp::LogicalVector write_gz,
		Rcpp::LogicalVector rmN,
		Rcpp::LogicalVector rmlow,
		Rcpp::IntegerVector min_qual,
		Rcpp::IntegerVector num_below_min,
		Rcpp::LogicalVector no_reverse_complement) {
  
	std::string c_outfq = Rcpp::as<std::string>(outfq);
	std::string c_r1 = Rcpp::as<std::string>(r1);
	std::vector<std::string> c_r2_list; 
	
	for( int i=0; i < r2_list.size(); i++ ){
		c_r2_list.push_back(Rcpp::as<std::string>(r2_list(i)));
	}
	
	std::string c_r3 = Rcpp::as<std::string>(r3);
	std::string c_valid_barcode = Rcpp::as<std::string>(valid_barcode_file);
	bool c_write_gz = Rcpp::as<bool>(write_gz);
	bool c_rmN = Rcpp::as<bool>(rmN);
	bool c_nRC = Rcpp::as<bool>(no_reverse_complement);
	
	bool c_rmlow = Rcpp::as<bool>(rmlow);
	int c_min_qual = Rcpp::as<int>(min_qual);
	int c_num_below_min = Rcpp::as<int>(num_below_min);
	
	Timer timer;
	timer.start();
	
	std::vector<int> out_vec = sc_atac_paired_fastq_to_fastq(
		(char *)c_r1.c_str(), 
		c_r2_list, 
		(char *)c_r3.c_str(), 
		(char *)c_valid_barcode.c_str(),
		(char *)c_outfq.c_str(), 
		c_write_gz, 
		c_rmN,
		c_rmlow,
		c_min_qual,
		c_num_below_min,
		c_nRC);
	
	Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
	
	return(out_vec);
}


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
void rcpp_sc_atac_bam_tagging(Rcpp::CharacterVector inbam,
                              Rcpp::CharacterVector outbam,
                              Rcpp::CharacterVector bc,
                              Rcpp::CharacterVector mb,
                              Rcpp::NumericVector nthreads)
{
  std::string c_outbam = Rcpp::as<std::string>(outbam);
  
  std::string c_bc = Rcpp::as<std::string>(bc);
  std::string c_mb = Rcpp::as<std::string>(mb);
  std::vector<std::string> c_inbam_vec = Rcpp::as<std::vector<std::string>>(inbam);
  int c_nthreads = Rcpp::as<int>(nthreads);
  
  Mapping a = Mapping();
  Timer timer;
  timer.start();
  a.sc_atac_parse_align_warpper(c_inbam_vec, c_outbam, c_bc, c_mb, c_nthreads);
  Rcpp::Rcout << "time elapsed: " << timer.time_elapsed() << "\n\n";
}






// [[Rcpp::export]]
void rcpp_fasta_bin_bed_file(std::string in_filename, std::string out_filename, int bin_size){
  // Reads a fasta file and writes a bed file with the according bin_size.
  // Function called in sc_atac_feature_counting()
  
  int read_lines_counter = 0;
  int read_sequences_counter = 0;
  int written_lines_counter = 0;
  int seq_length = 0;
  
  std::ifstream input(in_filename.c_str());
  std::string line, name, sequence;
  
  std::ofstream output;
  output.open(out_filename);
  
  double n_bins;
  
  while (std::getline(input, line)) {
    read_lines_counter++;
    
    // Ignore blank lines
    if(line.empty())
      continue;
    
    if (line[0] == '>') {
      // output previous line before overwriting name
      
      if(!name.empty()){ // Only do this if name actually contains something.
        read_sequences_counter++;
        
        seq_length = sequence.length();
        
        // divide the sequence in bins and write in file
        n_bins = floor((double)seq_length / (double)bin_size) + 1;
        int bin_counter = 0;
        while(bin_counter < n_bins - 1){
          output << name << "\t" << bin_counter*bin_size + 1 << "\t" << (bin_counter + 1)*bin_size << std::endl;
          written_lines_counter++;
          bin_counter++;
        }
        output << name << "\t" << bin_counter*bin_size + 1 << "\t" << seq_length << std::endl;
        written_lines_counter++;
        
      }
      
      //name = line.substr(1); old to remove > character, but keeps in any additional header information
      name = line.substr(1, std::min(line.find_first_of(' '), line.find_first_of("\t")));
      sequence.clear();
    } else { //  if (line[0] != '>')
      sequence += line;
    }
  } // end while
  
  // output final entry but only if name actually contains something
  if(!name.empty()){
    seq_length = sequence.length();
    
    read_sequences_counter++;
    
    // divide the sequence in bins and print
    n_bins = floor((double)seq_length / (double)bin_size) + 1;
    int bin_counter = 0;
    while(bin_counter < n_bins - 1){
      output << name << "\t" << bin_counter*bin_size + 1 << "\t" << (bin_counter + 1)*bin_size << std::endl;
      written_lines_counter++;
      bin_counter++;
    }
    output << name << "\t" << bin_counter*bin_size + 1 << "\t" << seq_length << std::endl;
    written_lines_counter++;
    
  }
  
  Rcout << "Number of read lines: " << read_lines_counter << std::endl;
  Rcout << "Number of read sequences: " << read_sequences_counter << std::endl;
  Rcout << "Number of written lines: " << written_lines_counter << std::endl;
  output.close();
  
}






// [[Rcpp::export]]
void rcpp_append_chr_to_bed_file(std::string in_filename, std::string out_filename){
  // Appends the string "chr" to each line of in_filename and saves it in out_filename
  // Function called in sc_atac_feature_counting()
  
  int read_lines_counter = 0;
  
  std::ifstream input(in_filename.c_str());
  std::string line, name, sequence;
  
  std::ofstream output;
  output.open(out_filename);
  
  
  while (std::getline(input, line)) {
    read_lines_counter++;
    
    output << "chr" << line << std::endl;
    
  } // end while
  
  
  Rcout << "Number of read and modified lines: " << read_lines_counter << std::endl;
  output.close();
  
}











