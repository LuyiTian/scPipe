#include <Rcpp.h>
#include "fq_simulator.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

void rcpp_generate_celseq2_data(Rcpp::CharacterVector r1fn,
                                Rcpp::CharacterVector r2fn,
                                Rcpp::CharacterVector annofn,
                                Rcpp::CharacterVector bc_anno,
                                Rcpp::CharacterVector fafn,
                                Rcpp::NumericVector UMI_len,
                                Rcpp::NumericVector r_len,
                                Rcpp::NumericVector frag_mean,
                                Rcpp::NumericVector dup_mean,
                                Rcpp::CharacterVector ran_dist,
                                Rcpp::NumericVector param,
                                Rcpp::NumericVector seed)
{
  std::string c_r1fn = Rcpp::as<std::string>(r1fn);
  std::string c_r2fn = Rcpp::as<std::string>(r2fn);
  std::string c_annofn = Rcpp::as<std::string>(annofn);
  std::string c_bc_anno = Rcpp::as<std::string>(bc_anno);
  std::string c_fafn = Rcpp::as<std::string>(fafn);
  std::string c_ran_dist = Rcpp::as<std::string>(ran_dist);
  std::vector<double> c_param = Rcpp::as<std::vector<double>>(param);

  int c_UMI_len = Rcpp::as<int>(UMI_len);
  int c_r_len = Rcpp::as<int>(r_len);
  int c_frag_mean = Rcpp::as<int>(frag_mean);
  int c_dup_mean = Rcpp::as<int>(dup_mean);

  unsigned c_seed = Rcpp::as<unsigned>(seed);
  if (c_seed==0)
  {
    Celseq2Simulator celseq2_sim(c_annofn, c_bc_anno);
    celseq2_sim.gen_gene_expression(c_ran_dist, c_param);
    celseq2_sim.makefq(c_r1fn, c_r2fn, c_fafn, c_UMI_len, c_r_len, c_frag_mean, c_dup_mean);
  }
  else{
    Celseq2Simulator celseq2_sim(c_annofn, c_bc_anno, c_seed);
    celseq2_sim.gen_gene_expression(c_ran_dist, c_param);
    celseq2_sim.makefq(c_r1fn, c_r2fn, c_fafn, c_UMI_len, c_r_len, c_frag_mean, c_dup_mean);
  }
}

