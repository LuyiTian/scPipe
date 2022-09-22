#include <Rcpp.h>
using namespace Rcpp;

// Given a DataFrame representing the chr, start and end of TSSes, a range to explore and a bin size
// Return a DataFrame representing the chr, start and end of all bins, ordered by TSS
// The resulting DataFrame will have n*n_bins where n is the number of TSSes

// [[Rcpp::export]]
DataFrame get_all_TSS_bins(DataFrame tss_df, int range, int bin_size) {
  NumericVector starts = tss_df["start"];
  StringVector chr = tss_df["chr"];
  int n = starts.size(), n_bins = (int)range/bin_size-1;
  
  NumericVector bin_starts(n*n_bins), bin_ends(n*n_bins);
  StringVector chrs(n*n_bins);
  
  int start;
  for (int i = 0; i < n; i++) {
    start = starts[i]-(int)range/2+(int)bin_size/2;
    for (int j = 0; j < n_bins; j++) {
      bin_starts[n_bins*i+j] = start + j*bin_size;
      bin_ends[n_bins*i+j] = bin_starts[n_bins*i+j] + bin_size-1;
      chrs[n_bins*i+j] = chr[i];
    }
  } 
  
  return DataFrame::create(Named("chr") = chrs,
                           Named("start") = bin_starts, 
                           Named("end") = bin_ends);
}