#include <Rcpp.h>
using namespace Rcpp;


// for each read, return 1 if it overlap with any region, otherwise 0;
// should sort both reads and regions
// [[Rcpp::export]]
IntegerVector sc_atac_getOverlaps_read2AnyRegion(DataFrame reads, DataFrame regions) {
  NumericVector start1 = reads["start"];
  NumericVector end1 = reads["end"];
  NumericVector start2 = regions["start"];
  NumericVector end2 = regions["end"];
  
  int n1 = start1.size(), n2 = start2.size();
  NumericVector midP1(n1), len1(n1), len2(n2), midP2(n2);
  IntegerVector over1(n1);
  
  len1 = (end1 - start1 + 1)/2;
  midP1 = (end1 + start1)/2;
  
  len2 = (end2 - start2 + 1)/2;
  midP2 = (end2 + start2)/2;
  int k = 0;
  for(int i=0; i<n1; i++){
    over1[i] = 0;
    for(int j=k; j<n2; j++){
      if((fabs(midP1[i] - midP2[j]) <= (len1[i]+len2[j]))){
        over1[i] = 1;
        k = j;
        break;
      }
    }
  }
  
  return(over1);
}

NumericVector getOverlaps_read2AllRegion(DataFrame reads, DataFrame regions) {
  NumericVector start1 = reads["start"];
  NumericVector end1 = reads["end"];
  
  NumericVector start2 = regions["start"];
  NumericVector end2 = regions["end"];
  
  int n1 = start1.size(), n2 = start2.size();
  NumericVector midP1(n1), len1(n1), len2(n2), midP2(n2);
  NumericVector over1(n1);
  
  len1 = (end1 - start1 + 1)/2;
  midP1 = (end1 + start1)/2;
  
  len2 = (end2 - start2 + 1)/2;
  midP2 = (end2 + start2)/2;
  int j = 0;
  int k = 0;
  
  if(start1[0] > end2[n2-1] || end1[n1-1] < start2[0]) {
    return(over1);
  }
  
  for(int i=0; i<n1; i++){
    while (k<n2 && fabs(midP1[i] - midP2[k]) > (len1[i]+len2[k])){
      k++;
    }
    // the kth element is the first element overlapped with the current fragment
    // take advantage of sorted fragments and regions, for next fragment,
    // it will start searching from kth region
    j=k;
    
    // current frag not overlapped any region then search from beginning
    if(j >= n2 -1) k = 0, j = 0;  
    while (j<n2 && fabs(midP1[i] - midP2[j]) <= (len1[i]+len2[j])){
      over1[i] = over1[i] + 1;
      j++;
    }
  }
  
  return(over1);
}

// [[Rcpp::export]]
NumericVector sc_atac_getOverlaps_tss2Reads(DataFrame regions, DataFrame left_flank, DataFrame reads) {
  NumericVector start1 = regions["start"];
  
  int n = start1.size();
  NumericVector over(n), left_over(n);
  NumericVector tss_enrich_score(n);
  over = getOverlaps_read2AllRegion(regions, reads);
  left_over = getOverlaps_read2AllRegion(left_flank, reads);
  tss_enrich_score = (over+1.0) / (left_over+1.0) ;
  return(tss_enrich_score);
}    