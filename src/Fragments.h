#ifndef BAM_UTIL_FUNCS_H
#define BAM_UTIL_FUNCS_H

// seqment->data is qname-cigar-seq-qual-aux
// qname is at bam1_qname(seqment) (which is (char*)seqment->data)
			
/*
Get the mapping quality of the alignment
@param b bam1_t pointer to an alignment
@return uint8_t mapping quality of alignment
*/
#define bam_mapping_qual(b) (b->core.qual)

/*
Get the 0-based leftmost coordinate of an alignment
@param b bam1_t pointer to an alignment
@return int32_t pos
*/
#define bam_alignment_start(b) (b->core.pos)

/*
Convert the uint8_t aux value to a string
Shorter version of htslib/sam.h bam_aux2Z
Removed return of 0 (null terminator) for bad input,
as Rcpp could not handle it and crashed
Now, valid input can be checked by the returning string's length
@param s uint8_t * pointer to tag data returned by bam_aux_get()
@return std::string of tag data, or empty string if tag type is not Z
*/
#define bam_aux2string(s) ( (s != NULL) ? std::string((*s++ == 'Z') ? (char *)s : "") : "")

#endif


#ifndef FRAGMENTS_HEADER
#define FRAGMENTS_HEADER
#include <string>
#include <map>

struct FragmentStruct {
	std::string chromosome;			// 0
	int32_t start;					// 1
	int32_t end;					// 2
	std::string cell_barcode;		// 3
	bool complete;					// 4
	int sum;
};

typedef std::map<std::string, FragmentStruct> FragmentMap;

#endif