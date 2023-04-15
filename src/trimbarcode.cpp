#include "trimbarcode.h"

#include <string>
#include <set>
#include <zlib.h> // for reading compressed .fq file
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <Rcpp.h>

#include "config_hts.h"
#include "utils.h"
#include "Trie.h"

using namespace Rcpp;

#ifndef INIT_KSEQ
#define INIT_KSEQ
KSEQ_INIT(gzFile, gzread)
#endif


std::vector<gzFile> open_gz_files(std::vector<std::string> files) {
	std::vector<gzFile> fq2_list;
    for(int i = 0; i < (int)files.size(); i++){
        char* fq2_fn = (char *)files[i].c_str();
        gzFile fq2 = gzopen(fq2_fn, "r");
        if (!fq2) {
            file_error(fq2_fn);
        }
        fq2_list.push_back(fq2);
    }

	return fq2_list;
}

// read the valid barcode file
// expects a text file or a text file with gz compression
// uses kseq stream processing to manage the gz file
std::vector<std::string> readBarcodes(const char *barcode_fn) {
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};

	fp = gzopen(barcode_fn, "r");
	ks = ks_init(fp);
	std::vector<std::string> out;

	while (ks_getuntil(ks, '\n', &str, 0) >= 0) {
		std::string line(str.s);
		size_t comma1 = line.find(',');
		size_t comma2 = line.find(',', comma1 + 1);
		std::string secondColumn = line.substr(
			comma1 != std::string::npos ? comma1 + 1 : 0, 
			(comma2 != std::string::npos ? comma2 - 1 : line.size()) - comma1);
		out.push_back(secondColumn);
	}

	ks_destroy(ks);
	gzclose(fp);
	free(str.s);
	return out;
}

// handle the building of a Trie using a barcode file
void preprocessBarcodes(const std::vector<std::string> &barcodes, Trie &trie) {
	// add barcodes to the trie, using the barcode index as both original and new seq id (because these barcodes aren't sorted)
	for (int i = 0; i < (int)barcodes.size(); i++) {
		trie.Add_String(barcodes[i], i, i);
	}
}

int count_unmatched_barcodes(const std::vector<std::string> &barcodes, const Trie &valid_barcodes, int total) {
	std::vector<gzFile> fq2_list = open_gz_files(barcodes);
	std::vector<kseq_t*> seq2_list;
	for (int i = 0; i < (int)fq2_list.size(); i++) {
		seq2_list.push_back(kseq_init(fq2_list[i]));
	}

	int failed = 0;
	for (int i = 0; i < total; i++) {
		bool noMatch = false;
		for (int j = 0; j < (int)seq2_list.size(); j++) {
			kseq_t* seq2 = seq2_list[j];
			if (kseq_read(seq2) >= 0) {
				std::string barcode (seq2->seq.s);
				std::vector<MismatchResult> matches = valid_barcodes.Locate_Seq_Mismatches(barcode, 0, barcode.size());
				if (matches.size() == 0) {
					noMatch = true;
				}
			}
		}
		if (noMatch) {
			failed++;
		}
	}

	for (int j = 0; j < (int)seq2_list.size(); j++) {
		kseq_t* seq2 = seq2_list[j];
        kseq_destroy(seq2);// free seq
        gzclose(fq2_list[j]);
	}

	return failed;
}

char charComplement(char c) {
  char out = 'N';
  if (c == 'A') out = 'T';
  if (c == 'T') out = 'A';
  if (c == 'C') out = 'G';
  if (c == 'G') out = 'C';
  return out;
}

std::string reverseComplement(const char *seq, size_t l) {
	std::string res;
	res.reserve(l);
	for (int i = l - 1; i >= 0; i--) {
		res.append(1, charComplement(seq[i]));
	}
	return res;
}

std::vector<std::string> processAndBuildTrie(const char *valid_barcode_fn, const std::vector<std::string> &fq2_fn_list, Trie &validBarcodeTrie, bool noReverseComplement, bool &useReverseComplement) {
	std::vector<std::string> validBarcodes = readBarcodes(valid_barcode_fn);
	preprocessBarcodes(validBarcodes, validBarcodeTrie);
	
	// check the validity of the barcodes in seq2_list
	// if more than 60% don't match (and mismatch) with valid barcode list, then we need to use the reverse complements
	if (!noReverseComplement) {
		int numToCheck = 10000;
		int unmatched = count_unmatched_barcodes(fq2_fn_list, validBarcodeTrie, numToCheck);
		if (unmatched >= 0.4f * (float)numToCheck) {
			Rcpp::Rcout << "Poor match between fastq barcodes and valid barcode file. Using reverse complement of found barcodes for error correction.\n";
			useReverseComplement = true;
		} else {
			Rcpp::Rcout << "High proportion of first " << numToCheck << " reads contain a valid barcode sequence.\n";
		}
	} else {
		Rcpp::Rcout << "Forced to use forward barcode sequences\n";
	}
	return validBarcodes;
}

bool check_qual(char *qual_s, int trim_n, int thr, int below_thr)
{
    int not_pass = 0;
    for (int i = 0; i < trim_n; i++)
    {
        if ((int)qual_s[i] <= thr)
        {
            not_pass++;
        }
    }
    return (not_pass > below_thr) ? false : true;
}



bool N_check(char *seq, int trim_n)
{
    bool pass = true;
    char *ptr = strchr(seq, 'N');
    if (ptr)
    {
        int index = ptr - seq;
        if (index < trim_n)
        {
            pass = false;
        }
    }
    return pass;
}



void kseq_t_to_bam_t(kseq_t *seq, bam1_t *b, int trim_n)
{
    int seq_l = seq->seq.l - trim_n; // seq length after trim the barcode
    b->l_data = seq->name.l + 1 + (int)(1.5 * seq_l + (seq_l % 2 != 0)); // +1 includes the tailing '\0'
    if ((int)b->m_data < (int)b->l_data)
    {
        b->m_data = b->l_data;
        kroundup32(b->m_data);
        b->data = (uint8_t*)realloc(b->data, b->m_data);
    }
    b->core.tid = -1;
    b->core.pos = -1;
    b->core.mtid = -1;
    b->core.mpos = -1;
    b->core.flag = BAM_FUNMAP;
    b->core.l_qname = seq->name.l + 1; // +1 includes the tailing '\0'
    b->core.l_qseq = seq_l;
    b->core.n_cigar = 0; // we have no cigar sequence
    memcpy(b->data, seq->name.s, sizeof(char) * seq->name.l); // first set qname
    b->data[seq->name.l] = '\0';
    uint8_t *s = bam_get_seq(b);
    int i = 0;
    for (i = 0; i < b->core.l_qseq; ++i) // set sequence
    {
        bam1_seq_seti(s, i, seq_nt16_table[(int)seq->seq.s[i + trim_n]]);
    }

    s = bam_get_qual(b);

    for (i = 0; i < b->core.l_qseq; ++i) // set quality
    {
        s[i] = seq->qual.s[i + trim_n]-33;
    }
}




void paired_fastq_to_bam(char *fq1_fn, char *fq2_fn, char *bam_out, const read_s read_structure, const filter_s filter_settings)
{
    // open files
    gzFile fq1 = gzopen(fq1_fn, "r"); // input fastq
    if (!fq1) {
        std::stringstream err_msg;
        err_msg << "Can't open file: " << fq1_fn << "\n";
        Rcpp::stop(err_msg.str());
    }
    gzFile fq2 = gzopen(fq2_fn, "r");
    if (!fq2) {
        std::stringstream err_msg;
        err_msg << "Can't open file: " << fq2_fn << "\n";
        Rcpp::stop(err_msg.str());
    }

    samFile *fp = sam_open(bam_out,"wb"); // output file

    // write header
    bam_hdr_t *hdr = bam_hdr_init();
    hdr->l_text = strlen(empty_header);
    hdr->text = strdup(empty_header);
    hdr->n_targets = 0;
    int hts_retcode;
    hts_retcode = sam_hdr_write(fp, hdr);
	if (hts_retcode < 0) Rcpp::Rcout << "Error in sam_hdr_write to " << bam_out << "\n";
    // get settings
    int id1_st = read_structure.id1_st;
    int id1_len = read_structure.id1_len;
    int id2_st = read_structure.id2_st;
    int id2_len = read_structure.id2_len;
    int umi_st = read_structure.umi_st;
    int umi_len = read_structure.umi_len;

    int bc1_end, bc2_end; // get total length of index + UMI for read1 and read2
    int state; // 0 for two index with umi, 1 for two index without umi, 2 for one index with umi, 3 for one index without umi

    // naming states
    const int TWO_INDEX_WITH_UMI = 0;
    const int TWO_INDEX_NO_UMI = 1;
    const int ONE_INDEX_WITH_UMI = 2;
    const int ONE_INDEX_NO_UMI = 3;

    if (id1_st >= 0) // if we have plate index
    {
        state = TWO_INDEX_WITH_UMI;
        bc1_end = id1_st + id1_len;
    }
    else // if no plate information, use id1_len to trim the read 1
    {
        state = ONE_INDEX_WITH_UMI;
        bc1_end = id1_len;
    }

    // set barcode end index
    if (umi_st >= 0)
    {
        if (id2_st + id2_len > umi_st + umi_len)
        {
            bc2_end = id2_st + id2_len;
        }
        else
        {
            bc2_end = umi_st + umi_len;
        }
    }
    else
    {
        state++; // no umi
        bc2_end = id2_st + id2_len;
    }

    // set offset for fastq header
    int name_offset;
    if (state == TWO_INDEX_WITH_UMI)
    {
        name_offset = id1_len + id2_len + umi_len + 2;
    }
    else if (state == TWO_INDEX_NO_UMI)
    {
        name_offset = id1_len + id2_len + 2;
    }
    else if (state == ONE_INDEX_WITH_UMI)
    {
        name_offset = id2_len + umi_len + 2;
    }
    else if (state == ONE_INDEX_NO_UMI)
    {
        name_offset = id2_len + 2;
    }

    // initialise fastq readers
    kseq_t *seq1;
    seq1 = kseq_init(fq1);
    kseq_t *seq2;
    seq2 = kseq_init(fq2);

    // filter tallies
    int passed_reads = 0;
    int removed_have_N = 0;
    int removed_low_qual = 0;

    int l1 = 0;
    int l2 = 0;
    // main loop, iterate through each fastq record
    // assume there are the name number of reads in read1 and read2 files, not checked.
    while (((l1 = kseq_read(seq1)) >= 0) && ((l2 = kseq_read(seq2)) >= 0))
    {
        // qual check before we do anything
        if (filter_settings.if_check_qual)
        { // Only check barcode/UMI quality
            if (!(check_qual(seq1->seq.s, bc1_end, filter_settings.min_qual, filter_settings.num_below_min) &&
                 check_qual(seq2->seq.s, bc2_end, filter_settings.min_qual, filter_settings.num_below_min)))
            {
                removed_low_qual++;
                continue;
            }
        }
        if (filter_settings.if_remove_N)
        {
            if (!(N_check(seq1->seq.s, bc1_end) && N_check(seq2->seq.s, bc2_end)))
            {
                removed_have_N++;
                continue;
            }
        }

        // begin processing valid read
        passed_reads++;

        bam1_t *b = bam_init1();

        // move original read name
        int new_name_length = name_offset + seq1->name.l;
        char* old_name_adr = seq1->name.s;
        char* new_name_adr = old_name_adr + name_offset;
        int n_char_copied = seq1->name.l * sizeof(char);
        seq1->name.s = (char*)realloc(old_name_adr, new_name_length);
        memcpy(new_name_adr, old_name_adr, n_char_copied);

        if (state == TWO_INDEX_WITH_UMI)
        {
            memcpy(seq1->name.s, seq1->seq.s + id1_st, id1_len*sizeof(char)); // copy index one
            memcpy(seq1->name.s + id1_len, seq2->seq.s + id2_st, id2_len*sizeof(char)); // copy index two
            seq1->name.s[id1_len + id2_len] = '_'; // add separator
            memcpy(seq1->name.s + id1_len + id2_len + 1, seq2->seq.s + umi_st, umi_len*sizeof(char)); // copy umi

        }
        else if (state == TWO_INDEX_NO_UMI)
        {
            memcpy(seq1->name.s, seq1->seq.s + id1_st, id1_len*sizeof(char)); // copy index one
            memcpy(seq1->name.s + id1_len, seq2->seq.s + id2_st, id2_len*sizeof(char)); // copy index two
            seq1->name.s[id1_len + id2_len] = '_'; // add separator
        }
        else if (state == ONE_INDEX_WITH_UMI)
        {
            memcpy(seq1->name.s, seq2->seq.s + id2_st, id2_len*sizeof(char)); // copy index two
            seq1->name.s[id2_len] = '_'; // add separator
            memcpy(seq1->name.s + id2_len + 1, seq2->seq.s + umi_st, umi_len*sizeof(char)); // copy umi
        }
        else if (state == ONE_INDEX_NO_UMI)
        {
            memcpy(seq1->name.s, seq2->seq.s + id2_st, id2_len*sizeof(char)); // copy index two
            seq1->name.s[id2_len] = '_'; // add separator
        }
        seq1->name.s[name_offset-1] = '#';
        seq1->name.l = name_offset + seq1->name.l;

        kseq_t_to_bam_t(seq1, b, bc1_end);
        // write bam file
        int ret = sam_write1(fp, hdr, b);
        if (ret < 0)
        {
            std::stringstream err_msg;
            err_msg << "fail to write the bam file: " << seq1->name.s << "\n";
            err_msg << "return code: " << ret << "\n";
            Rcpp::stop(err_msg.str());
        }
        bam_destroy1(b);
    }

    // cleanup
    kseq_destroy(seq1); kseq_destroy(seq2); // free seq
    gzclose(fq1); gzclose(fq2); // close fastq file
    sam_close(fp); // close bam file

    // print stats
    Rcpp::Rcout << "pass QC: " << passed_reads << "\n";
    Rcpp::Rcout << "removed_have_N: " << removed_have_N << "\n";
    Rcpp::Rcout << "removed_low_qual: " << removed_low_qual << "\n";
}




void fq_write(std::ofstream& o_stream, kseq_t *seq, int trim_n)
{
    o_stream << "@" << seq->name.s << "\n" << 
        (seq->seq.s+trim_n) << "\n" << 
        "+" << "\n" << 
        (seq->qual.s+trim_n) << "\n";
}




void fq_gz_write(gzFile out_file, kseq_t *seq, int trim_n) {
    std::stringstream stream;
    stream << "@" << seq->name.s << "\n" << 
        (seq->seq.s+trim_n) << "\n" << 
        "+" << "\n" << 
        (seq->qual.s+trim_n) << "\n";
    gzputs(out_file, stream.str().c_str());
}




void paired_fastq_to_fastq(
    char *fq1_fn,
    char *fq2_fn,
    char *fq_out,
    const read_s read_structure,
    const filter_s filter_settings,
    const bool write_gz
)
{
    int passed_reads = 0;
    int removed_have_N = 0;
    int removed_low_qual = 0;

    int l1 = 0;
    int l2 = 0;
    gzFile fq1 = gzopen(fq1_fn, "r"); // input fastq
    if (!fq1) {
        file_error(fq1_fn);
    }
    gzFile fq2 = gzopen(fq2_fn, "r");
    if (!fq2) {
        file_error(fq2_fn);
    }

    gzFile o_stream_gz;
    std::ofstream o_stream;
    if (write_gz) {
        o_stream_gz = gzopen(fq_out, "wb2"); // open gz file
        if (!o_stream_gz) {
            file_error(fq_out);
        }
    } else {
        o_stream.open(fq_out); // output file
        if (!o_stream.is_open()) {
            file_error(fq_out);
        }
    }

    // get settings
    int id1_st = read_structure.id1_st;
    int id1_len = read_structure.id1_len;
    int id2_st = read_structure.id2_st;
    int id2_len = read_structure.id2_len;
    int umi_st = read_structure.umi_st;
    int umi_len = read_structure.umi_len;

    int bc1_end, bc2_end; // get total length of index + UMI for read1 and read2
    int state; // 0 for two index with umi, 1 for two index without umi, 2 for one index with umi, 3 for one index without umi

    // naming states
    const int TWO_INDEX_WITH_UMI = 0;
    const int TWO_INDEX_NO_UMI = 1;
    const int ONE_INDEX_WITH_UMI = 2;
    const int ONE_INDEX_NO_UMI = 3;

    if (id1_st >= 0) // if we have plate index
    {
        state = TWO_INDEX_WITH_UMI;
        bc1_end = id1_st + id1_len;
    }
    else // if no plate information, use id1_len to trim the read 1
    {
        state = ONE_INDEX_WITH_UMI;
        bc1_end = id1_len;
    }

    // set barcode end index
    if (umi_st >= 0)
    {
        if (id2_st + id2_len > umi_st + umi_len)
        {
            bc2_end = id2_st + id2_len;
        }
        else
        {
            bc2_end = umi_st + umi_len;
        }
    }
    else
    {
        state++; // no umi
        bc2_end = id2_st + id2_len;
    }

    // set offset for fastq header
    int name_offset;
    if (state == TWO_INDEX_WITH_UMI)
    {
        name_offset = id1_len + id2_len + umi_len + 2;
    }
    else if (state == TWO_INDEX_NO_UMI)
    {
        name_offset = id1_len + id2_len + 2;
    }
    else if (state == ONE_INDEX_WITH_UMI)
    {
        name_offset = id2_len + umi_len + 2;
    }
    else if (state == ONE_INDEX_NO_UMI)
    {
        name_offset = id2_len + 2;
    }

    kseq_t *seq1;
    seq1 =  kseq_init(fq1);
    kseq_t *seq2;
    seq2 =  kseq_init(fq2);

    size_t _interrupt_ind = 0;
    // main loop, iter through each fastq records
    // ideally there should be equal number of reads in fq1 and fq2. we dont check this.
    while (((l1 = kseq_read(seq1)) >= 0) && ((l2 = kseq_read(seq2)) >= 0))
    {
        if (++_interrupt_ind % 4096 == 0) checkUserInterrupt();

        // validity of input parameters against length of read
        if (id2_st + id2_len > l2) continue; // check for barcode in read 2
        // check and test if barcode in read 1 is beyond read 1 length
        if ((state == TWO_INDEX_NO_UMI || state == TWO_INDEX_WITH_UMI) && (id1_st + id1_len > l1)) continue; 
        // check for a UMI if we have it
        if ((state == TWO_INDEX_WITH_UMI || state == ONE_INDEX_WITH_UMI) && (umi_st + umi_len > l2)) continue;

        // qual check before we do anything
        if (filter_settings.if_check_qual)
        {
            if (!(check_qual(seq1->qual.s, bc1_end, filter_settings.min_qual, filter_settings.num_below_min) \
                && check_qual(seq2->qual.s, bc2_end, filter_settings.min_qual, filter_settings.num_below_min)))
            {
                removed_low_qual ++;
                continue;
            }
        }
        if (filter_settings.if_remove_N)
        {
            if (!(N_check(seq1->seq.s, bc1_end) && N_check(seq2->seq.s, bc2_end)))
            {
                removed_have_N ++;
                continue;
            }
        }

        passed_reads++;

        const int new_name_length = name_offset + seq1->name.l + 1;
        seq1->name.s = (char*)realloc(seq1->name.s, new_name_length); // allocate additional memory

        const int name_size = seq1->name.l + 1; // +1 for the null byte
        char * const seq1_name = seq1->name.s;
        char * const seq1_seq = seq1->seq.s;
        char * const seq2_seq = seq2->seq.s;
        memmove(seq1_name + name_offset, seq1_name, name_size * sizeof(char)); // move original read name
        if (state == TWO_INDEX_WITH_UMI)
        {
            memcpy(seq1_name, seq1_seq + id1_st, id1_len * sizeof(char)); // copy index one
            memcpy(seq1_name + id1_len, seq2_seq + id2_st, id2_len * sizeof(char)); // copy index two
            seq1_name[id1_len+id2_len] = '_'; // add separator
            memcpy(seq1_name + id1_len + id2_len + 1, seq2_seq + umi_st, umi_len * sizeof(char)); // copy umi
        }
        else if (state == TWO_INDEX_NO_UMI)
        {
            memcpy(seq1_name, seq1_seq + id1_st, id1_len * sizeof(char)); // copy index one
            memcpy(seq1_name + id1_len, seq2_seq + id2_st, id2_len * sizeof(char)); // copy index two
            seq1_name[id1_len+id2_len] = '_'; // add separator
        }
        else if (state == ONE_INDEX_WITH_UMI)
        {
            memcpy(seq1_name, seq2_seq + id2_st, id2_len * sizeof(char)); // copy index two
            seq1_name[id2_len] = '_'; // add separator
            memcpy(seq1_name + id2_len + 1, seq2_seq + umi_st, umi_len * sizeof(char)); // copy umi
        }
        else if (state == ONE_INDEX_NO_UMI)
        {
            memcpy(seq1_name, seq2_seq + id2_st, id2_len * sizeof(char)); // copy index two
            seq1_name[id2_len] = '_'; // add separator
        }
        seq1_name[name_offset - 1] = '#';
        seq1->name.l = name_offset + seq1->name.l;

        if (write_gz) {
            fq_gz_write(o_stream_gz, seq1, bc1_end); // write to gzipped fastq file
        } else {
            fq_write(o_stream, seq1, bc1_end); // write to fastq file
        }
    }

    kseq_destroy(seq1); kseq_destroy(seq2); // free seq
    gzclose(fq1); gzclose(fq2); // close fastq file
    if (write_gz) gzclose(o_stream_gz);
    Rcpp::Rcout << "pass QC: " << passed_reads << "\n";
    Rcpp::Rcout << "removed_have_N: " << removed_have_N << "\n";
    Rcpp::Rcout << "removed_low_qual: " << removed_low_qual << "\n";
}


// Find whether kseq_t has an N in the sequence starting at startPos, for len bases
// Checks the fastq sequence for this, so this check must happen before
// the barcode and UMI sequences are removed from the fastq line
bool find_N(const kseq_t *read, size_t startPos, size_t len)
{
    std::string seq = read->seq.s;
    
    // We only care for the string before the pound sign
    std::string substr_of_interest = seq.substr(startPos, len);
    
    // If it finds an "N" in the substring of interest,
    // then the value of ".find" is different from std::string::npos
    return(substr_of_interest.find("N") != std::string::npos);
    
}

// Very similar to check_qual(). Should probably merge eventually.
bool sc_atac_check_qual(const char *qual_s, int trim_n, int thr, int below_thr){
    int not_pass = 0;
    for (int i = 0; i < trim_n; i++){
        // https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm
        if ((int)qual_s[i] - 33 <= thr){
            not_pass++;
        }
    }

    return (not_pass > below_thr) ? false : true;
}

// copy a sequence into the name part of a kseq_t
// used for taking a barcode sequence from a fastq file and prefixing each name in another fastq file
void copySequenceIntoKseqName(kseq_t *seqNameDestination, const char *seqToCopy, size_t copyLength) {
	const int new_name_length = seqNameDestination->name.l + copyLength + 1;
	// allocate additional memory and move original read name

	std::string newName;
	newName.reserve(new_name_length);
	newName += seqToCopy;
	newName += "#";
	newName += seqNameDestination->name.s;

	// copy the new name back into the kseq_t
	free(seqNameDestination->name.s);
	seqNameDestination->name.l = new_name_length;
	seqNameDestination->name.s = (char *)malloc((newName.size() + 1) * sizeof(char));
	strcpy(seqNameDestination->name.s, newName.c_str());
	
	// old code
	// seqNameDestination->name.s = (char *)realloc(seqNameDestination->name.s, new_name_length);
	// memmove(seqNameDestination->name.s + copyLength + 1, seqNameDestination->name.s, seqNameDestination->name.l * sizeof(char));
	// // copy the sequence into the new name array
	// memcpy(seqNameDestination->name.s, seqToCopy, copyLength * sizeof(char));
	// // add separator
	// seqNameDestination->name.s[copyLength] = '#';
	// // is this necessary?
	// seqNameDestination->name.s[new_name_length] = '\0';
}

// sc_atac_paired_fastq_to_fastq ------------------

std::vector<int> sc_atac_paired_fastq_to_fastq(
        const char *fq1_fn,
        std::vector<std::string> fq2_fn_list,
        const char *fq3_fn,
		const char *valid_barcode_fn,
        const char *fq_out,
        const bool write_gz,
        const bool rmN,
        const bool rmlow,
        int min_qual,
        int num_below_min,
		bool no_reverse_complement
) {
    // // get rid of kseq.h warnings
    // REMOVE_KSEQ_WARNINGS();

    std::vector<int> out_vect(7, 0);  // output vector of length 7 filled with zeroes
    // Input parameters when rmlow is true
    // int min_qual = 20; // minq: the minimum base pair quality that we allowed (from scPipe wrapper_scPipeCPP.R)
    // int num_below_min = 2; // numbq: the maximum number of base pair that have quality (from scPipe wrapper_scPipeCPP.R)
	// (these are unused)
    // int id1_st = 0; // bs1: starting position of barcode in read one. -1 if no barcode in read one.
    // int id1_len = 50; // bl1: length of barcode in read one, if there is no barcode in read one this number is used for trimming beginning of read one.
    // int id2_st = 0; // bs2: starting position of barcode in read two
    // int id2_len = 49; // bl2: length of barcode in read two
    // int umi_st = 0; // us: starting position of UMI // umi_start
    // int umi_len = 0; // ul: length of UMI // umi_length
    
    int passed_reads = 0, removed_Ns = 0, removed_low_qual = 0;
    
    bool R3 = strcmp(fq3_fn, "") != 0; // R3 exists if fq3_fn is not empty
    
    int l1 = 0, l2 = 0, l3 = 0;
    gzFile fq1 = gzopen(fq1_fn, "r"); // input fastq
	gzFile fq3;
    if (!fq1) {
        file_error(fq1_fn);
    }
    
	// preprocess the fastq barcode file(s) so that we have a vector of gzFiles to deal with
    std::vector<gzFile> fq2_list = open_gz_files(fq2_fn_list);
	std::vector<kseq_t*> seq2_list;
	for (int i = 0; i < (int)fq2_list.size(); i++) {
		seq2_list.push_back(kseq_init(fq2_list[i]));
	}
    

	gzFile o_stream_gz_R1;
    std::ofstream o_stream_R1;
	gzFile o_stream_gz_R3;
    std::ofstream o_stream_R3;

	gzFile o_stream_gz_R1_Partial;
    std::ofstream o_stream_R1_Partial;
    gzFile o_stream_gz_R3_Partial;
    std::ofstream o_stream_R3_Partial;

	gzFile o_stream_gz_R1_No;
    std::ofstream o_stream_R1_No;
    gzFile o_stream_gz_R3_No;
    std::ofstream o_stream_R3_No;
    
	kseq_t *seq1 = kseq_init(fq1);;
    kseq_t *seq3;
	
    const char* appendCompleteMatch = "/demux_completematch_";
    const char* appendPartialMatch = "/demux_partialmatch_";
    const char* appendNoMatch = "/demux_nomatch_";

    char *fqoutR1 = createFileWithAppend(fq_out, appendCompleteMatch, fq1_fn);
    openFile(o_stream_gz_R1, o_stream_R1,fqoutR1, write_gz);
    
    char *fqoutR1Partial = createFileWithAppend(fq_out, appendPartialMatch, fq1_fn);
    openFile(o_stream_gz_R1_Partial, o_stream_R1_Partial, fqoutR1Partial, write_gz);
    
    char *fqoutR1No = createFileWithAppend(fq_out, appendNoMatch, fq1_fn);
    openFile(o_stream_gz_R1_No, o_stream_R1_No,fqoutR1No, write_gz);

 	if (R3) {
        fq3 = gzopen(fq3_fn, "r");
        if (!fq3) {
            file_error(fq3_fn);
        }
        
        char *fqoutR3 = createFileWithAppend(fq_out,appendCompleteMatch,fq3_fn);
        openFile(o_stream_gz_R3,o_stream_R3,fqoutR3, write_gz);
        
        char *fqoutR3Partial = createFileWithAppend(fq_out,appendPartialMatch,fq3_fn);
        openFile(o_stream_gz_R3_Partial,o_stream_R3_Partial,fqoutR3Partial, write_gz);
        
        char *fqoutR3No = createFileWithAppend(fq_out,appendNoMatch,fq3_fn);
        openFile(o_stream_gz_R3_No,o_stream_R3_No,fqoutR3No, write_gz);
        
        seq3 =  kseq_init(fq3);
    }
    
	// preprocess the valid_barcode_file (csv) reads into a Trie, allowing for matching and mismatching.
	bool checkBarcodeMismatch = false;
	std::vector<std::string> validBarcodes;
	Trie validBarcodeTrie;
	bool useReverseComplement = false;
	if (std::string(valid_barcode_fn) != "") {
		checkBarcodeMismatch = true;
		validBarcodes = processAndBuildTrie(valid_barcode_fn, fq2_fn_list, validBarcodeTrie, no_reverse_complement, useReverseComplement);
	} else {
		Rcpp::Rcout << "No valid_barcode_file provided; no barcode error correction will occur.\n";
	}

    std::set<std::string> seq_2_set; // Set that will include the unique barcode sequences
    MatchType barcodeMatchType = Exact;
	int exactMatches = 0, partialMatches = 0, noMatches = 0;
    size_t _interrupt_ind = 0;
	// for output files;
    gzFile *R1_gz_outfile;
    gzFile  *R3_gz_outfile;
    std::ofstream *R1_outfile;
    std::ofstream *R3_outfile;
    // Assuming R1, R2, R3 all are of equal lengths.
	// int tmpC = 0;
    while ((l1 = kseq_read(seq1)) >= 0) {
		// if (tmpC++ > 10000) break;
		// Rcpp::Rcout << tmpC << ": \n\t" << seq1->name.s << "\n";
        if (++_interrupt_ind % 4096 == 0) checkUserInterrupt();
        passed_reads++;
        
        if (R3){
            if((l3 = kseq_read(seq3)) < 0){
				Rcpp::Rcout << "read2 file is not of the same length as the barcode fastq file:.\n";
				break;
            }
        } 
        
		// append every barcode in the barcode fastq list 
		// to the correct position in the name of the current read (in both R1 and R3)
		int barcodeExact = 0, barcodePartial = 0, barcodeNo = 0;
		bool removedAnN = false, removedLowQual = false;
        for (int i = 0; i < (int)seq2_list.size(); i++) {
			// grab the current kseq_t fastq read
            kseq_t* seq2 = seq2_list[i];
            if ((l2 = kseq_read(seq2)) < 0) {
				Rcpp::Rcout << "read1 file is not the same length as the barcode fastq file.\n";
				break;
			}
			
			// nonconst pointer to const char, important as we redirect the pointer
			const char * seq2_seq = seq2->seq.s; // if we match to a new barcode when doing barcode mismatch
			
			seq_2_set.insert(std::string{seq2_seq});

			// quality check before the heavy work of shifting barcodes around is done. Early break
			if(rmlow) {
				// Check quality for the entire read (using kstring_t length)
				if (!sc_atac_check_qual(seq2->qual.s, seq2->qual.l, min_qual, num_below_min)) {
					removedLowQual = true;
					break;
				} 
			}

			// check for N values in this current barcode from the fastq
			// Check for N values before doing any processing to save time
			if(rmN){
				// If find_N is TRUE, then there is an N in the sequence
				if(find_N(seq2, 0, seq2->seq.l)) {
					removedAnN = true;
					break; // ignore the entire read
				} 
			} 
			
			// verify that we have the correct  barcode sequence, by checking against the valid barcode trie
			if (checkBarcodeMismatch) {
				std::string seq2_seq_str = 
					useReverseComplement ? 
						reverseComplement(seq2_seq, seq2->seq.l) :
						seq2_seq;
				
				std::vector<MismatchResult> possibleBarcodes = validBarcodeTrie.Locate_Seq_Mismatches(seq2_seq_str, 0, seq2_seq_str.size());
				// // find the barcode which matches best (highest chance of matching perfectly based quality score)
				int highestScore = 0;
				int barcodePosition = -1;
				bool exactMatch = false;
				for (const MismatchResult &match : possibleBarcodes) {
					if (match.mismatchPosition == -1) {
						// we've found a perfect match
						// ignore all other matches
						barcodePosition = match.sequenceIndex;
						exactMatch = true;
						break;
					}
					int thisQual = (int)seq2->qual.s[match.mismatchPosition] - 33; // fastq qscores are ASCII 33 to 126
					if (thisQual > highestScore) {
						highestScore = thisQual;
						barcodePosition = match.sequenceIndex;
					}
				}

				if (barcodePosition != -1) {
					// we've found a better barcode match
					seq2_seq = validBarcodes[barcodePosition].c_str();
					exactMatch ?
						barcodeExact++ :
						barcodePartial++;
				} else {
					barcodeNo++;
				}
			} else { 
				// if there's no valid barcode file
				barcodeExact++;
			}

			copySequenceIntoKseqName(seq1, seq2_seq, seq2->seq.l);

			if (R3) {
				copySequenceIntoKseqName(seq3, seq2_seq, seq2->seq.l);
			}
            
		}

		// if we've removed a low quality barcode sequence,
		// we must ignore the entire read
		if (removedLowQual) {
			removed_low_qual++;
			continue;
		}
		// if we've removed an N value from a barcode sequence,
		// ignore the entire read as part of quality checking
		if (removedAnN) {
			// If there was an N in the sequence, then
			// we add 1 to the counter of reads deleted and 
			// skip the rest of the code in the loop.
			removed_Ns++; // Add 1 to the counter of reads deleted
			continue;
		}

		// decide the type of match
		// if all n barcode files give an exact match, match type is exact
		// if 1<=x<n are exact or partial matches, match type is partial
		// if all n barcodes are no match, no match
		if (barcodeExact == (int)seq2_list.size()) {
			barcodeMatchType = Exact;
			exactMatches++;
		} else if (barcodePartial + barcodeExact >= 1) {
			barcodeMatchType = Partial;
			partialMatches++;
		} else {
			barcodeMatchType = NoMatch;
			noMatches++;
		}

		// write to R1 output file	
		switch (barcodeMatchType) {
            case Exact:
                R1_outfile = &o_stream_R1;
                R1_gz_outfile = &o_stream_gz_R1;
				R3_outfile = &o_stream_R3;
				R3_gz_outfile = &o_stream_gz_R3;
                break;
            case Partial:
                R1_outfile = &o_stream_R1_Partial;
                R1_gz_outfile = &o_stream_gz_R1_Partial;
				R3_outfile = &o_stream_R3_Partial;
				R3_gz_outfile = &o_stream_gz_R3_Partial;
                break;
            case NoMatch:
            default:
                R1_outfile = &o_stream_R1_No;
                R1_gz_outfile = &o_stream_gz_R1_No;
				R3_outfile = &o_stream_R3_No;
				R3_gz_outfile = &o_stream_gz_R3_No;
                break;
        }

        if (write_gz) {
            fq_gz_write(*R1_gz_outfile, seq1, 0); // write to gzipped fastq file
			if (R3) {
				fq_gz_write(*R3_gz_outfile, seq3, 0);
			}
        } else {
            fq_write(*R1_outfile, seq1, 0); // write to fastq file
			if (R3) {
				fq_write(*R3_outfile, seq3, 0);
			}
        }
    } // end while
    
    // free subStr1
    kseq_destroy(seq1); 
    for(int i=0;i<(int)seq2_list.size();i++){
        kseq_t* seq2 = seq2_list[i];
        kseq_destroy(seq2);// free seq
        gzclose(fq2_list[i]);
    }
    
	gzclose(fq1); // close fastq file
    if(R3){
        kseq_destroy(seq3);
        gzclose(fq3);
    }
    if (write_gz){
        gzclose(o_stream_gz_R1);
        gzclose(o_stream_gz_R1_Partial);
        gzclose(o_stream_gz_R1_No);
        
        if(R3){
            gzclose(o_stream_gz_R3);
            gzclose(o_stream_gz_R3_Partial);
            gzclose(o_stream_gz_R3_No);
		}
    } else {
		o_stream_R1.close();
		o_stream_R1_Partial.close();
		o_stream_R1_No.close();
		
		if (R3) {
			o_stream_R3.close();
			o_stream_R3_Partial.close();
			o_stream_R3_No.close();
		}
	}

    Rcpp::Rcout << "Total reads: " << passed_reads << "\n";
    Rcpp::Rcout << "Total reads removed due to N's in barcodes: " << removed_Ns << "\n";
    Rcpp::Rcout << "Total reads removed due to low quality barcodes: " << removed_low_qual << "\n";
    Rcpp::Rcout << "Total barcodes provided in FASTQ file: " << seq_2_set.size() << "\n";
	if (checkBarcodeMismatch) {
		Rcpp::Rcout << "Exact barcode matches: " << exactMatches << "\n";
		Rcpp::Rcout << "Matched after barcode correction : " << partialMatches << "\n";
		Rcpp::Rcout << "No barcode matches: " << noMatches << "\n";
	}
    
    out_vect[0] = passed_reads;
    out_vect[1] = removed_Ns;
    out_vect[2] = removed_low_qual;
    out_vect[3] = seq_2_set.size();
    out_vect[4] = exactMatches;
    out_vect[5] = partialMatches;
    out_vect[6] = noMatches;
    
    return(out_vect);
}




// sc_atac_paired_fastq_to_csv ------------------

std::vector<int> sc_atac_paired_fastq_to_csv(
        char *fq1_fn,
        char *fq3_fn,
        char *fq_out, 
        char *bc_fn, // barcode file must be a file where each line is a barcode (not comma separated)
		char *valid_barcode_fn,
        int umi_start,
        int umi_length,
        char *umi_in,
        const bool write_gz,
        const bool rmN,
        const bool rmlow,
        int min_qual,
        int num_below_min,
        int id1_st,
        int id1_len,
        int id2_st,
        int id2_len
)
{
    std::vector<int> out_vect(6, 0);  // output vector of length 6 filled with zeroes
    
    int passed_reads = 0;
    int removed_Ns = 0;
    int removed_low_qual = 0;
    int exact_match = 0;
    int approx_match = 0;
    
    bool R3 = false;
    bool isUMIR1 = (umi_length > 0 && (strcmp(umi_in,"both") == 0 || strcmp(umi_in,"R1") == 0));
    bool isUMIR2 = (umi_length > 0 && (strcmp(umi_in,"both") == 0 || strcmp(umi_in,"R2") == 0));
    //bool write_line;
    
    if(!strcmp(fq3_fn,"")){
        R3 = false;
    }else{
        R3 = true;
    }
    
    int l1 = 0;
    //int l2 = 0;
    int l3 = 0;
    gzFile fq1 = gzopen(fq1_fn, "r"); // input fastq
    if (!fq1) {
        file_error(fq1_fn);
    }
    
    std::map<std::string, int> barcode_map; 
    std::ifstream bc(bc_fn);
    std::string line;
    
	while(std::getline(bc, line))
	{
		std::string substr = line.substr(0,std::min(id1_len, id2_len));
		barcode_map.insert(std::pair<std::string, int>(substr, 1)); 
	}
	
	if(barcode_map.empty()){
		std::stringstream err_msg;
		err_msg << "Error in retrieving barcodes from the barcode File. Please check the barcode file format. " << bc_fn << std::endl;
		Rcpp::stop(err_msg.str());
	}
    
	// preprocess the valid_barcode_file (csv) reads into a Trie, allowing for matching and mismatching.
	bool checkBarcodeMismatch = false;
	std::vector<std::string> validBarcodes;
	Trie validBarcodeTrie;
	int noValidBcMatch = 0;
	if (std::string(valid_barcode_fn) != "") {
		checkBarcodeMismatch = true;
		validBarcodes = readBarcodes(valid_barcode_fn);
		preprocessBarcodes(validBarcodes, validBarcodeTrie);
	} else {
		Rcpp::Rcout << "No valid_barcode_file provided; no barcode error correction will occur.\n";
	}
    
    
    gzFile fq3;
    gzFile o_stream_gz_R3;
    std::ofstream o_stream_R3;
    gzFile o_stream_gz_R3_Partial;
    std::ofstream o_stream_R3_Partial;
    gzFile o_stream_gz_R3_No;
    std::ofstream o_stream_R3_No;
    
    kseq_t *seq3;
    const char* appendCompleteMatch = "/demux_completematch_";
    const char* appendPartialMatch = "/demux_partialmatch_";
    const char* appendNoMatch = "/demux_nomatch_";
    
    if(R3){
        fq3 = gzopen(fq3_fn, "r");
        if (!fq3) {
            file_error(fq3_fn);
        }
        
        char *fqoutR3 = createFileWithAppend(fq_out,appendCompleteMatch,fq3_fn);
        openFile(o_stream_gz_R3,o_stream_R3,fqoutR3, write_gz);
        
        char *fqoutR3Partial = createFileWithAppend(fq_out,appendPartialMatch,fq3_fn);
        openFile(o_stream_gz_R3_Partial,o_stream_R3_Partial,fqoutR3Partial, write_gz);
        
        char *fqoutR3No = createFileWithAppend(fq_out,appendNoMatch,fq3_fn);
        openFile(o_stream_gz_R3_No,o_stream_R3_No,fqoutR3No, write_gz);
        
        seq3 =  kseq_init(fq3);
        
    }
    
    kseq_t *seq1;
    seq1 =  kseq_init(fq1);
    
    
    gzFile o_stream_gz_R1;
    std::ofstream o_stream_R1;
    char *fqoutR1 = createFileWithAppend(fq_out,appendCompleteMatch,fq1_fn);
    openFile(o_stream_gz_R1,o_stream_R1,fqoutR1, write_gz);
    
    gzFile o_stream_gz_R1_Partial;
    std::ofstream o_stream_R1_Partial;
    char *fqoutR1Partial = createFileWithAppend(fq_out,appendPartialMatch,fq1_fn);
    openFile(o_stream_gz_R1_Partial,o_stream_R1_Partial,fqoutR1Partial, write_gz);
    
    
    gzFile o_stream_gz_R1_No;
    std::ofstream o_stream_R1_No;
    char *fqoutR1No = createFileWithAppend(fq_out,appendNoMatch,fq1_fn);
    openFile(o_stream_gz_R1_No,o_stream_R1_No,fqoutR1No, write_gz);
    
    
    // Define some variables.
    int bc1_end, bc2_end; // get end position in the read of barcode and umi
    
    if (isUMIR1) {
        // 
        bc2_end = std::max(id1_st + id2_len, umi_start + umi_length);
    } else {
        if (id1_st >= 0) { // if we have plate index
            bc1_end = id1_st + id1_len;
        }
        else { // if no plate information, use id1_len to trim the read 1
            bc1_end = id1_len;
        }
    }
    
    if (isUMIR2) {
        // set barcode end index
        // This is basically doing the following: bc2_end = max(id2_st + id2_len, umi_start + umi_length)
        bc2_end = std::max(id2_st + id2_len, umi_start + umi_length);
    } else {
        bc2_end = id2_st + id2_len;
    }
    
    std::set<std::string> seq_2_set; // Set that will include the unique barcode sequences
    size_t _interrupt_ind = 0;

    // for output files;
    gzFile *R1_gz_outfile;
    gzFile  *R3_gz_outfile;
    std::ofstream *R1_outfile;
    std::ofstream *R3_outfile;
    // Assuming R1, R2, R3 all are of equal lengths.
    while (((l1 = kseq_read(seq1)) >= 0))
    {
        if (++_interrupt_ind % 50 == 0) checkUserInterrupt();
        passed_reads++;

        if(passed_reads % 1000==0) {
            Rcout << passed_reads << " lines have been read..." << std::endl; 
        }
        
        // check if this read is long enough for input params
        if (bc1_end > (int)seq1->seq.l) {
            Rcpp::stop("Read not long enough to support bc1 length");
        }

        if (R3){
            if((l3 = kseq_read(seq3)) < 0){
                Rcpp::Rcout << "R3 file is not of same length as R2: " << std::endl;
            } else {
                if (bc2_end > (int)seq3->seq.l) {
                    Rcpp::stop("Read not long enough to support bc2 length");
                }
            }
        } 


		// quality control checks for this read 1 and read 2 (if R3)
        if(rmlow) { // Only check barcode/UMI quality
            // Check quality
            if (!sc_atac_check_qual(seq1->qual.s, bc1_end, min_qual, num_below_min)) {
                removed_low_qual++;
                continue; // the rest of the lines in the while loop are ignored
            } 

            if (R3) {
                // Check quality
                if (!sc_atac_check_qual(seq3->qual.s, bc2_end, min_qual, num_below_min)) {
                    removed_low_qual++;
                    continue; // the rest of the lines in the while loop are ignored
                } 
            }
        } 
                    
        if(rmN){
            // If find_N is TRUE, then there is an N in the sequence
            if(find_N(seq1, id1_st, id1_len)){
                // If there was an N in the sequence, then
                // we add 1 to the counter of reads deleted and
                // skip the rest of the code in the loop.
                removed_Ns++; // Add 1 to the counter of reads deleted
                continue; // the rest of the lines in the while loop are ignored
            } else if (isUMIR1 && find_N(seq1, umi_start, umi_length)) {
				removed_Ns++;
				continue;
			}

            if (R3) {
                if(find_N(seq3, id2_st, id2_len)){
                    // If there was an N in the sequence, then
                    // we add 1 to the counter of reads deleted and
                    // skip the rest of the code in the loop.
                    removed_Ns++; // Add 1 to the counter of reads deleted
                    continue; // the rest of the lines in the while loop are ignored
                } else if (isUMIR2 && find_N(seq3, umi_start, umi_length)) {
					// remove this read if there's Ns in the UMI sequence
					removed_Ns++;
					continue;
				}
            }
        } // end if(rmN)

		// allocate and copy just the barcode, to check against the barcodes in the barcode map
        char *barcode = (char *)malloc((id1_len + 1) * sizeof(char));
        memcpy( barcode, seq1->seq.s + id1_st, id1_len); 
        barcode[id1_len] = '\0'; 
        
        MatchType r1_match_type = Exact; // 0 for exact, 1 for partial, 2 for no
		// we want to check for an exact match first
		if (barcode_map.find(barcode) != barcode_map.end()) {
			// exact match
			exact_match ++;   
			r1_match_type = Exact;
		} else {
			// inexact match, we need to iterate over all barcodes
			r1_match_type = NoMatch; // if we never find a barcode, no match
			for (std::map<std::string,int>::iterator it=barcode_map.begin(); it!=barcode_map.end(); ++it){
				if(hamming_distance(it->first, barcode) <2){
					approx_match++;
					r1_match_type = Partial;
					break;
				} 
			}
		}

		// if we need to check barcodes against a valid barcode file
		if (checkBarcodeMismatch) {
			std::vector<MismatchResult> possibleBarcodes = validBarcodeTrie.Locate_Seq_Mismatches(barcode, 0, id1_len);
			int barcodePosition = -1;
			bool exactMatch = false;
			if (possibleBarcodes.size() > 0) {
				for (const MismatchResult &match : possibleBarcodes) {
					if (match.mismatchPosition == -1) {
						// we've found a perfect match
						// ignore all other matches
						barcodePosition = match.sequenceIndex;
						exactMatch = true;
						break;
					}
				}
				// if we didn't find an exact match, just use the first mismatch
				if (!exactMatch && possibleBarcodes.size() > 0) {
					barcodePosition = possibleBarcodes[0].sequenceIndex;
				}

				if (barcodePosition != -1) {
					// we've found a better barcode match
					strcpy(barcode, validBarcodes[barcodePosition].c_str());
				} 
			} 

			if (barcodePosition == -1) {
				r1_match_type = NoMatch;
				noValidBcMatch++;
			}
			
		}

        // allocate space and copy the sequence from the read containing the barcode and UMI (if applicable)
        char* subStr1;
        int bcUMIlen1 = 0;
        if(isUMIR1){ // create subStr1 of barcode concatenated with UMI, seperated by '_'
            bcUMIlen1 = id1_len + umi_length + 1;
            subStr1 = (char *)malloc(bcUMIlen1 + 1); // additional space for zero terminator
            memcpy(subStr1, barcode, id1_len);
            subStr1[id1_len] = '_'; 
            memcpy(subStr1 + id1_len + 1, seq1->seq.s + umi_start, umi_length); 
            subStr1[bcUMIlen1] = '\0';
        } else { // create subStr1 of barcode sequence
            bcUMIlen1 = id1_len;
            subStr1 = (char *)malloc(bcUMIlen1 + 1); // additional space for zero terminator
            memcpy(subStr1, barcode, id1_len); 
            subStr1[bcUMIlen1] = '\0';
        }
        


        // if the barcode matches exactly or inexactly, write the modifiyed sequence lines to the output file
        const int new_name_length1 = seq1->name.l + bcUMIlen1 + 1;
        seq1->name.s = (char*)realloc(seq1->name.s, new_name_length1 + 1); // allocate additional memory
        memmove(seq1->name.s + bcUMIlen1 + 1, seq1->name.s, seq1->name.l);// move original read name from second arg to first arg
        memcpy(seq1->name.s, subStr1, bcUMIlen1 * sizeof(char)); // copy index one
        seq1->name.s[bcUMIlen1] = '#'; // add separator
        seq1->name.s[new_name_length1] = '\0'; // should we really add a 0 terminator?? this will get added every time, for every barcode
        
        // chop read to only be what has not been copied to header
        memmove(seq1->seq.s, seq1->seq.s + bc1_end, seq1->name.l - bc1_end);   
        
        if(R3) {
			char *barcode3 = (char *)malloc((id2_len + 1) * sizeof(char));
            memcpy( barcode3, seq3->seq.s + id2_st, id2_len);
            barcode3[id2_len] = '\0';

			MatchType r2_match_type = Exact; // 0 for exact, 1 for partial, 2 for no
			if (barcode_map.find(barcode3) != barcode_map.end()) {
				r2_match_type = Exact;
			} else {
				r2_match_type = NoMatch; // assume there is no match
				for (std::map<std::string,int>::iterator it=barcode_map.begin(); it!=barcode_map.end(); ++it){
					if(hamming_distance(it->first, barcode3) < 2){
						r2_match_type = Partial;
						break;
					}
				}
			}

			// if we need to check barcodes against a valid barcode file
			if (checkBarcodeMismatch) {
				std::vector<MismatchResult> possibleBarcodes3 = validBarcodeTrie.Locate_Seq_Mismatches(barcode3, 0, id1_len);
				// // find the barcode which matches best (highest chance of matching perfectly based on quality score)
				int barcodePosition3 = -1;
				bool exactMatch3 = false;
				if (possibleBarcodes3.size() > 0) {
					for (const MismatchResult &match3 : possibleBarcodes3) {
						if (match3.mismatchPosition == -1) {
							// we've found a perfect match
							// ignore all other matches
							barcodePosition3 = match3.sequenceIndex;
							exactMatch3 = true;
							break;
						}
					}

					if (!exactMatch3 && possibleBarcodes3.size() > 0) {
						barcodePosition3 = possibleBarcodes3[0].sequenceIndex;
					}
					if (barcodePosition3 != -1) {
						// we've found a better barcode match
						strcpy(barcode3, validBarcodes[barcodePosition3].c_str());
					} 
				} 
				if (barcodePosition3 == -1) {
					r2_match_type = NoMatch;
					noValidBcMatch++;
				}
			}

            // allocate and copy second barcode and umi (if applicable) to subStr3, to copy to header
            char* subStr3;
            int bcUMIlen3 = 0;
            if (isUMIR2) {
                bcUMIlen3 = id2_len + umi_length +1;
                subStr3 = (char*)malloc(bcUMIlen3 + 1); // additional space for zero terminator
                memcpy(subStr3, barcode3, id2_len);
                subStr3 [id2_len] = '_';
                memcpy(subStr3 + id2_len + 1, seq3->seq.s + umi_start, umi_length);
                subStr3[bcUMIlen3] = '\0';          
            } else {
                bcUMIlen3 = id2_len;
                subStr3 = (char*)malloc(bcUMIlen3 + 1); // additional space for zero terminator
                memcpy( subStr3, barcode3, id2_len );
                subStr3[bcUMIlen3] = '\0';
            }
            
            

            const int new_name_length1 = seq3->name.l + bcUMIlen3 + 1;
            seq3->name.s = (char*)realloc(seq3->name.s, new_name_length1 + 1); // allocate additional memory
            memmove(seq3->name.s + bcUMIlen3+1, seq3->name.s, seq3->name.l);// move original read name
            memcpy(seq3->name.s, subStr3, bcUMIlen3 * sizeof(char)); // copy index one
            seq3->name.s[bcUMIlen3] = '#'; // add separator
            seq3->name.s[new_name_length1] = '\0';
            
            memmove(seq3->seq.s, seq3->seq.s + bc2_end, seq3->seq.l - bc2_end);           
            //seq3->seq.s[seq3->seq.l - bc2_end];

			// both reads must belong to the same match type
			// so both reads get output to the same category
			if (r1_match_type > r2_match_type) {
				r2_match_type = r1_match_type;
			} else { // r2_match_type >= r1_match_type
				r1_match_type = r2_match_type;
			}

            switch (r2_match_type) {
                case Exact:
                    R3_outfile = &o_stream_R3;
                    R3_gz_outfile = &o_stream_gz_R3;
                    break;
                case Partial:
                    R3_outfile = &o_stream_R3_Partial;
                    R3_gz_outfile = &o_stream_gz_R3_Partial;
                    break;
                case NoMatch:
                default:
                    R3_outfile = &o_stream_R3_No;
                    R3_gz_outfile = &o_stream_gz_R3_No;
                    break;
            }

            if (write_gz) {
                fq_gz_write(*R3_gz_outfile, seq3, 0); // write to gzipped fastq file
            } else {
                fq_write(*R3_outfile, seq3, 0); // write to fastq file
            }
            free(barcode3);
            free(subStr3);
            
        }
    
    
		// write to R1 output file
		// if R3 exists, r1_match_type will have been adjusted to ensure both reads are
		// output to the same file, which is determined by the best match
		switch (r1_match_type) {
            case Exact:
                R1_outfile = &o_stream_R1;
                R1_gz_outfile = &o_stream_gz_R1;
                break;
            case Partial:
                R1_outfile = &o_stream_R1_Partial;
                R1_gz_outfile = &o_stream_gz_R1_Partial;
                break;
            case NoMatch:
            default:
                R1_outfile = &o_stream_R1_No;
                R1_gz_outfile = &o_stream_gz_R1_No;
                break;
        }

        if (write_gz) {
            fq_gz_write(*R1_gz_outfile, seq1, 0); // write to gzipped fastq file
        } else {
            fq_write(*R1_outfile, seq1, 0); // write to fastq file
        }

        free(barcode);
        free(subStr1);
        
    }
    
    kseq_destroy(seq1); 
    
    bc.close();
    gzclose(fq1); // close fastq file
    if(R3){
        kseq_destroy(seq3);
        gzclose(fq3);
    }
    if (write_gz){ // why do we only close the file if writing gz? shouldn't we still close if it's not writing gz?
        gzclose(o_stream_gz_R1);
        gzclose(o_stream_gz_R1_Partial);
        gzclose(o_stream_gz_R1_No);
        
        if(R3){
            gzclose(o_stream_gz_R3);
            gzclose(o_stream_gz_R3_Partial);
            gzclose(o_stream_gz_R3_No);
            
        }
    }
    
    out_vect[0] = passed_reads;
    out_vect[1] = removed_Ns;
    out_vect[2] = removed_low_qual;
    out_vect[3] = exact_match;
    out_vect[4] = approx_match;
    out_vect[5] = (int)barcode_map.size();
    
    Rcpp::Rcout << "Total Reads: " << passed_reads << std::endl;
    Rcpp::Rcout << "Total reads removed due to N's in barcodes: " << removed_Ns << std::endl;
    Rcpp::Rcout << "Total reads removed due to low quality barcodes: " << removed_low_qual << std::endl;
    Rcpp::Rcout << "Exact match Reads: " << exact_match << std::endl;
    Rcpp::Rcout << "Matched after barcode corrections: " << approx_match << std::endl;
	if (checkBarcodeMismatch) {
		Rcpp::Rcout << "Total reads removed due to no matching barcode in valid_barcode_file: " << noValidBcMatch << std::endl;
	}
    Rcpp::Rcout << "Total barcodes provided in CSV file: " << (int)barcode_map.size() << std::endl;
    
    return(out_vect);
}
