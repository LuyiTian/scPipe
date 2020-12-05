//trim_barcode
#include "trimbarcode.h"
#include <string.h>

using namespace Rcpp;



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
    if (b->m_data < b->l_data)
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










// Find whether kseq_t has an N before the pound (#) symbol
// Very similar to N_check(). Should probably merge eventually.
bool find_N(kseq_t *seq)
{
    std::string name = seq->name.s;
    
    // We only care for the string before the pound sign
    std::string substr_of_interest = name.substr(0, name.find("#"));
    
    // If it finds an "N" in the substring of interest,
    // then the value of ".find" is different from std::string::npos
    return(substr_of_interest.find("N") != std::string::npos);
    
}




// Very similar to check_qual(). Should probably merge eventually.
bool sc_atac_check_qual(char *qual_s, int trim_n, int thr, int below_thr){
    int not_pass = 0;
    for (int i = 0; i < trim_n; i++){
        // https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm
        if ((int)qual_s[i] - 33 <= thr){
            not_pass++;
        }
    }
    return (not_pass > below_thr) ? false : true;
}


// sc_atac_paired_fastq_to_fastq ------------------

void sc_atac_paired_fastq_to_fastq(
        char *fq1_fn,
        std::vector<std::string> fq2_fn_list,
        char *fq3_fn,
        char *fq_out,
        const bool write_gz,
        const bool rmN,
        const bool rmlow,
        int min_qual,
        int num_below_min,
        int id1_st,
        int id1_len,
        int id2_st,
        int id2_len,
        int umi_st,
        int umi_len
) {
    
    
    // Input parameters when rmlow is true
    // int min_qual = 20; // minq: the minimum base pair quality that we allowed (from scPipe wrapper_scPipeCPP.R)
    // int num_below_min = 2; // numbq: the maximum number of base pair that have quality (from scPipe wrapper_scPipeCPP.R)
    // int id1_st = 0; // bs1: starting position of barcode in read one. -1 if no barcode in read one.
    // int id1_len = 50; // bl1: length of barcode in read one, if there is no barcode in read one this number is used for trimming beginning of read one.
    // int id2_st = 0; // bs2: starting position of barcode in read two
    // int id2_len = 49; // bl2: length of barcode in read two
    // int umi_st = 0; // us: starting position of UMI // umi_start
    // int umi_len = 0; // ul: length of UMI // umi_length
    
    int passed_reads = 0;
    int removed_Ns = 0;
    int removed_low_qual = 0;
    
    bool R3 = false;
    
    if(!strcmp(fq3_fn,"")){
        R3 = false;
    }else{
        R3 = true;
    }
    
    int l1 = 0;
    int l2 = 0;
    int l3 = 0;
    gzFile fq1 = gzopen(fq1_fn, "r"); // input fastq
    if (!fq1) {
        file_error(fq1_fn);
    }
    
    std::vector<kseq_t*> seq2_list;
    std::vector<gzFile> fq2_list;
    for(int i=0;i<fq2_fn_list.size();i++){
        char* fq2_fn = (char *)fq2_fn_list[i].c_str();
        gzFile fq2 = gzopen(fq2_fn, "r");
        if (!fq2) {
            file_error(fq2_fn);
        }
        fq2_list.push_back(fq2);
        seq2_list.push_back(kseq_init(fq2));
    }
    
    gzFile o_stream_gz_R1;
    std::ofstream o_stream_R1;
    
    gzFile fq3;
    gzFile o_stream_gz_R3;
    std::ofstream o_stream_R3;
    kseq_t *seq3;
    
    if(R3){
        fq3 = gzopen(fq3_fn, "r");
        if (!fq3) {
            file_error(fq3_fn);
        }
        
        const char* appendR3 = "/demux_" ;
        char *fqoutR3 = (char*)malloc(strlen(fq_out)+ strlen(appendR3) + strlen(getFileName(fq3_fn)) +  1 );
        strcpy(fqoutR3, fq_out);
        strcat(fqoutR3, appendR3);
        strcat(fqoutR3, getFileName(fq3_fn));
        seq3 =  kseq_init(fq3);
        if (write_gz){
            o_stream_gz_R3 = gzopen(fqoutR3, "wb2"); // open gz file
            if (!o_stream_gz_R3) {
                file_error(fqoutR3);
            }
        }else{
            o_stream_R3.open(fqoutR3); // output file
            if (!o_stream_R3.is_open()) {
                file_error(fqoutR3);
            }
        }
    }
    
    kseq_t *seq1;
    seq1 =  kseq_init(fq1);
    
    
    const char* appendR1 = "/demux_";
    char *fqoutR1 = (char*)malloc(strlen(fq_out) + strlen(appendR1) + strlen(getFileName(fq1_fn)) + 1);
    strcpy(fqoutR1, fq_out);
    strcat(fqoutR1,appendR1);
    strcat(fqoutR1, getFileName(fq1_fn));
    
    
    if (write_gz) {
        o_stream_gz_R1 = gzopen(fqoutR1, "wb2"); // open gz file
        if (!o_stream_gz_R1) {
            file_error(fqoutR1);
        }
        
    } else {
        o_stream_R1.open(fqoutR1); // output file
        if (!o_stream_R1.is_open()) {
            file_error(fqoutR1);
        }
    }
    
    
    
    // id1_st: bs1: starting position of barcode in read one. -1 if no barcode in read one.
    // id1_len: bl1: length of barcode in read one, if there is no barcode in read one this number is used for trimming beginning of read one.
    // id2_st: bs2: starting position of barcode in read two
    // id2_len: bl2: length of barcode in read two
    // umi_st: us: starting position of UMI
    // umi_len: ul: length of UMI
    
    
    // Define some variables. These are only used if rmlow = T
    int bc1_end, bc2_end; // get total length of index + UMI for read1 and read2
    
    if (id1_st >= 0) // if we have plate index
    {
        bc1_end = id1_st + id1_len;
    }
    else // if no plate information, use id1_len to trim the read 1
    {
        bc1_end = id1_len;
    }
    
    // set barcode end index
    if (umi_st >= 0)
    {
        // This is basically doing the following: bc2_end = max(id2_st + id2_len, umi_st + umi_len)
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
        bc2_end = id2_st + id2_len;
    }
    
    
    
    
    size_t _interrupt_ind = 0;
    // Assuming R1, R2, R3 all are of equal lengths.
    while (((l1 = kseq_read(seq1)) >= 0))
    {
        if (++_interrupt_ind % 4096 == 0) checkUserInterrupt();
        passed_reads++;
        
        
        char * const seq1_name = seq1->name.s;
        char * const seq1_seq = seq1->seq.s;
        int seq1_namelen = seq1->name.l;
        int seq1_seqlen = seq1->seq.l;
        
        char * seq3_name;
        char * seq3_seq;
        int seq3_namelen;
        int seq3_seqlen;
        if (R3){
            if(l3 = kseq_read(seq3) >= 0){
                seq3_name = seq3->name.s;
                seq3_seq = seq3->seq.s;
                seq3_namelen = seq3->name.l;
                seq3_seqlen = seq3->seq.l;
            }
            else{
                Rcpp::Rcout << "read2 file is not of the same length as the barcode fastq file: " << "\n";
            }
        } 
        
        for(int i=0;i<seq2_list.size();i++){
            kseq_t* seq2 = seq2_list[i];
            if (l2 = kseq_read(seq2) >= 0){
                char * const seq2_name = seq2->name.s;
                char * const seq2_seq = seq2->seq.s;
                int seq2_namelen = seq2->name.l;
                int seq2_seqlen = seq2->seq.l;
                
                const int new_name_length1 = seq1_namelen + seq2_seqlen+1;
                seq1->name.s = (char*)realloc(seq1->name.s, new_name_length1); // allocate additional memory
                memmove(seq1_name + seq2_seqlen+1, seq1_name, seq1_namelen * sizeof(char));// move original read name
                memcpy(seq1_name, seq2_seq, seq2_seqlen * sizeof(char)); // copy index one
                seq1_name[seq2_seqlen] = '#'; // add separator
                seq1_name[new_name_length1] = '\0';
                
                if(R3){
                    const int new_name_length2 = seq3_namelen + seq2_seqlen+1;
                    seq3->name.s = (char*)realloc(seq3->name.s, new_name_length2); // allocate additional memory
                    memmove(seq3_name + seq2_seqlen+1, seq3_name, seq3_namelen * sizeof(char));// move original read name
                    memcpy(seq3_name, seq2_seq, seq2_seqlen * sizeof(char)); // copy index one
                    seq3_name[seq2_seqlen] = '#'; // add separator  
                    seq3_name[new_name_length1] = '\0';
                    
                }
            }else{
                Rcpp::Rcout << "read1 file is not the same length as the barcode fastq file: " << "\n";
            }
        }
        
        
        
        // If rmlow parameter is TRUE:
        if(rmlow) { // Only check barcode/UMI quality
            
            // Check quality
            if (!sc_atac_check_qual(seq1->qual.s, bc1_end, min_qual, num_below_min)) {
                removed_low_qual++;
                // Rcout << "R1: "<< std::endl;
                // Rcout << seq1->qual.s << std::endl << std::endl;
                continue;
            } else{
                if(R3){
                    if (!sc_atac_check_qual(seq3->qual.s, bc2_end, min_qual, num_below_min)) {
                        removed_low_qual++;
                        // Rcout << "R3: "<< std::endl;
                        // Rcout << seq3->qual.s << std::endl << std::endl;
                        continue;
                    }
                }
            }
        } // end if(rmlow){
        
        
        // If the rmN parameter is TRUE:
        if(rmN){
            // If find_N is TRUE, then there is an N in the sequence
            if(find_N(seq1)){
                // If there was an N in the sequence, then
                // we add 1 to the counter of reads deleted and 
                // skip the rest of the code in the loop.
                removed_Ns++; // Add 1 to the counter of reads deleted
                continue; // the rest of the lines in the while loop are ignored
            } else{
                if(R3){
                    // If there wasn't an N in the sequence, then we check in R3
                    if(find_N(seq3)){
                        // If there was an N in the sequence, then
                        // we add 1 to the counter of reads deleted and 
                        // skip the rest of the code in the loop.
                        removed_Ns++; // Add 1 to the counter of reads deleted
                        continue; // the rest of the lines in the while loop are ignored
                    }
                } // end if(R3)
            } // end else from if(findN(seq1)) 
        } // end if(rmN)
        
        
        if (write_gz) {
            fq_gz_write(o_stream_gz_R1, seq1, 0); // write to gzipped fastq file
            if(R3){
                fq_gz_write(o_stream_gz_R3, seq3, 0); // write to gzipped fastq file
            }
        } else {
            fq_write(o_stream_R1, seq1, 0); // write to fastq file
            if(R3){
                fq_write(o_stream_R3, seq3, 0); // write to fastq file
            }
        }
        
        
        
    } // end while
    
    // free subStr1
    kseq_destroy(seq1); 
    for(int i=0;i<seq2_list.size();i++){
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
        if(R3){
            gzclose(o_stream_gz_R3);
        }
    }
    Rcpp::Rcout << "Total Reads: " << passed_reads << "\n";
    Rcpp::Rcout << "Total N's removed: " << removed_Ns << "\n";
    Rcpp::Rcout << "removed_low_qual: " << removed_low_qual << "\n";
}














// sc_atac_paired_fastq_to_csv ------------------

void sc_atac_paired_fastq_to_csv(
        char *fq1_fn,
        char *fq3_fn,
        char *fq_out, 
        char *bc_fn, 
        int start,
        int length, 
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
    int passed_reads = 0;
    int removed_Ns = 0;
    int removed_low_qual = 0;
    int exact_match = 0;
    int approx_match = 0;
    
    bool R3 = false;
    bool isUMIR1 = (umi_length > 0 && (strcmp(umi_in,"both") == 0 || strcmp(umi_in,"R1") == 0));
    bool isUMIR2 = (umi_length > 0 && (strcmp(umi_in,"both") == 0 || strcmp(umi_in,"R2") == 0));
    bool write_line;
    
    if(!strcmp(fq3_fn,"")){
        R3 = false;
    }else{
        R3 = true;
    }
    
    int l1 = 0;
    int l2 = 0;
    int l3 = 0;
    gzFile fq1 = gzopen(fq1_fn, "r"); // input fastq
    if (!fq1) {
        file_error(fq1_fn);
    }
    
    
    std::map<std::string, int> barcode_map; 
    std::ifstream bc(bc_fn);
    std::string line;
    if(!bc.is_open()) throw std::runtime_error("Could not open file");
    while(std::getline(bc, line))
    {
        std::stringstream s_stream(line);
        std::string bcode;
        if(s_stream.good()) {
            getline(s_stream, bcode, ','); //get first string delimited by comma
        }
        
        if(s_stream.good()) {
            getline(s_stream, bcode, ','); //get second string delimited by comma
            std::string substr = bcode.substr(0,length);
            barcode_map.insert(std::pair<std::string, int>(substr, 1)); 
        }
        
    }
    
    if(barcode_map.empty()){
        std::stringstream err_msg;
        err_msg << "Error in retrieving barcodes from the barcode File. Please check the barcode file format. " << bc_fn << "\n";
        Rcpp::stop(err_msg.str());
    }
    
    
    
    gzFile fq3;
    gzFile o_stream_gz_R3;
    std::ofstream o_stream_R3;
    gzFile o_stream_gz_R3_Partial;
    std::ofstream o_stream_R3_Partial;
    gzFile o_stream_gz_R3_No;
    std::ofstream o_stream_R3_No;
    
    kseq_t *seq3;
    const char* appendCompleteMatch = "/demultiplexed_completematch_";
    const char* appendPartialMatch = "/demultiplexed_partialmatch_";
    const char* appendNoMatch = "/demultiplexed_nomatch_";
    
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
    
    
    // Define some variables. These are only used if rmlow = T
    int bc1_end, bc2_end; // get total length of index + UMI for read1 and read2
    
    if (id1_st >= 0) // if we have plate index
    {
        bc1_end = id1_st + id1_len;
    }
    else // if no plate information, use id1_len to trim the read 1
    {
        bc1_end = id1_len;
    }
    
    // set barcode end index
    if (umi_start >= 0)
    {
        // This is basically doing the following: bc2_end = max(id2_st + id2_len, umi_start + umi_length)
        if (id2_st + id2_len > umi_start + umi_length)
        {
            bc2_end = id2_st + id2_len;
        }
        else
        {
            bc2_end = umi_start + umi_length;
        }
    }
    else
    {
        bc2_end = id2_st + id2_len;
    }
    
    
    size_t _interrupt_ind = 0;
    // Assuming R1, R2, R3 all are of equal lengths.
    while (((l1 = kseq_read(seq1)) >= 0))
    {
        if (++_interrupt_ind % 50 == 0) checkUserInterrupt();
        passed_reads++;
        
        
        char *seq1_name = seq1->name.s;
        char *seq1_seq = seq1->seq.s;
        int seq1_namelen = seq1->name.l;
        int seq1_seqlen = seq1->seq.l;
        
        char * seq3_name;
        char * seq3_seq;
        int seq3_namelen;
        int seq3_seqlen;
        if (R3){
            if(l3 = kseq_read(seq3) >= 0){
                seq3_name = seq3->name.s;
                seq3_seq = seq3->seq.s;
                seq3_namelen = seq3->name.l;
                seq3_seqlen = seq3->seq.l;
            }
            else{
                Rcpp::Rcout << "R3 file is not of same length as R2: " << "\n";
            }
        } 
        
        
        char* subStr1;
        int bcUMIlen1 = 0;
        if(umi_length > 0 ){
            subStr1 = (char*)malloc(length+umi_length+1);
            memcpy(subStr1, &seq1_seq[start], length );
            subStr1 [length] = '_';
            memcpy(subStr1+length+1, &seq1_seq[umi_start], umi_length);
            subStr1[length + umi_length +1] = '\0';
            bcUMIlen1 = length + umi_length +1;
            
        }else{
            subStr1 = (char*)malloc(length);
            memcpy( subStr1, &seq1_seq[start], length );
            subStr1[length] = '\0';
            bcUMIlen1 =length;
        }
        
        char barcode[length];
        memcpy( barcode, &seq1_seq[start], length );
        barcode[length] = '\0';
        
        
        
        if ( barcode_map.find(barcode) == barcode_map.end() ) {
            for (std::map<std::string,int>::iterator it=barcode_map.begin(); it!=barcode_map.end(); ++it){
                if(hamming_distance(it->first, barcode) <2){
                    approx_match++;
                    const int new_name_length1 = seq1_namelen + bcUMIlen1 + 1;
                    seq1->name.s = (char*)realloc(seq1->name.s, new_name_length1); // allocate additional memory
                    memmove(seq1_name + bcUMIlen1+1, seq1_name, seq1_namelen );// move original read name
                    memcpy(seq1_name, subStr1, bcUMIlen1 * sizeof(char)); // copy index one
                    seq1_name[bcUMIlen1] = '#'; // add separator
                    seq1_name[new_name_length1] = '\0';
                    
                    if ( (umi_start + umi_length) > (start + length) && umi_length > 0){
                        memmove (seq1_seq, seq1_seq + umi_start + umi_length, seq1_seqlen-umi_start-umi_length+1); 
                    }else{
                        memmove (seq1_seq, seq1_seq + start + length, seq1_seqlen-start-length+1); 
                    }
                    
                    
                    // If rmlow parameter is TRUE:
                    if(rmlow) { // Only check barcode/UMI quality
                        // Check quality
                        if (!sc_atac_check_qual(seq1->qual.s, bc1_end, min_qual, num_below_min)) {
                            removed_low_qual++;
                            continue; // the rest of the lines in the while loop are ignored
                        } 
                    } 
                    
                    
                    
                    // If the rmN parameter is TRUE:
                    if(rmN){
                        // If find_N is TRUE, then there is an N in the sequence
                        if(find_N(seq1)){
                            // If there was an N in the sequence, then
                            // we add 1 to the counter of reads deleted and
                            // skip the rest of the code in the loop.
                            removed_Ns++; // Add 1 to the counter of reads deleted
                            continue; // the rest of the lines in the while loop are ignored
                        }
                    } // end if(rmN)
                    
                    
                    if (write_gz) {
                        fq_gz_write(o_stream_gz_R1_Partial, seq1, 0); // write to gzipped fastq file
                    } else {
                        fq_write(o_stream_R1_Partial, seq1, 0); // write to fastq file
                    }
                    
                    
                }else{
                    
                    // If rmlow parameter is TRUE:
                    if(rmlow) { // Only check barcode/UMI quality
                        // Check quality
                        if (!sc_atac_check_qual(seq1->qual.s, bc1_end, min_qual, num_below_min)) {
                            removed_low_qual++;
                            continue; // the rest of the lines in the while loop are ignored
                        } 
                    } 
                    
                    // If the rmN parameter is TRUE:
                    if(rmN){
                        // If find_N is TRUE, then there is an N in the sequence
                        if(find_N(seq1)){
                            // If there was an N in the sequence, then
                            // we add 1 to the counter of reads deleted and
                            // skip the rest of the code in the loop.
                            removed_Ns++; // Add 1 to the counter of reads deleted
                            continue; // the rest of the lines in the while loop are ignored
                        }
                    } // end if(rmN)
                    
                    
                    
                    if (write_gz) {
                        fq_gz_write(o_stream_gz_R1_No, seq1, 0); // write to gzipped fastq file
                    } else {
                        fq_write(o_stream_R1_No, seq1, 0); // write to fastq file
                    }
                    
                }
            }
        } else {
            
            exact_match ++;   
            const int new_name_length1 = seq1_namelen + bcUMIlen1 + 1;
            seq1->name.s = (char*)realloc(seq1->name.s, new_name_length1); // allocate additional memory
            memmove(seq1_name + bcUMIlen1+1, seq1_name, seq1_namelen );// move original read name
            memcpy(seq1_name, subStr1, bcUMIlen1 * sizeof(char)); // copy index one
            seq1_name[bcUMIlen1] = '#'; // add separator
            seq1_name[new_name_length1] = '\0';
            
            if ( (umi_start + umi_length) > (start + length) && umi_length > 0){
                memmove (seq1_seq, seq1_seq + umi_start + umi_length, seq1_seqlen-umi_start-umi_length+1); 
            }else{
                memmove (seq1_seq, seq1_seq + start + length, seq1_seqlen-start-length+1); 
            }
            
            // If rmlow parameter is TRUE:
            if(rmlow) { // Only check barcode/UMI quality
                // Check quality
                if (!sc_atac_check_qual(seq1->qual.s, bc1_end, min_qual, num_below_min)) {
                    removed_low_qual++;
                    continue; // the rest of the lines in the while loop are ignored
                } 
            } 
            
            
            // If the rmN parameter is TRUE:
            if(rmN){
                // If find_N is TRUE, then there is an N in the sequence
                if(find_N(seq1)){
                    // If there was an N in the sequence, then
                    // we add 1 to the counter of reads deleted and
                    // skip the rest of the code in the loop.
                    removed_Ns++; // Add 1 to the counter of reads deleted
                    continue; // the rest of the lines in the while loop are ignored
                }
            } // end if(rmN)
            
            
            if (write_gz) {
                fq_gz_write(o_stream_gz_R1, seq1, 0); // write to gzipped fastq file
            } else {
                fq_write(o_stream_R1, seq1, 0); // write to fastq file
            }
            
        }
        
        
        if(R3){
            char* subStr3;
            int bcUMIlen3 = 0;
            if(umi_length > 0 ){
                subStr3 = (char*)malloc(length+umi_length+1);
                memcpy(subStr3, &seq3_seq[start], length );
                subStr3 [length] = '_';
                memcpy(subStr3+length+1, &seq3_seq[umi_start], umi_length);
                subStr3[length + umi_length +1] = '\0';
                bcUMIlen3 = length + umi_length +1;
                
            }else{
                subStr3 = (char*)malloc(length);
                memcpy( subStr3, &seq3_seq[start], length );
                subStr3[length] = '\0';
                bcUMIlen3 =length;
            }
            
            char barcode[length];
            memcpy( barcode, &seq3_seq[start], length );
            barcode[length] = '\0';
            
            if ( barcode_map.find(subStr3) == barcode_map.end() ) {
                for (std::map<std::string,int>::iterator it=barcode_map.begin(); it!=barcode_map.end(); ++it){
                    if(hamming_distance(it->first, barcode) < 2){
                        const int new_name_length1 = seq3_namelen + bcUMIlen3+1;
                        seq3->name.s = (char*)realloc(seq3->name.s, new_name_length1); // allocate additional memory
                        memmove(seq3_name + bcUMIlen3+1, seq3_name, seq3_namelen * sizeof(char) );// move original read name
                        memcpy(seq3_name, subStr3, bcUMIlen3 * sizeof(char)); // copy index one
                        seq3_name[bcUMIlen3] = '#'; // add separator
                        seq3_name[new_name_length1] = '\0';
                        
                        if ( (umi_start + umi_length) > (start + length) && umi_length > 0){
                            memmove (seq3_seq, seq3_seq + umi_start + umi_length, seq3_seqlen-umi_start-umi_length+1); 
                        }else{
                            memmove (seq3_seq, seq3_seq + start + length, seq3_seqlen-start-length+1); 
                        }
                        
                        // If rmlow parameter is TRUE:
                        if(rmlow) { // Only check barcode/UMI quality
                            // Check quality
                            if (!sc_atac_check_qual(seq3->qual.s, bc2_end, min_qual, num_below_min)) {
                                removed_low_qual++;
                                continue; // the rest of the lines in the while loop are ignored
                            } 
                        } 
                        
                        // If the rmN parameter is TRUE:
                        if(rmN){
                            // If find_N is TRUE, then there is an N in the sequence
                            if(find_N(seq3)){
                                // If there was an N in the sequence, then
                                // we add 1 to the counter of reads deleted and
                                // skip the rest of the code in the loop.
                                removed_Ns++; // Add 1 to the counter of reads deleted
                                continue; // the rest of the lines in the while loop are ignored
                            }
                        } // end if(rmN)
                        
                        
                        if (write_gz) {
                            fq_gz_write(o_stream_gz_R3_Partial, seq3, 0); // write to gzipped fastq file
                        } else {
                            fq_write(o_stream_R3_Partial, seq3, 0); // write to fastq file
                        }
                    }else{
                        
                        // If rmlow parameter is TRUE:
                        if(rmlow) { // Only check barcode/UMI quality
                            // Check quality
                            if (!sc_atac_check_qual(seq3->qual.s, bc2_end, min_qual, num_below_min)) {
                                removed_low_qual++;
                                continue; // the rest of the lines in the while loop are ignored
                            } 
                        }
                        
                        // If the rmN parameter is TRUE:
                        if(rmN){
                            // If find_N is TRUE, then there is an N in the sequence
                            if(find_N(seq3)){
                                // If there was an N in the sequence, then
                                // we add 1 to the counter of reads deleted and
                                // skip the rest of the code in the loop.
                                removed_Ns++; // Add 1 to the counter of reads deleted
                                continue; // the rest of the lines in the while loop are ignored
                            }
                        } // end if(rmN)
                        
                        
                        if (write_gz) {
                            fq_gz_write(o_stream_gz_R3_No, seq3, 0); // write to gzipped fastq file
                        } else {
                            fq_write(o_stream_R3_No, seq3, 0); // write to fastq file
                        }
                    }
                }
            } else {
                const int new_name_length1 = seq3_namelen + bcUMIlen3+1;
                seq3->name.s = (char*)realloc(seq3->name.s, new_name_length1); // allocate additional memory
                memmove(seq3_name + bcUMIlen3+1, seq3_name, seq3_namelen * sizeof(char) );// move original read name
                memcpy(seq3_name, subStr3, bcUMIlen3 * sizeof(char)); // copy index one
                seq3_name[bcUMIlen3] = '#'; // add separator
                seq3_name[new_name_length1] = '\0';
                
                if ( (umi_start + umi_length) > (start + length) && umi_length > 0){
                    memmove (seq3_seq, seq3_seq + umi_start + umi_length, seq3_seqlen-umi_start-umi_length+1); 
                }else{
                    memmove (seq3_seq, seq3_seq + start + length, seq3_seqlen-start-length+1); 
                }
                
                // If the rmN parameter is TRUE:
                if(rmN){
                    // If find_N is TRUE, then there is an N in the sequence
                    if(find_N(seq3)){
                        // If there was an N in the sequence, then
                        // we add 1 to the counter of reads deleted and
                        // skip the rest of the code in the loop.
                        removed_Ns++; // Add 1 to the counter of reads deleted
                        continue; // the rest of the lines in the while loop are ignored
                    }
                } // end if(rmN)
                
                
                if (write_gz) {
                    fq_gz_write(o_stream_gz_R3, seq3, 0); // write to gzipped fastq file
                } else {
                    fq_write(o_stream_R3, seq3, 0); // write to fastq file
                }
                
                
            }
            
            free(subStr3);
            
        }
        free(subStr1);
        
    }
    
    kseq_destroy(seq1); 
    
    bc.close();
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
    }
    Rcpp::Rcout << "Total Reads: " << passed_reads << "\n";
    Rcpp::Rcout << "Exact match Reads: " << exact_match << "\n";
    Rcpp::Rcout << "Approx Match Reads: " << approx_match << "\n";
    Rcpp::Rcout << "Total N's removed: " << removed_Ns << "\n";
    //Rcpp::Rcout << "removed_low_qual: " << removed_low_qual << "\n";
}
















