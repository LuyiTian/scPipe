//trim_barcode
#include "trimbarcode.h"

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
