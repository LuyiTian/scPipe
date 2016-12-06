//trim_barcode
#include "trimbarcode.h"

/*
void dump_read(bam1_t* b) {
    printf("->core.tid:(%d)\n", b->core.tid);
    printf("->core.pos:(%d)\n", b->core.pos);
    printf("->core.bin:(%d)\n", b->core.bin);
    printf("->core.qual:(%d)\n", b->core.qual);
    printf("->core.l_qname:(%d)\n", b->core.l_qname);
    printf("->core.flag:(%d)\n", b->core.flag);
    printf("->core.n_cigar:(%d)\n", b->core.n_cigar);
    printf("->core.l_qseq:(%d)\n", b->core.l_qseq);
    printf("->core.mtid:(%d)\n", b->core.mtid);
    printf("->core.mpos:(%d)\n", b->core.mpos);
    printf("->core.isize:(%d)\n", b->core.isize);
    if (b->data) {
        printf("->data:");
        int i;
        for (i = 0; i < b->l_data; ++i) {
            printf("%x ", b->data[i]);
        }
        printf("\n");
    }
    if (b->core.l_qname) {
        printf("qname: %s\n",bam_get_qname(b));
    }
    if (b->core.l_qseq) {
        printf("qseq:");
        int i;
        for (i = 0; i < b->core.l_qseq; ++i) {
            printf("%c", seq_nt16_str[bam_seqi(bam_get_seq(b),i)]);
        }
        printf("\n");
        printf("qual:");
        for (i = 0; i < b->core.l_qseq; ++i) {
            printf("%c",bam_get_qual(b)[i]);
        }
        printf("\n");

    }

    if (bam_get_l_aux(b)) {
        int i = 0;
        uint8_t* aux = bam_get_aux(b);

        while (i < bam_get_l_aux(b)) {
            printf("%.2s:%c:",aux+i,*(aux+i+2));
            i += 2;
            switch (*(aux+i)) {
                case 'Z':
                    while (*(aux+1+i) != '\0') { putc(*(aux+1+i), stdout); ++i; }
                    break;
            }
            putc('\n',stdout);
            ++i;++i;
        }
    }
    printf("\n");
}
*/


void kseq_t_to_bam_t(kseq_t *seq, bam1_t *b, int trim_n)
{
    int seq_l = seq->seq.l-trim_n; // seq length after trim the barcode
    b->l_data = seq->name.l+1+(int)(1.5*seq_l+(seq_l % 2 != 0)); // +1 includes the tailing '\0'
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
    b->core.l_qname = seq->name.l+1; // +1 includes the tailing '\0'
    b->core.l_qseq = seq_l;
    b->core.n_cigar = 0; // we have no cigar sequence
    memcpy(b->data, seq->name.s, sizeof(char)*seq->name.l); // first set qname
    b->data[seq->name.l] = '\0';
    uint8_t *s = bam_get_seq(b);
    int i = 0;
    for (i = 0; i < b->core.l_qseq; ++i) // set sequence
    {

        bam1_seq_seti(s, i, seq_nt16_table[seq->seq.s[i+trim_n]]);
    }
    s = bam_get_qual(b);
    for (i = 0; i < b->core.l_qseq; ++i) // set quality
    {
        s[i] = seq->qual.s[i+trim_n]-33;
    }
}


bool check_qual(char *qual_s, int trim_n, int thr, int below_thr)
{
    int not_pass = 0;
    for (int i = 0; i < trim_n; i++)
    {
        if ((int)qual_s[i] <= thr){
            not_pass++;
        }
    }
    return not_pass>below_thr?false:true;
}


bool N_check(char *seq, int trim_n){
    bool pass = true;
    char *ptr = strchr(seq, 'N');
    if (ptr)
    {
        int index = ptr - seq;
        if (index <= trim_n)
        {
            pass = false;
        }
    }
    return pass;
}

// most aligners do not support optional fields in bam files.
/*
void paired_fastq_to_bam(char *fq1_fn, char *fq2_fn, char *bam_out, read_s r, filter_s fl, std::string cellular_tag, std::string molecular_tag)
{
    const char * c_ptr = cellular_tag.c_str();
    const char * m_ptr = molecular_tag.c_str();

    int passed_reads = 0;
    int removed_have_N = 0;
    int removed_low_qual = 0;

    int l1 = 0;
    int l2 = 0;
    gzFile *fq1 = (gzFile *)gzopen(fq1_fn, "r"); // input fastq
    if (!fq1){fprintf(stderr, "cant open file: %s\n", fq1_fn); exit(EXIT_FAILURE);}
    gzFile *fq2 = (gzFile *)gzopen(fq2_fn, "r");
    if (!fq2){fprintf(stderr, "cant open file: %s\n", fq2_fn); exit(EXIT_FAILURE);}

    samFile *fp = sam_open(bam_out,"wb"); // output file

    // write header
    bam_hdr_t *hdr = bam_hdr_init();
    hdr->l_text = strlen(empty_header);
    hdr->text = strdup(empty_header);
    hdr->n_targets = 0;
    sam_hdr_write(fp, hdr);

    uint8_t *idx = new uint8_t[r.id1_len+r.id2_len+1];
    uint8_t *umi = new uint8_t[r.umi_len+1];

    int bc1_end, bc2_end; // get total length of index+UMI for read1 and read2
    if (r.id1_st>=0) // if we have plate index
    {
        bc1_end = r.id1_st+r.id1_len;
        idx[r.id1_len+r.id2_len] = '\0'; // add null terminater
    }
    else // if dont have plate, use r.id1_len to trim the read 1
    {
        bc1_end = r.id1_len;
        uint8_t *idx = new uint8_t[r.id2_len+1];
        idx[r.id2_len] = '\0';
    }
    if (r.umi_st >= 0)
    {
        if (r.id2_st+r.id2_len > r.umi_st+r.umi_len)
        {
            bc2_end = r.id2_st+r.id2_len;
        }
        else
        {
            bc2_end = r.umi_st+r.umi_len;
        }
        umi[r.umi_len] = '\0';
    }
    else
    {
        bc2_end = r.id2_st+r.id2_len;
    }

    kseq_t *seq1;
    seq1 =  kseq_init(fq1);
    kseq_t *seq2;
    seq2 =  kseq_init(fq2);
    // main loop, iter through each fastq records
    // ideally there should be equal number of reads in fq1 and fq2. we dont check this.
    while (((l1 = kseq_read(seq1)) >= 0) && ((l2 = kseq_read(seq2)) >= 0))
    {  
        // qual check before we do anything
        if (fl.if_check_qual)
        {
            if(!(check_qual(seq1->seq.s, bc1_end, fl.min_qual, fl.num_below_min) && check_qual(seq2->seq.s, bc2_end, fl.min_qual, fl.num_below_min)))
            {
                removed_low_qual ++;
                continue;
            }
        }
        if (fl.if_remove_N)
        {
            if(!(N_check(seq1->seq.s, bc1_end) && N_check(seq2->seq.s, bc2_end)))
            {
                removed_have_N ++;
                continue;
            }
        }
        
        passed_reads++;

        bam1_t *b = bam_init1();

        kseq_t_to_bam_t(seq1, b, bc1_end);

        if (r.id1_st>=0) // with plate index (index in read 1)
        {
            memcpy(idx, seq1->seq.s+r.id1_st, r.id1_len*sizeof(uint8_t));
            memcpy(idx+r.id1_len, seq2->seq.s+r.id2_st, r.id2_len*sizeof(uint8_t));
            bam_aux_append(b, c_ptr, 'Z', r.id1_len+r.id2_len+1, (uint8_t*)idx);
        }
        else // without plate index
        {
            memcpy(idx, seq2->seq.s+r.id2_st, r.id2_len*sizeof(uint8_t));
            bam_aux_append(b, c_ptr, 'Z', r.id2_len+1, (uint8_t*)idx);
        }

        if (r.umi_st >= 0) // if we have UMI
        {
            memcpy(umi, seq2->seq.s+r.umi_st, r.umi_len*sizeof(uint8_t));
            bam_aux_append(b, m_ptr, 'Z', r.umi_len+1, umi);
        }

        int ret = sam_write1(fp, hdr, b);
        if (ret < 0)
        {
            std::cout << "fail to write the bam file: " << seq1->name.s << std::endl;
            std::cout << "return code: " << ret << std::endl;
            exit(EXIT_FAILURE);
        }
        dump_read(b);
        bam_destroy1(b);
    }

    std::cout << passed_reads << "\t" << removed_low_qual << "\t" << removed_have_N << std::endl;

    delete idx;
    delete umi;
    kseq_destroy(seq1); kseq_destroy(seq2); // free seq 
    gzclose(fq1); gzclose(fq2); // close fastq file
    sam_close(fp); // close bam file

}
*/



void paired_fastq_to_bam(char *fq1_fn, char *fq2_fn, char *bam_out, read_s r, filter_s fl)
{

    int passed_reads = 0;
    int removed_have_N = 0;
    int removed_low_qual = 0;

    int l1 = 0;
    int l2 = 0;
    gzFile fq1 = gzopen(fq1_fn, "r"); // input fastq
    if (!fq1){fprintf(stderr, "cant open file: %s\n", fq1_fn); exit(EXIT_FAILURE);}
    gzFile fq2 = gzopen(fq2_fn, "r");
    if (!fq2){fprintf(stderr, "cant open file: %s\n", fq2_fn); exit(EXIT_FAILURE);}

    samFile *fp = sam_open(bam_out,"wb"); // output file

    // write header
    bam_hdr_t *hdr = bam_hdr_init();
    hdr->l_text = strlen(empty_header);
    hdr->text = strdup(empty_header);
    hdr->n_targets = 0;
    sam_hdr_write(fp, hdr);

    uint8_t *se = new uint8_t[r.id1_len+r.id2_len+1];
    uint8_t *umi = new uint8_t[r.umi_len+1];

    int bc1_end, bc2_end; // get total length of index+UMI for read1 and read2
    int state; // 0 for two index with umi, 1 for two index without umi, 2 for one index with umi, 3 for one index without umi
    if (r.id1_st>=0) // if we have plate index
    {
        state = 0;
        bc1_end = r.id1_st+r.id1_len;
    }
    else // if dont have plate, use r.id1_len to trim the read 1
    {
        state = 2;
        bc1_end = r.id1_len;
        uint8_t *idx = new uint8_t[r.id2_len+1];
    }
    if (r.umi_st >= 0)
    {
        if (r.id2_st+r.id2_len > r.umi_st+r.umi_len)
        {
            bc2_end = r.id2_st+r.id2_len;
        }
        else
        {
            bc2_end = r.umi_st+r.umi_len;
        }
    }
    else
    {
        state++; // no umi
        bc2_end = r.id2_st+r.id2_len;
    }

    int name_offset;
    if (state == 0)
    {
        name_offset = r.id1_len + r.id2_len + r.umi_len + 2;
    }
    else if (state == 1)
    {
        name_offset = r.id1_len + r.id2_len  + 2;
    }
    else if (state == 2)
    {
        name_offset = r.id2_len + r.umi_len + 2;
    }
    else if (state == 3)
    {
        name_offset = r.id2_len + 2;
    }



    kseq_t *seq1;
    seq1 =  kseq_init(fq1);
    kseq_t *seq2;
    seq2 =  kseq_init(fq2);
    // main loop, iter through each fastq records
    // ideally there should be equal number of reads in fq1 and fq2. we dont check this.
    while (((l1 = kseq_read(seq1)) >= 0) && ((l2 = kseq_read(seq2)) >= 0))
    {  
        // qual check before we do anything
        if (fl.if_check_qual)
        {
            if(!(check_qual(seq1->seq.s, bc1_end, fl.min_qual, fl.num_below_min) && check_qual(seq2->seq.s, bc2_end, fl.min_qual, fl.num_below_min)))
            {
                removed_low_qual ++;
                continue;
            }
        }
        if (fl.if_remove_N)
        {
            if(!(N_check(seq1->seq.s, bc1_end) && N_check(seq2->seq.s, bc2_end)))
            {
                removed_have_N ++;
                continue;
            }
        }
        
        passed_reads++;

        bam1_t *b = bam_init1();

        seq1->name.s = (char*)realloc(seq1->name.s, name_offset + seq1->name.l);
        memcpy(seq1->name.s+name_offset, seq1->name.s, seq1->name.l*sizeof(char)); // move original read name
        if (state == 0)
        {
            memcpy(seq1->name.s, seq1->seq.s+r.id1_st, r.id1_len*sizeof(char)); // copy index one
            memcpy(seq1->name.s+r.id1_len, seq2->seq.s+r.id2_st, r.id2_len*sizeof(char)); // copy index two
            seq1->name.s[r.id1_len+r.id2_len] = '_'; // add separater
            memcpy(seq1->name.s+r.id1_len+r.id2_len+1, seq2->seq.s+r.umi_st, r.umi_len*sizeof(char)); // copy umi

        }
        else if (state == 1)
        {
            memcpy(seq1->name.s, seq1->seq.s+r.id1_st, r.id1_len*sizeof(char)); // copy index one
            memcpy(seq1->name.s+r.id1_len, seq2->seq.s+r.id2_st, r.id2_len*sizeof(char)); // copy index two
            seq1->name.s[r.id1_len+r.id2_len] = '_'; // add separater
        }
        else if (state == 2)
        {
            memcpy(seq1->name.s, seq2->seq.s+r.id2_st, r.id2_len*sizeof(char)); // copy index two
            seq1->name.s[r.id2_len] = '_'; // add separater
            memcpy(seq1->name.s+r.id2_len+1, seq2->seq.s+r.umi_st, r.umi_len*sizeof(char)); // copy umi
        }
        else if (state == 3)
        {
            memcpy(seq1->name.s, seq2->seq.s+r.id2_st, r.id2_len*sizeof(char)); // copy index two
            seq1->name.s[r.id2_len] = '_'; // add separater
        }
        seq1->name.s[name_offset-1] = '#';
        seq1->name.l = name_offset + seq1->name.l;

        kseq_t_to_bam_t(seq1, b, bc1_end);
        int ret = sam_write1(fp, hdr, b);
        if (ret < 0)
        {
            std::cout << "fail to write the bam file: " << seq1->name.s << std::endl;
            std::cout << "return code: " << ret << std::endl;
            exit(EXIT_FAILURE);
        }
        //dump_read(b);
        bam_destroy1(b);
    }

    //std::cout << passed_reads << "\t" << removed_low_qual << "\t" << removed_have_N << std::endl;

    kseq_destroy(seq1); kseq_destroy(seq2); // free seq 
    gzclose(fq1); gzclose(fq2); // close fastq file
    sam_close(fp); // close bam file
    std::cout << "pass QC: " << passed_reads << std::endl;
    std::cout << "removed_have_N: " << removed_have_N << std::endl;
    std::cout << "removed_low_qual: " << removed_low_qual << std::endl;
}