// parse aligned bam files with gene tags.

#include "parsebam.h"

Bamdemultiplex::Bamdemultiplex(std::string odir, Barcode b, std::string cellular_tag, std::string molecular_tag, std::string gene_tag, std::string map_tag, std::string MT_tag)
{
    bar = b;
    c_tag = cellular_tag;
    m_tag = molecular_tag;
    g_tag = gene_tag;
    a_tag = map_tag;
    out_dir = odir;
    mt_tag = MT_tag;

    for (const auto& n : bar.cellid_list)
    {
        cell_mapped_exon[n] = 0;
        cell_mapped_intron[n] = 0;
        cell_mapped_ambigious[n] = 0;
        cell_align_unmapped[n] = 0;
        cell_unaligned[n] = 0;
        cell_ERCC[n] = 0;
        cell_MT[n] = 0;
    }

    overall_count_stat["barcode_match"] = 0;
    overall_count_stat["barcode_unmatch_unaligned"] = 0;
    overall_count_stat["barcode_unmatch_aligned"] = 0;
    overall_count_stat["barcode_unmatch_mapped_to_exon"] = 0;
    overall_count_stat["barcode_unmatch_mapped_to_intron"] = 0;
    overall_count_stat["barcode_unmatch_ambigious_mapping"] = 0;
}

Bamdemultiplex::~Bamdemultiplex()
{
    // do sth
}


int Bamdemultiplex::write_statistics(std::string overall_stat_f, std::string chr_stat_f, std::string cell_stat_f)
{
    std::string stat_dir = join_path(out_dir, "stat");
    std::ofstream overall_stat(join_path(stat_dir, overall_stat_f+".csv"));
    std::ofstream chr_stat(join_path(stat_dir, chr_stat_f+".csv"));
    std::ofstream cell_stat(join_path(stat_dir, cell_stat_f+".csv"));
    overall_stat << "status,count" << std::endl;
    for (const auto& n : overall_count_stat)
    {
        overall_stat << n.first << "," << n.second << std::endl;
    }

    chr_stat << "chromosome name,count" << std::endl;
    for (const auto& n : chr_aligned_stat)
    {
        chr_stat << n.first << "," << n.second << std::endl;
    }

    cell_stat << "cell_id,unaligned,aligned_unmapped,mapped_to_exon,mapped_to_intron,ambigious_mapping,mapped_to_ERCC,mapped_to_MT" << std::endl;
    for (const auto& n : bar.cellid_list)
    {
        cell_stat << n << "," << cell_unaligned[n] << "," << \
            cell_align_unmapped[n] << "," << \
            cell_mapped_exon[n] << "," << \
            cell_mapped_intron[n] << "," << \
            cell_mapped_ambigious[n] << "," << \
            cell_ERCC[n] << ","<< \
            cell_MT[n] << std::endl;
    }
}

int Bamdemultiplex::barcode_demultiplex(std::string bam_path, int max_mismatch)
{
    check_file_exists(bam_path); // htslib does not check if file exist so we do it manually
    bam1_t *b = bam_init1();
    BGZF *fp = bgzf_open(bam_path.c_str(), "r");
    bam_hdr_t *header = bam_hdr_read(fp);
    int mt_idx = -1;
    for (int i = 0; i < header->n_targets; ++i)
    {
        chr_aligned_stat[header->target_name[i]] = 0;
        if (strcmp(header->target_name[i], mt_tag.c_str()) == 0)
        {
            mt_idx = i;
        }
    }
    if (mt_idx == -1)
    {
        std::cout << "Warning: mitachondral chromosome not found using chromosome name `"<< mt_tag << "`.\n";
    }
    std::unordered_map<std::string, std::ofstream> outfn_dict = bar.get_count_file_w(join_path(out_dir, "count"));
    const char * c_ptr = c_tag.c_str();
    const char * m_ptr = m_tag.c_str();
    const char * g_ptr = g_tag.c_str();
    const char * a_ptr = a_tag.c_str();
    std::string bc_seq;
    std::string match_res;
    int tmp_cnt = 0;
    int ret_code = 999999; // see `transcriptmapping.h` for different code
    while (bam_read1(fp, b) >= 0)
    {
        tmp_cnt++;
        //match barcode
        bc_seq = (char*)(bam_aux_get(b, c_ptr)+1); // +1 to skip `Z`
        match_res = bar.get_closest_match(bc_seq, max_mismatch);
        if ((b->core.flag&BAM_FUNMAP) > 0)
        {
            if (match_res.empty())
            {
                overall_count_stat["barcode_unmatch_unaligned"] ++;
            }
            else
            {
                overall_count_stat["barcode_match"] ++;
                cell_unaligned[bar.barcode_dict[match_res]] ++;
            }

        }
        else
        {
            chr_aligned_stat[header->target_name[b->core.tid]] ++;
            if (bam_aux_get(b, g_ptr)) // found a gene; read mapped to transcriptome
            {

                if (match_res.empty())
                {
                    overall_count_stat["barcode_unmatch_mapped_to_exon"] ++;
                }
                else
                {
                    overall_count_stat["barcode_match"] ++;
                    if (std::strncmp (header->target_name[b->core.tid],"ERCC",4) == 0)
                    {
                        cell_ERCC[bar.barcode_dict[match_res]] ++;
                    }
                    else
                    {
                        cell_mapped_exon[bar.barcode_dict[match_res]] ++;
                    }
                    if (b->core.tid == mt_idx)
                    {
                        cell_MT[bar.barcode_dict[match_res]] ++;
                    }
                    if (a_tag.empty())
                    {
                        outfn_dict[bar.barcode_dict[match_res]] <<\
                         (bam_aux_get(b, g_ptr)+1) << "," <<\
                          (bam_aux_get(b, m_ptr)+1) << "," <<\
                           b->core.pos << std::endl;                        
                    }
                    else
                    {
                        outfn_dict[bar.barcode_dict[match_res]] <<\
                         (bam_aux_get(b, g_ptr)+1) << "," <<\
                          (bam_aux_get(b, m_ptr)+1) << "," <<\
                           (-std::atoi((char*)bam_aux_get(b, a_ptr)+1)) << std::endl;                          
                    }

                    
                }
            }
            else
            {
                // return:
                //  <=0 - unique map to exon, number indicate the distance to transcript end pos
                //  1 - ambigious map to multiple exon
                //  2 - map to intron
                //  3 - unmapped
                //  4 - unaligned
                if (match_res.empty())
                {
                    if (a_tag.empty())
                    {
                        overall_count_stat["barcode_unmatch_aligned"] ++;
                    }
                    else
                    {
                        tmp_cnt = std::atoi((char*)bam_aux_get(b, a_ptr)+1);
                        if (tmp_cnt == 1)
                        {
                            overall_count_stat["barcode_unmatch_ambigious_mapping"] ++;
                        }
                        else if (tmp_cnt == 2)
                        {
                            overall_count_stat["barcode_unmatch_mapped_to_intron"] ++;
                        }
                        else
                        {
                            overall_count_stat["barcode_unmatch_aligned"] ++;
                        }
                    }
                    
                }
                else
                {
                    overall_count_stat["barcode_match"] ++;
                    if (a_tag.empty())
                    {
                        cell_align_unmapped[bar.barcode_dict[match_res]] ++;                    
                    }
                    else
                    {
                        tmp_cnt = std::atoi((char*)bam_aux_get(b, a_ptr)+1);
                        if (tmp_cnt == 1)
                        {
                            cell_mapped_ambigious[bar.barcode_dict[match_res]] ++; 
                        }
                        else if (tmp_cnt == 2)
                        {
                            cell_mapped_intron[bar.barcode_dict[match_res]] ++; 
                        }
                        else
                        {
                            cell_align_unmapped[bar.barcode_dict[match_res]] ++;  
                        }
                    }

                }
                
            }
        }

    }
    bgzf_close(fp);
    return 0;
}

/*
int main(int argc, char const *argv[])
{
    bam1_t *q = bam_init1();
    BGZF *fp = bgzf_open("/Users/luyi/Downloads/star_gene_exon_tagged_clean.bam","r");
    bam_hdr_t *header = bam_hdr_read(fp);
    std::cout << header->target_name[20] << std::endl;
    int ii = 0;
    while (bam_read1(fp, q) >= 0 && ii <100){
        printf("%d\t", q->core.flag&BAM_FMUNMAP);
        printf("%s\t", bam_get_qname(q));
        printf("%s\t", header->target_name[q->core.tid]);
        printf("%d\t", q->core.pos);
        if (bam_aux_get(q, "GE")){
            printf(":::%s:::\n", bam_aux_get(q, "GE"));
            ii++;
        }
    }
    //printf("%s\n\n",  header->target_name[20]);
    return 0;
}
*/