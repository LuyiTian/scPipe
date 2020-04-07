// parse aligned bam files with gene tags.

#include "parsebam.h"

using std::string;
using std::unordered_map;
using std::ofstream;
using namespace Rcpp;

Bamdemultiplex::Bamdemultiplex(string odir, Barcode b, string cellular_tag, string molecular_tag, string gene_tag, string map_tag, string MT_tag)
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
        cell_mapped_ambiguous[n] = 0;
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
    overall_count_stat["barcode_unmatch_ambiguous_mapping"] = 0;
}

void Bamdemultiplex::write_statistics(string overall_stat_f, string chr_stat_f, string cell_stat_f)
{
    string stat_dir = join_path(out_dir, "stat");
    ofstream overall_stat(join_path(stat_dir, overall_stat_f + ".csv"));
    ofstream chr_stat(join_path(stat_dir, chr_stat_f + ".csv"));
    ofstream cell_stat(join_path(stat_dir, cell_stat_f + ".csv"));
    overall_stat << "status,count" << "\n";

    for (const auto& n : overall_count_stat)
    {
        overall_stat << n.first << "," << n.second << "\n";
    }

    chr_stat << "chromosome name,count" << "\n";
    for (const auto& n : chr_aligned_stat)
    {
        chr_stat << n.first << "," << n.second << "\n";
    }

    cell_stat << "cell_id,unaligned,aligned_unmapped,mapped_to_exon,mapped_to_intron,ambiguous_mapping,mapped_to_ERCC,mapped_to_MT" << "\n";
    for (const auto& n : bar.cellid_list)
    {
        cell_stat << n << "," << cell_unaligned[n] << "," << \
            cell_align_unmapped[n] << "," << \
            cell_mapped_exon[n] << "," << \
            cell_mapped_intron[n] << "," << \
            cell_mapped_ambiguous[n] << "," << \
            cell_ERCC[n] << ","<< \
            cell_MT[n] << "\n";
    }
}

int Bamdemultiplex::clean_bam_barcode(string bam_path, string out_bam, int max_mismatch, int nthreads)
{
    check_file_exists(bam_path); // htslib does not check if file exist so we do it manually
    bam1_t *b = bam_init1();
    BGZF *fp = bgzf_open(bam_path.c_str(), "r");
    bam_hdr_t *header = bam_hdr_read(fp);
    
    int hts_retcode;

    samFile *of = sam_open(out_bam.c_str(), "wb"); // output file
    hts_retcode = sam_hdr_write(of, header); // (void) explicitly discard return value

    // Early benchmarking shows BAM reading doesn't saturate even 2 cores
    // so capped reading threads to 2
    int out_threads = std::min(nthreads - 1, 2);
    // make sure we're not using 0 or negative cores
    out_threads = std::max(nthreads, 1);

    hts_tpool *p;
    const int queue_size = 64;

    if (out_threads > 1) {
        p = hts_tpool_init(out_threads);
        bgzf_thread_pool(fp, p, queue_size);
    }

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
        Rcpp::Rcout << "Warning: mitochondrial chromosome not found using chromosome name `"<< mt_tag << "`.\n";
    }

    const char * c_ptr = c_tag.c_str();

    string bc_seq;
    string match_res;
    int map_status;

    size_t _interrupt_ind = 0;

    while (bam_read1(fp, b) >= 0)
    {
        if (++_interrupt_ind % 1024 == 0) checkUserInterrupt();
        //match barcode
        bc_seq = (char*)(bam_aux_get(b, c_ptr) + 1); // +1 to skip `Z`
        match_res = bar.get_closest_match(bc_seq, max_mismatch);

        bool is_unmapped = (b->core.flag & BAM_FUNMAP) > 0;
        // if the read is aligned and with matched barcode.
        if ((!is_unmapped) & (!match_res.empty())) 
        {
            bam_aux_update_str(b, c_ptr, match_res.length()+1, match_res.c_str());
            hts_retcode = sam_write1(of, header, b); // (void) discards return value
        }
    }

    sam_close(of);
    bgzf_close(fp);
    return 0;
}

int Bamdemultiplex::barcode_demultiplex(string bam_path, int max_mismatch, bool has_UMI, int nthreads)
{
    check_file_exists(bam_path); // htslib does not check if file exist so we do it manually
    bam1_t *b = bam_init1();
    BGZF *fp = bgzf_open(bam_path.c_str(), "r");
    bam_hdr_t *header = bam_hdr_read(fp);

    // Early benchmarking shows BAM reading doesn't saturate even 2 cores
    // so capped reading threads to 2
    int out_threads = std::min(nthreads - 1, 2);
    // make sure we're not using 0 or negative cores
    out_threads = std::max(nthreads, 1);

    hts_tpool *p;
    const int queue_size = 64;

    if (out_threads > 1) {
        p = hts_tpool_init(out_threads);
        bgzf_thread_pool(fp, p, queue_size);
    }

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
        Rcpp::Rcout << "Warning: mitochondrial chromosome not found using chromosome name `"<< mt_tag << "`.\n";
    }

    string output_dir = join_path(out_dir, "count");
    unordered_map<string, string> out_fn_path = bar.get_count_file_path(output_dir);
    unordered_map<string, std::vector<string>> out_reads;
    const char * c_ptr = c_tag.c_str();
    const char * m_ptr = m_tag.c_str();
    const char * g_ptr = g_tag.c_str();
    const char * a_ptr = a_tag.c_str();

    string bc_seq;
    string match_res;
    int map_status;

    size_t _interrupt_ind = 0;

    while (bam_read1(fp, b) >= 0)
    {
        if (++_interrupt_ind % 1024 == 0) checkUserInterrupt();
        //match barcode
        bc_seq = (char*)(bam_aux_get(b, c_ptr) + 1); // +1 to skip `Z`
        match_res = bar.get_closest_match(bc_seq, max_mismatch);

        bool is_unmapped = (b->core.flag & BAM_FUNMAP) > 0;
        if (is_unmapped)
        {
            if (match_res.empty())
            {
                overall_count_stat["barcode_unmatch_unaligned"]++;
            }
            else
            {
                overall_count_stat["barcode_match"]++;
                cell_unaligned[bar.barcode_dict[match_res]]++;
            }

        }
        else
        {
            map_status = bam_aux2i(bam_aux_get(b, a_ptr));
            chr_aligned_stat[header->target_name[b->core.tid]]++;
            if (bam_aux_get(b, g_ptr)) // found a gene; read mapped to transcriptome
            {

                if (match_res.empty())
                {
                    overall_count_stat["barcode_unmatch_mapped_to_exon"]++;
                }
                else
                {
                    overall_count_stat["barcode_match"]++;
                    if (std::strncmp (header->target_name[b->core.tid],"ERCC",4) == 0)
                    {
                        cell_ERCC[bar.barcode_dict[match_res]]++;
                    }
                    else
                    {
                        cell_mapped_exon[bar.barcode_dict[match_res]]++;
                    }
                    if (b->core.tid == mt_idx)
                    {
                        cell_MT[bar.barcode_dict[match_res]]++;
                    }
                    if (a_tag.empty())
                    {
                        if(has_UMI)
                        {
                            out_reads[bar.barcode_dict[match_res]].push_back(string(bam_aux2Z(bam_aux_get(b, g_ptr)))+","+
                            string(bam_aux2Z(bam_aux_get(b, m_ptr)))+","+
                            std::to_string(b->core.pos));
                        }
                        else
                        {
                            out_reads[bar.barcode_dict[match_res]].push_back(string(bam_aux2Z(bam_aux_get(b, g_ptr)))+","+
                            string(bam_get_qname(b))+","+
                            std::to_string(b->core.pos));
                        }
                        
                    }
                    else
                    {
                        if(has_UMI)
                        {
                            out_reads[bar.barcode_dict[match_res]].push_back(string(bam_aux2Z(bam_aux_get(b, g_ptr)))+","+
                            string(bam_aux2Z(bam_aux_get(b, m_ptr)))+","+
                            std::to_string(-map_status));
                        }
                        else
                        {
                            out_reads[bar.barcode_dict[match_res]].push_back(string(bam_aux2Z(bam_aux_get(b, g_ptr)))+","+
                            string(bam_get_qname(b))+","+
                            std::to_string(-map_status));
                        }
                    }

                }
            }
            else
            {
                // return:
                //  <=0 - unique map to exon, number indicate the distance to transcript end pos
                //  1 - ambiguous map to multiple exon
                //  2 - map to intron
                //  3 - unmapped
                //  4 - unaligned
                if (match_res.empty())
                {
                    if (a_tag.empty())
                    {
                        overall_count_stat["barcode_unmatch_aligned"]++;
                    }
                    else
                    {
                        if (map_status == 1)
                        {
                            overall_count_stat["barcode_unmatch_ambiguous_mapping"]++;
                        }
                        else if (map_status == 2)
                        {
                            overall_count_stat["barcode_unmatch_mapped_to_intron"]++;
                        }
                        else
                        {
                            overall_count_stat["barcode_unmatch_aligned"]++;
                        }
                    }

                }
                else
                {
                    overall_count_stat["barcode_match"]++;
                    if (a_tag.empty())
                    {
                        cell_align_unmapped[bar.barcode_dict[match_res]]++;
                    }
                    else
                    {
                        if (map_status == 1)
                        {
                            cell_mapped_ambiguous[bar.barcode_dict[match_res]]++;
                        }
                        else if (map_status == 2)
                        {
                            cell_mapped_intron[bar.barcode_dict[match_res]]++;
                        }
                        else
                        {
                            cell_align_unmapped[bar.barcode_dict[match_res]]++;
                        }
                    }
                }
            }
        }
    }

    for (auto const& fn: out_fn_path)
    {
        ofstream ofile(fn.second);
        ofile << "gene_id,UMI,position\n";
        for (auto const& rd: out_reads[fn.first])
        {
            ofile << rd << "\n";
        }
        ofile.close();
    }



    bgzf_close(fp);
    return 0;
}
