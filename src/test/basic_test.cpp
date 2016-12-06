// test_utils.cpp
#include <string>
#include <iostream>
#include "../utils.h"
#include "../cellbarcode.h"
#include "../parsebam.h"
#include "../parsecount.h"
#include "../trimbarcode.h"
#include "../transcriptmapping.h"


void print_result(bool is_good, int &passed_test, int &failed_test){
    if(is_good)
    {
        std::cout << " passed" << std::endl;
        passed_test++;
    }
    else
    {
        std::cout << " failed" << std::endl;
        failed_test++;
    }
}


bool test_join_path()
{
    std::cout << "\ttest join_path() ...";
    std::string a = "aa/bb/cc";
    std::string b = "aa/bb/cc/";
    std::string c = "mm.csv";
    return join_path(a,c) == join_path(b,c);
}


bool test_hamming_distance()
{
    std::cout << "\ttest hamming_distance() ...";
    std::string a = "ATCGTAAC";
    std::string b = "ATGCTAAC";
    int dist = hamming_distance(a,b);
    return dist==2;
}


bool test_vector_counter()
{
    std::cout << "\ttest hamming_distance() ...";
    std::vector<std::string> v = {"ATGCTAAC", "ATCTGCCC", "ATGCTAAC", "GTAGTAG"};
    std::unordered_map<std::string, int> tmp_res = vector_counter(v);
    return tmp_res["ATGCTAAC"] == 2;
}


bool test_overlap()
{
    std::cout << "\ttest Interval::overlap() ...";
    Interval a = Interval(500,700);
    Interval b = Interval(300,480);
    bool tmp = a.overlap(400,510) == 0 && \
        a.overlap(690,790) == 0 && \
        a.overlap(300,900) == 0  && \
        a > b;
    return tmp;
}


bool test_read_anno()
{
    std::cout << "\ttest Barcode::read_anno() ...";
    Barcode bar;
    std::string fn = "test/test_data/barcode_anno.csv";
    bar.read_anno(fn);
    bool tmp = bar.barcode_dict.count("ATCTGCCC")>0 && \
        bar.barcode_dict["ACGATCGA"] == "CELL003" && \
        bar.cellid_list.size() == 3 && \
        bar.barcode_list.size() == 4;
    return tmp;
}


bool test_get_closest_match()
{
    std::cout << "\ttest Barcode::get_closest_match() ...";
    Barcode bar;
    std::string fn = "test/test_data/barcode_anno.csv";
    bar.read_anno(fn);
    std::string bc1 = "ATGATAAA";
    std::string bc2 = "AAAAAAAA";
    bool tmp = bar.get_closest_match(bc1, 2).compare("ATGATAAT") == 0 && \
        bar.get_closest_match(bc2, 2).empty();
    return tmp;

}


bool test_barcode_demultiplex1()
{
    std::cout << "\ttest Bamdemultiplex::barcode_demultiplex() check_file_exists()...";
    //std::string bamfn = "/home/users/allstaff/tian.l/public_data/GSM1544799/GSM1544799_SpeciesMix_HundredSTAMPs.bam";
    //std::string annofn = "/home/users/allstaff/tian.l/public_data/GSM1544799/GSM1544799_smp_list.csv";
    std::string tmp_out = "/Users/luyi/Downloads/dropseq";
    std::string bamfn = "/Users/luyi/Downloads/star_gene_exon_tagged_clean.bam";
    std::string annofn = "/Users/luyi/Downloads/anno_tmp.csv";
    std::string bc = "XC";
    std::string mb = "XM";
    std::string gb = "GE";
    std::string mt = "MT";
    std::string am;
    int max_mismatch = 1;
    bool tmp;
    Barcode bar;
    bar.read_anno(annofn);
    try 
    {
        std::string wrong_fn = bamfn + "c1veesvwr"; // this is a bad file path
        Bamdemultiplex bam_de = Bamdemultiplex(tmp_out, bar, bc, mb, gb, am, mt);
        bam_de.barcode_demultiplex(wrong_fn, max_mismatch);
    }
    catch(const std::invalid_argument& e) 
    {
        tmp = true;
    }
    catch (...)
    {
        tmp = false;
    }
    return tmp;
}


bool test_barcode_demultiplex2()
{
    std::cout << "\ttest Bamdemultiplex::barcode_demultiplex() ...";
    //std::string tmp_out = "/home/users/allstaff/tian.l/public_data/GSM1544799/count";
    //std::string bamfn = "/home/users/allstaff/tian.l/public_data/GSM1544799/GSM1544799_SpeciesMix_HundredSTAMPs.bam";
    //std::string annofn = "/home/users/allstaff/tian.l/public_data/GSM1544799/GSM1544799_smp_list.csv";
    std::string tmp_out = "/Users/luyi/Downloads/dropseq";
    std::string bamfn = "/Users/luyi/Downloads/star_gene_exon_tagged_clean.bam";
    std::string annofn = "/Users/luyi/Downloads/anno_tmp.csv";
    std::string bc = "XC";
    std::string mb = "XM";
    std::string gb = "GE";
    std::string mt = "MT";
    std::string am;
    int max_mismatch = 1;
    Barcode bar;
    bar.read_anno(annofn);
    Bamdemultiplex bam_de = Bamdemultiplex(tmp_out, bar, bc, mb, gb, am, mt);
    bam_de.barcode_demultiplex(bamfn, max_mismatch);
    bam_de.write_statistics("overall_stat", "chr_stat", "cell_stat");
    // TODO: to be finished
    return true;
}


bool test_read_count()
{
    std::cout << "\ttest read_count() ...";
    char sep = ',';
    std::ifstream in_file("/Users/luyi/Downloads/dropseq/TCCGGGCTTAC.csv");
    std::unordered_map<std::string, std::vector<std::string>> tmp_res = read_count(in_file, sep);
    // TODO: to be finished
    return true;
}


bool test_UMI_correct1()
{
    std::cout << "\ttest UMI_correct1() ...";
    std::unordered_map<std::string, int> test_count;
    test_count["ATAATTA"] = 9;
    test_count["GTAGTAG"] = 6;
    test_count["ATAATTT"] = 1;
    int tmp_res = UMI_correct1(test_count);
    bool tmp = test_count["ATAATTA"] == 10 && \
        tmp_res == 1;
    return tmp;
}


bool test_UMI_dedup()
{
    std::cout << "\ttest UMI_dedup() ...";
    std::vector<std::string> v = {"ATGCTAAC", "ATCTGCCC", "ATGCTAAC", "GTAGTAG"};
    std::vector<std::string> v1 = {"ATGCTAAC", "ATGCTAAT", "ATGCTAAC", "GTAGTAG"};
    std::unordered_map<std::string, std::vector<std::string>> gene_read;
    gene_read["GENE01"] = v;
    gene_read["GENE02"] = v;
    gene_read["GENE03"] = v1;
    std::vector<int> UMI_dup_count(MAX_UMI_DUP+1);
    UMI_dedup_stat s = {};
    std::unordered_map<std::string, int> tmp_res;
    tmp_res = UMI_dedup(gene_read, UMI_dup_count, s, 1, true);
    bool tmp = tmp_res["GENE01"]== 3 && \
        tmp_res["GENE03"] == 2 && \
        UMI_dup_count[2] == 1 && \
        s.corrected_UMI == 1;
    return tmp;

}


bool test_get_counting_matrix()
{
    std::cout << "\ttest get_counting_matrix() ...";
    std::string in_dir = "/Users/luyi/Downloads/dropseq";
    std::string annofn = "/Users/luyi/Downloads/anno_tmp.csv";
    Barcode bar;
    bar.read_anno(annofn);
    get_counting_matrix(bar, in_dir, 1, true);
    // TODO: to be finished
    return true;

}


bool test_paired_fastq_to_bam()
{
    std::cout << "\ttest paired_fastq_to_bam() ...";
    read_s s = {};
    filter_s fl = {};

    s.id1_st = 0;
    s.id1_len = 8;
    s.id2_st = 6;
    s.id2_len = 8;
    s.umi_st = 0;
    s.umi_len = 6;

    fl.if_check_qual = true;
    fl.if_remove_N = true;
    fl.min_qual = 60;
    fl.num_below_min = 1;

    std::string bc = "XC";
    std::string mb = "XM";

    char *fq1_fn = (char *)"/Users/luyi/git/luyi_script/c_code/test_data/rd1.fq.gz";
    char *fq2_fn = (char *)"/Users/luyi/git/luyi_script/c_code/test_data/rd2.fq.gz";
    char *bam_out = (char *)"/Users/luyi/git/luyi_script/c_code/test_data/test_out.bam";


    //paired_fastq_to_bam(fq1_fn, fq2_fn, bam_out, s, fl, bc, mb);
    paired_fastq_to_bam(fq1_fn, fq2_fn, bam_out, s, fl);

    return true;

}


bool test_in_exon()
{
    std::cout << "\ttest Gene::in_exon() ...";
    Gene a = Gene("GENE001", 100,900, 1);
    a.add_exon(Interval(100,200,0));
    a.add_exon(Interval(400,500,0));
    a.add_exon(Interval(800,900,0));
    return a.in_exon(Interval(350,450,0));

}


bool test_parse_annotation()
{
    std::cout << "\ttest GeneAnnotation::parse_annotation() ...";
    std::string fn = "/Users/luyi/Downloads/Mus_musculus.GRCm38.83.gff3";
    GeneAnnotation anno = GeneAnnotation();
    anno.parse_gff3_annotation(fn, false);
    std::cout << anno << std::endl;
    //TODO: to be finished
    return true;

}


bool test_parse_align()
{
    std::cout << "\ttest Mapping::parse_align() ...";
    std::string gff3_fn = "/Users/luyi/Downloads/Mus_musculus.GRCm38.83.gff3";
    //std::string fn = "/Users/luyi/Downloads/10.bam";
    //std::string fn_out = "/Users/luyi/Downloads/10_mapped.bam";
    std::string fn = "/Users/luyi/Downloads/test_out.bam";
    std::string fn_out = "/Users/luyi/Downloads/test_mapped.bam";
    std::string bc = "YC";
    std::string mb = "YM";
    std::string ge = "GE";
    std::string am = "YE";
    Mapping a = Mapping();
    a.add_annotation(gff3_fn, false);
    a.parse_align(fn, fn_out, false, am, ge, bc, mb, 8, 6);
    //TODO: to be finished
    return true;

}


int main(int argc, char const *argv[])
{
    int passed_test = 0, failed_test = 0;
    bool test_result = false;

    // test join_path()
    test_result = test_join_path();
    print_result(test_result, passed_test, failed_test);
    
    // test hamming_distance()
    test_result = test_hamming_distance();
    print_result(test_result, passed_test, failed_test);

    // test vector_counter()
    test_result = test_vector_counter();
    print_result(test_result, passed_test, failed_test);

    // test Interval::overlap()
    test_result = test_overlap();
    print_result(test_result, passed_test, failed_test);

    // test Barcode::read_anno()
    test_result = test_read_anno();
    print_result(test_result, passed_test, failed_test);

    // test Barcode::get_closest_match()
    test_result = test_get_closest_match();
    print_result(test_result, passed_test, failed_test);

    // test Bamdemultiplex::Bamdemultiplex()
    test_result = test_barcode_demultiplex1();
    print_result(test_result, passed_test, failed_test);

    //test_result = test_barcode_demultiplex2();
    //print_result(test_result, passed_test, failed_test);

    // test read_count()
    test_result = test_read_count();
    print_result(test_result, passed_test, failed_test);

    // test test_UMI_correct1() method one for UMI correct
    test_result = test_UMI_correct1();
    print_result(test_result, passed_test, failed_test);

    // test UMI_dedup()
    test_result = test_UMI_dedup();
    print_result(test_result, passed_test, failed_test);

    // test get_counting_matrix()
    test_result = test_get_counting_matrix();
    print_result(test_result, passed_test, failed_test);

    test_result = test_paired_fastq_to_bam();
    print_result(test_result, passed_test, failed_test);

    // test Gene::in_exon()
    test_result = test_in_exon();
    print_result(test_result, passed_test, failed_test);

    //test_result = test_parse_annotation();
    //print_result(test_result, passed_test, failed_test);

    test_result = test_parse_align();
    print_result(test_result, passed_test, failed_test);

    // print all test result
    std::cout << "passed_test: " << passed_test << "; " 
    << "failed_test: " << failed_test << std::endl;
    return 0;
}