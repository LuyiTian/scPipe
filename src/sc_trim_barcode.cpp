//sc_trim_barcode.cpp

#include <iostream>
#include "trimbarcode.h"


int main(int argc, char* argv[]) {
    if (argc < 4) 
    {
        std::cout << "Usage: sc_trim_barcode \n" <<\
            "\t-O <outfile> the output bam file (required)\n"<<\
            "\t-R1 <read1_file> read 1(required)\n"<<\
            "\t-R2 <read2_file> read 2 for paired read (default: none)\n"<<\
            "\t-BS1 <id1_st> starting position of index in read one (default: 0)\n"<<\
            "\t-BL1 <id1_len> length of index in read one (default: 8)\n"<<\
            "\t-BS2 <id2_st> starting position of index in read two (default: 6)\n"<<\
            "\t-BL2 <id2_len> length of index in read two (default: 8)\n"<<\
            "\t-US <umi_st> starting position of UMI (default: 0)\n"<<\
            "\t-UL <umi_len> length of UMI (default: 6)\n"<<\
            "\t-QC <if_check_qual> whether removing low quality reads (default: 1)\n"<<\
            "\t-N <if_remove_N> remove reads that contains N in barcode and UMI (default: 1)\n"<<\
            "\t-MQ <min_qual> minimal quality of base (default: 30)\n"<<\
            "\t-MN <num_below_min> maximium number of base below quality (default: 1)\n"
            ;
        exit(0);
    } 
    else 
    { // if we got enough parameters...
    read_s s = {};
    filter_s fl = {};
    std::string fq1_fn, fq2_fn, bam_out;
    s.id1_st = 0;
    s.id1_len = 8;
    s.id2_st = 6;
    s.id2_len = 8;
    s.umi_st = 0;
    s.umi_len = 6;

    fl.if_check_qual = false;
    fl.if_remove_N = false;
    fl.min_qual = 30;
    fl.num_below_min = 1;
        for (int i = 1; i < argc; i++) 
        {
            std::string arg = argv[i];

            if (i + 1 != argc) 
            {
                if (arg == "-O") 
                {
                    bam_out = argv[i + 1];
                } 
                else if (arg == "-R1") 
                {
                    fq1_fn = argv[i + 1];
                }
                else if (arg == "-R2")
                {
                    fq2_fn = argv[i + 1];
                }
                else if (arg == "-BS1")
                {
                    s.id1_st = std::atoi(argv[i + 1]);
                }
                else if (arg == "-BL1")
                {
                    s.id1_len = std::atoi(argv[i + 1]);
                }
                else if (arg == "-BS2")
                {
                    s.id2_st = std::atoi(argv[i + 1]);
                }
                else if (arg == "-BL2")
                {
                    s.id2_len = std::atoi(argv[i + 1]);
                }
                else if (arg == "-US")
                {
                    s.umi_st = std::atoi(argv[i + 1]);
                }
                else if (arg == "-UL")
                {
                    s.umi_len = std::atoi(argv[i + 1]);
                }
                else if (arg == "-MQ")
                {
                    fl.min_qual = std::atoi(argv[i + 1]);
                }
                else if (arg == "-MN")
                {
                    fl.num_below_min = std::atoi(argv[i + 1]);
                }
            }
            if (arg == "-QC")
            {
                fl.if_check_qual = true;
            }
            else if (arg == "-N")
            {
                fl.if_remove_N = true;       
            }
        }
        std::cout << "######### trim barcode:" << std::endl;
        std::cout << "parameters:" << std::endl;
        std::cout << "\tout file: " << bam_out << std::endl;
        std::cout << "\tread1 file: " << fq1_fn << std::endl;
        std::cout << "\tread2 file: " << fq2_fn << std::endl;
        std::cout << "\tindex one start / length: " << s.id1_st <<" / " <<s.id1_len<<std::endl;
        std::cout << "\tindex two start / length: " << s.id2_st <<" / " <<s.id2_len<<std::endl;
        std::cout << "\tUMI start / length: " << s.umi_st <<" / " << s.umi_len<<std::endl;
        std::cout << "\tquality filter: " << fl.if_check_qual << std::endl;
        std::cout << "\tremove reads with N: " << fl.if_remove_N << std::endl;
        std::cout << "\tminimal quality: " << fl.min_qual << std::endl;
        std::cout << "\tmaximium number of base below quality: " << fl.num_below_min << std::endl;

        paired_fastq_to_bam((char *)fq1_fn.c_str(), (char *)fq2_fn.c_str(), (char *)bam_out.c_str(), s, fl);
        return 0;
    }
}