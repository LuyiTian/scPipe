#include <iostream>
#include "parsebam.h"
#include "cellbarcode.h"


int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cout << "Usage: sc_demultiplex \n" <<\
            "\t-I <infile> the input bam file (required)\n"<<\
            "\t-O <outdir> the output dir (required)\n"<<\
            "\t-A <barcode_annotation_file> annotate barcode and cell_id, \n for Drop-seq data you should run `./sc_detect_bc` to get the annotation file (required)\n"<<\
            "\t-M <max_mismatch> maximum mismatch allowed when matching barcode (default: 1)\n"<<\
            "\t-GE <gene_tag> two characters gene tag used in bam file (default: `GE`)\n"<<\
            "\t-MB <molecular_tag> two characters UMI tag used in bam file (default: `XM`)\n"<<\
            "\t-BC <cellular_tag> two characters cell barcode tag used in bam file (defalue: `XC`)\n"<<\
            "\t-MP <cellular_tag> two characters mapping status tag used in bam file (defalue: `YE`)\n"<<\
            "\t-MI <mitachondral_chromosome_name> should be consistant with the chromosome name in bam file.(defalue: `chrM`)\n"; 
        exit(0);
    } 
    else 
    { // if we got enough parameters...
        std::string bamfn, out_dir, annofn;
        std::string bc = "YC";
        std::string mb = "YM";
        std::string gb = "GE";
        std::string am = "YE";
        std::string mt = "chrM";
        int max_mismatch = 1;
        for (int i = 1; i < argc; i++) 
        {
            std::string arg = argv[i];

            if (i + 1 != argc) 
            {
                if (arg == "-I") 
                {
                    bamfn = argv[i + 1];
                } 
                else if (arg == "-O") 
                {
                    out_dir = argv[i + 1];
                } 
                else if (arg == "-A") 
                {
                    annofn = argv[i + 1];
                }
                else if (arg == "-M")
                {
                    max_mismatch = argv[i + 1][0]-'0';
                }
                else if (arg == "-GE")
                {
                    gb = argv[i + 1];
                }
                else if (arg == "-MB")
                {
                    mb = argv[i + 1];
                }
                else if (arg == "-BC")
                {
                    bc = argv[i + 1];
                }
                else if (arg == "-MI")
                {
                    mt = argv[i + 1];
                }
            }
        }
        std::cout << "######### demultiplexing:" << std::endl;
        std::cout << "parameters:" << std::endl;
        std::cout << "\tbam file: " << bamfn << std::endl;
        std::cout << "\tout dir: " << out_dir << std::endl;
        std::cout << "\tannotation file: " << annofn << std::endl;
        std::cout << "\tmax mismatch: " << max_mismatch << std::endl;
        std::cout << "\tcell/molecular/gene tag: " << bc << "/"<< mb <<"/" << gb << std::endl;

        Barcode bar;
        bar.read_anno(annofn);
        Bamdemultiplex bam_de = Bamdemultiplex(out_dir, bar, bc, mb, gb, am, mt);
        bam_de.barcode_demultiplex(bamfn, max_mismatch);
        bam_de.write_statistics("overall_stat", "chr_stat", "cell_stat");
        return 0;
    }
}