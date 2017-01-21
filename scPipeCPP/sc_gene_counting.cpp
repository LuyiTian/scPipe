// gene_counting.cpp

#include <iostream>
#include "parsecount.h"
#include "cellbarcode.h"


int main(int argc, char* argv[]) {
    if (argc < 4) 
    {
        std::cout << "Usage: sc_demultiplex \n" <<\
            "\t-O <outdir> the output dir (required)\n"<<\
            "\t-A <barcode_annotation_file> annotate barcode and cell_id, \n for Drop-seq data you should run `./sc_detect_bc` to get the annotation file (required)\n"<<\
            "\t-U <UMI_correction> UMI correction methods, 0 mean no correction. \n for the detail of each methods see the documentation (default: 1)\n"<<\
            "\t-F <gene_filter> whether to filter low coverage genes (default: false)\n";
        exit(0);
    } 
    else 
    { // if we got enough parameters...
        std::string out_dir, annofn;
        int UMI_correction = 1;
        bool gene_filter = false;
        for (int i = 1; i < argc; i++) 
        {
            std::string arg = argv[i];

            if (i + 1 != argc) 
            {
                if (arg == "-O") 
                {
                    out_dir = argv[i + 1];
                } 
                else if (arg == "-A") 
                {
                    annofn = argv[i + 1];
                }
                else if (arg == "-U")
                {
                    UMI_correction = argv[i + 1][0]-'0';
                }
            }
            if (arg == "-F")
            {
                gene_filter = true;                    
            }
        }
        std::cout << "######### molecular counting:" << std::endl;
        std::cout << "parameters:" << std::endl;
        std::cout << "\tout dir: " << out_dir << std::endl;
        std::cout << "\tannotation file: " << annofn << std::endl;
        std::cout << "\tUMI_correction: " << UMI_correction << std::endl;
        std::cout << "\tgene_filter " << gene_filter << std::endl;

        Barcode bar;
        bar.read_anno(annofn);
        get_counting_matrix(bar, out_dir, UMI_correction, gene_filter);
        return 0;
    }
}