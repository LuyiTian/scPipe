
# file can be *.gz or *.fastq
fq1="/wehisan/general/user_managed/grpu_naik.s_2/2016 Sequencing Runs/NN34/Daniela_nn34_NEXS1286_120816/Undetermined_S0_R2.fastq"
fq2="/wehisan/general/user_managed/grpu_naik.s_2/2016 Sequencing Runs/NN34/Daniela_nn34_NEXS1286_120816/Undetermined_S0_R1.fastq"

# no index1
index1_start=-1 
# trim 2bp of read1
index1_len=2

index2_start=6
index2_len=8
umi_start=0
umi_len=6

# experiment name
expr_name=NN34_15.08.2016

out_dir=/wehisan/general/user_managed/grpu_naik.s_1/$expr_name

mkdir -p $out_dir
mkdir -p $out_dir/count
mkdir -p $out_dir/stat

unaligned_bam=$out_dir/$expr_name.bam
aligned_bam=$out_dir/$expr_name.aligned.bam
mapped_bam=$out_dir/$expr_name.aligned.mapped.bam



########### pipeline:

#### trim read. read qc
bin/sc_trim_barcode -O $unaligned_bam -R1 $fq1 -R2 $fq2 -BS1 $index1_start -BL1 $index1_len \
    -BS2 $index2_start -BL2 $index2_len -US $umi_start -UL $umi_len -QC -N

#### alignment using HISAT
index_anno=/home/users/allstaff/tian.l/MM.GRCm38/subRead_MM_GRCm38_ERCC_phiX_WPRE
subread-align -T 12 -i $index_anno -r $unaligned_bam -o $aligned_bam

#### map to transcriptome
exon_anno=/home/users/allstaff/tian.l/MM.GRCm38/Mus_musculus.GRCm38.84.gff3
ERCC_anno=/wehisan/general/user_managed/grpu_naik.s_1/Annotation/ercc_anno.txt

# note: -BL should be inex1_len+index2_len. we dont have index1 at this experiment
# if the chr name in .bam file does not consistant with chr name in .gff3. add -C to fix the name issue
bin/sc_exon_mapping -O $mapped_bam -I $aligned_bam -A $exon_anno,$ERCC_anno -BL $index2_len -UL $umi_len -S

#### demultiplexing

index_anno="/home/users/allstaff/tian.l/git/scQC/old_py_pipeline/NN34 sample info.csv"

bin/sc_demultiplex -I $mapped_bam -O $out_dir -A $index_anno

#### molecular counting

bin/sc_gene_counting -O $out_dir -A $index_anno


