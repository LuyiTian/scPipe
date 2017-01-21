
# file can be *.gz or *.fastq
fq1=/home/users/allstaff/tian.l/daniela_nn29_040716/Undetermined_S0_R1.fastq
fq2=/home/users/allstaff/tian.l/daniela_nn29_040716/Undetermined_S0_R2.fastq

# no index1
index1_start=-1 
# trim 8bp of read1
index1_len=8

index2_start=6
index2_len=8
umi_start=0
umi_len=6

# experiment name
expr_name=NN29_16.07.2016

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
index_anno=/home/users/allstaff/tian.l/hisat_index/hisat_mm10_ERCC

samtools bam2fq $unaligned_bam > $out_dir/$expr_name.fastq
/home/users/allstaff/tian.l/hisat-master/hisat -x $index_anno -U $out_dir/$expr_name.fastq -S $out_dir/$expr_name.aligned.sam -5 5 -p 8
samtools view -bS $out_dir/$expr_name.aligned.sam -o $aligned_bam
#### map to transcriptome
exon_anno=/home/users/allstaff/tian.l/MM.GRCm38/Mus_musculus.GRCm38.84.gff3

# note: -BL should be inex1_len+index2_len. we dont have index1 at this experiment
# if the chr name in .bam file does not consistant with chr name in .gff3. add -C to fix the name issue
bin/sc_exon_mapping -O $mapped_bam -I $aligned_bam -A $exon_anno -BL $index2_len -UL $umi_len -S -C

#### demultiplexing

index_anno=/home/users/allstaff/tian.l/git/scQC/old_py_pipeline/NN29/NN29_sample_index_info.csv

bin/sc_demultiplex -I $mapped_bam -O $out_dir -A $index_anno

#### molecular counting

bin/sc_gene_counting -O $out_dir -A $index_anno


