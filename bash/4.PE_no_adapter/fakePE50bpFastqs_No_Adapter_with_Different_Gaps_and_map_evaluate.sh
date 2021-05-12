
# Run this script in this way to redirect output to a log file with the day it is operated.
# bash /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/fakePE50bpFastqs_No_Adapter_with_Different_Gaps_and_map_evaluate.sh > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/fakePE50bpFastqs_No_Adapter_with_Different_Gaps_and_map_evaluate\_$(date +"%m-%d-%Y").log 2>&1


# Whole genome sampled DNA pieces with 2bp i7 adpter
cd /data/guang/SAM_evaluation_Paper/
# 50bp reads without Adapter but different gaps
gaps1=$(seq -50 20 50)
gap0=0
gaps2=$(seq 100 100 900)

for gap in $gaps1
do
./code/CreatePairedEndFastqs ./MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap 50 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTG 0 $gap
done

for gap in $gaps2
do
./code/CreatePairedEndFastqs ./MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap 50 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTG 0 $gap
done

./code/CreatePairedEndFastqs ./MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap 50 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTG 0 $gap0

# Alignment with 'bwa aln'
for gap in $gaps1
do
echo "before bwa aln mapping GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.sai

bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R2.sai

bwa sampe /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R2.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.aln.pe.sam
echo "After bwa aln mapping GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
done

for gap in $gaps2
do
echo "before bwa aln mapping GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.sai

bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R2.sai

bwa sampe /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R2.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.aln.pe.sam
echo "After bwa aln mapping GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
done

echo "before bwa aln mapping GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R1.fq and R2.fq $(date)"
bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R1.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R1.sai

bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R2.sai

bwa sampe /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R1.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R2.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.aln.pe.sam
echo "After bwa aln mapping GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R1.fq and R2.fq $(date)"



# Try to restrict insert length
gaps1=$(seq -50 20 50)
gap0=0
gaps2=$(seq 100 100 900)

for gap in $gaps1
do
echo "before bwa sampe mapping with -a 10000 for GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
bwa sampe -a 10000 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R2.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.aln.a10k.pe.sam
echo "before bwa sampe mapping with -a 10000 for GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
done


for gap in $gaps2
do
echo "before bwa sampe mapping with -a 10000 for GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
bwa sampe -a 10000 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R2.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.aln.a10k.pe.sam
echo "After bwa sampe mapping with -a 10000 for GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
done

echo "before bwa sampe mapping with -a 10000 for GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R1.fq and R2.fq $(date)"
bwa sampe -a 10000 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R1.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R2.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.aln.a10k.pe.sam
echo "After bwa sampe mapping with -a 10000 for GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap0\_bpGap.R1.fq and R2.fq $(date)"


# Check -25bp gap
gap_minus25=-25
./code/CreatePairedEndFastqs ./MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap_minus25\_bpGap 50 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTG 0 $gap_minus25
bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap_minus25\_bpGap.R1.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap_minus25\_bpGap.R1.sai

echo "Before bwa aln mapping GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap_minus25\_bpGap.R1.fq and R2.fq $(date)"
bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap_minus25\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap_minus25\_bpGap.R2.sai

bwa sampe /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap_minus25\_bpGap.R1.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap_minus25\_bpGap.R2.sai /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap_minus25\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap_minus25\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap/GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap_minus25\_bpGap.aln.pe.sam
echo "After bwa aln mapping GRCm38.100.genome.reads.50bpwithoutAdapterAnd$gap_minus25\_bpGap.R1.fq and R2.fq $(date)"

# Count mapping tags
cd /data/guang/SAM_evaluation_Paper/FakePE50noAdapter_differentGap
sam_files=$(ls *.sam)
for sam in $sam_files
do
# first get gap length
gap=$(sed 's/.*And\(.*\)\_bpGap.*/\1/' <<< $sam)
/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/PE_no_adapter_count $sam 50 $gap
done

