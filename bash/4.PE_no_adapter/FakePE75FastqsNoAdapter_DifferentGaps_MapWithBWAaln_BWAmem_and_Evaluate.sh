# Run this script in this way to redirect output to a log file with the day it is operated.
# bash /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/FakePE75FastqsNoAdapter_DifferentGaps_MapWithBWAaln_BWAmem_and_Evaluate.sh > /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/FakePE75FastqsNoAdapter_DifferentGaps_MapWithBWAaln_BWAmem_and_Evaluate\_$(date +"%m-%d-%Y").log 2>&1


# Whole genome sampled DNA pieces with 2bp i7 adpter
cd /data/guang/SAM_evaluation_Paper/
# 50bp reads without Adapter but different gaps
gaps1=$(seq -75 25 75)
gaps2=$(seq 100 100 800)

for gap in $gaps1
do
./code/CreatePairedEndFastqs ./MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap 75 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTG 0 $gap
done

for gap in $gaps2
do
./code/CreatePairedEndFastqs ./MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap 75 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTG 0 $gap
done

# Alignment with 'bwa aln'
for gap in $gaps1
do
echo "before bwa aln mapping GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq > /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.sai

bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R2.sai

bwa sampe /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.sai /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R2.sai /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.aln.pe.sam
echo "After bwa aln mapping GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
done

for gap in $gaps2
do
echo "before bwa aln mapping GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq > /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.sai

bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R2.sai

bwa sampe /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.sai /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R2.sai /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.aln.pe.sam
echo "After bwa aln mapping GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
done




# BWA mem

for gap in $gaps1
do
echo "before bwa mem mapping of GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
bwa mem -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.bwa.mem.pe.sam
echo "before bwa mem mapping of GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
done


for gap in $gaps2
do
echo "before bwa mem mapping of GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
bwa mem -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps/GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.bwa.mem.pe.sam
echo "before bwa mem mapping of GRCm38.100.genome.reads.75bpwithoutAdapterAnd$gap\_bpGap.R1.fq and R2.fq $(date)"
done


# Count mapping tags
cd /data/guang/SAM_evaluation_Paper/FakePE75NoAdapter_DifferentGaps
sam_files=$(ls *.sam)
for sam in $sam_files
do
# first get gap length
gap=$(sed 's/.*And\(.*\)\_bpGap.*/\1/' <<< $sam)
/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/PE_no_adapter_count $sam 75 $gap
done

