# Whole genome sampled DNA pieces with 2bp i7 adpter
cd /data/guang/SAM_evaluation_Paper/
# 50bp reads without Adapter but different gaps
gaps1=$(seq -100 25 100)
gap0=0
gaps2=$(seq 200 100 800)

for gap in $gaps1
do
./code/CreatePairedEndFastqs ./MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps/GRCm38.100.genome.reads.100bpwithoutAdapterAnd$gap\_bpGap 100 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTG 0 $gap
done

for gap in $gaps2
do
./code/CreatePairedEndFastqs ./MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps/GRCm38.100.genome.reads.100bpwithoutAdapterAnd$gap\_bpGap 100 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTG 0 $gap
done

./code/CreatePairedEndFastqs ./MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps/GRCm38.100.genome.reads.100bpwithoutAdapterAnd$gap0\_bpGap 100 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTG 0 $gap0

# Alignment with 'bwa mem'
for gap in $gaps1
do

bwa mem /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps/GRCm38.100.genome.reads.100bpwithoutAdapterAnd$gap\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps/GRCm38.100.genome.reads.100bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps/GRCm38.100.genome.reads.100bpwithoutAdapterAnd$gap\_bpGap.bwamem.pe.sam
done

# gaps2 fqs
for gap in $gaps2
do
bwa mem /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps/GRCm38.100.genome.reads.100bpwithoutAdapterAnd$gap\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps/GRCm38.100.genome.reads.100bpwithoutAdapterAnd$gap\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps/GRCm38.100.genome.reads.100bpwithoutAdapterAnd$gap\_bpGap.bwamem.pe.sam
done

# Align gap0 fq
bwa mem /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps/GRCm38.100.genome.reads.100bpwithoutAdapterAnd$gap0\_bpGap.R1.fq /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps/GRCm38.100.genome.reads.100bpwithoutAdapterAnd$gap0\_bpGap.R2.fq > /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps/GRCm38.100.genome.reads.100bpwithoutAdapterAnd$gap0\_bpGap.bwamem.pe.sam

# Count mapping tags
cd /data/guang/SAM_evaluation_Paper/FakePE100_NoAdapter_DifferentGaps
sam_files=$(ls *.sam)
for sam in $sam_files
do
# first get gap length
gap=$(sed 's/.*And\(.*\)\_bpGap.*/\1/' <<< $sam)
/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/PE_no_adapter_count $sam 100 $gap
done

