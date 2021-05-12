# Get into right folder
cd /data/guang/SAM_evaluation_Paper/

# fake fastq files with 30~50 bp reads without adapter
for length in {51..75}
do
# generate fastq files without adapter
./code/CreateFastqRead1File ./MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.fq $length AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTG 0


# bwa aln
# Align with "aln" algorithm
bwa aln -t 4 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.fq > /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.bwa

# Make single-end alignment to SAM file
bwa samse /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.bwa /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.fq > /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.Read1.bwa_aln.sam

/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_no_adapter_count /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.Read1.bwa_aln.sam

/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_no_adapter_sum /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.Read1.bwa_aln.sam

# bwa mem
bwa mem /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.fq > /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.Read1.bwa_mem.sam

/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_no_adapter_count /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.Read1.bwa_mem.sam

/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_no_adapter_sum /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.Read1.bwa_mem.sam

done


# Merge SumUp tsv files together
cd /data/guang/SAM_evaluation_Paper/FakeSEReads51-100bpNoAdapter
sum_tsv_files=$(ls *.different_mapping_count_SumUp.tsv)
cat $sum_tsv_files > merged_different_mapping_count_SumUp.tsv

