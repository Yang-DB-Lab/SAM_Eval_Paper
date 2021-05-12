
# Run the code in this way:
#
# cd /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters
# bash FakeSE100bpReads_Include0-40bpAdapter_and_Align_and_Evaluate.sh > FakeSE100bpReads_Include0-40bpAdapter_and_Align_and_Evaluate\_$(date +"%m-%d-%Y").log 2>&1
#
# Get into the right folder
cd /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters

# fake fastq file of 100bp reads Include 0~40bp adapter
for adapterLength in {0..40}
do
# Generate fastq file Include adapters of different length
# generate fastq files of 100bp reads Include 1-40bp adapter
# Use "Index Adapter 5" from "TruSeq DNA and RNA Index Adapters"(TruSeq Single Indexes)
/data/guang/SAM_evaluation_Paper/code/CreateFastqRead1File /data/guang/SAM_evaluation_Paper/MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters/GRCm38.100.genome.reads.100bpInclude_$adapterLength\_bp_Adapter.fq 100 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG $adapterLength
# do bwa mem alignment
echo "Before bwa mem: $(date)"
echo "Now processing: GRCm38.100.genome.reads.100bpInclude_$adapterLength\_bp_Adapter.fq"
bwa mem -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters/GRCm38.100.genome.reads.100bpInclude_$adapterLength\_bp_Adapter.fq > /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters/GRCm38.100.genome.reads.100bpInclude_$adapterLength\_bp_Adapter.read1.bwamem.se.sam
echo "After bwa mem: $(date)"

done


# Evaluate mapping quality of each SAM
cd /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters
for sam in $(ls *.sam)
do
/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_with_adapter_count $sam AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG mem
done



