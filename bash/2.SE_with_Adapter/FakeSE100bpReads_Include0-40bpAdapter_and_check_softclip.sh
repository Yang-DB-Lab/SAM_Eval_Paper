
# Run the code in this way:
#
# cd /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters
# bash FakeSE100bpReads_Include0-40bpAdapter_and_check_softclip.sh > FakeSE100bpReads_Include0-40bpAdapter_and_check_softclip\_$(date +"%m-%d-%Y").log 2>&1 

#
# Get into the right folder
cd /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters

# fake fastq file of 100bp reads Include 0~40bp adapter
for adapterLength in {0..40}
do


# do bwa mem alignment
echo "Before bwa mem -L 1: $(date)"
echo "Now processing: GRCm38.100.genome.reads.100bpInclude_$adapterLength\_bp_Adapter.fq with -L 1"
bwa mem -t 6 -L 1 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters/GRCm38.100.genome.reads.100bpInclude_$adapterLength\_bp_Adapter.fq > /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters/GRCm38.100.genome.reads.100bpInclude_$adapterLength\_bp_Adapter.read1.bwamem_L1.se.sam
echo "After bwa mem -L 1: $(date)"

done

for adapterLength in {0..40}
do


# do bwa mem alignment
echo "Before bwa mem -L 100: $(date)"
echo "Now processing: GRCm38.100.genome.reads.100bpInclude_$adapterLength\_bp_Adapter.fq with -L 1"
bwa mem -t 6 -L 100 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters/GRCm38.100.genome.reads.100bpInclude_$adapterLength\_bp_Adapter.fq > /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters/GRCm38.100.genome.reads.100bpInclude_$adapterLength\_bp_Adapter.read1.bwamem_L100.se.sam
echo "After bwa mem -L 100: $(date)"

done


# Evaluate mapping quality of each SAM
cd /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-40bpAdapters
for clip_sam in $(ls *.sam | grep "_L")
do
/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_with_adapter_count $clip_sam AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG mem
done



