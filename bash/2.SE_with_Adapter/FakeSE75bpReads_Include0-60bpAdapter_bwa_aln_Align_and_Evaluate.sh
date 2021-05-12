

# Run the code in this way:
#
# cd /data/guang/SAM_evaluation_Paper/FakeSEReads75bpInclude0-60bpAdapters
# bash FakeSE75bpReads_Include0-60bpAdapter_bwa_aln_Align_and_Evaluate.sh > FakeSE75bpReads_Include0-60bpAdapter_bwa_aln_Align_and_Evaluate\_$(date +"%m-%d-%Y").log 2>&1
#
# Get into the right folder
cd /data/guang/SAM_evaluation_Paper/FakeSEReads75bpInclude0-60bpAdapters


for adapterLength in {0..60}
do

# do bwa mem alignment
echo "Before bwa aln: $(date)"
echo "Now processing: GRCm38.100.genome.reads.75bpInclude_$adapterLength\_bp_Adapter.fq using bwa aln"
bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakeSEReads75bpInclude0-60bpAdapters/GRCm38.100.genome.reads.75bpInclude_$adapterLength\_bp_Adapter.fq > /data/guang/SAM_evaluation_Paper/FakeSEReads75bpInclude0-60bpAdapters/GRCm38.100.genome.reads.75bpInclude_$adapterLength\_bp_Adapter.bwa
# Make single-end alignment to SAM file
bwa samse /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex \
          /data/guang/SAM_evaluation_Paper/FakeSEReads75bpInclude0-60bpAdapters/GRCm38.100.genome.reads.75bpInclude_$adapterLength\_bp_Adapter.bwa \
          /data/guang/SAM_evaluation_Paper/FakeSEReads75bpInclude0-60bpAdapters/GRCm38.100.genome.reads.75bpInclude_$adapterLength\_bp_Adapter.fq \
          > /data/guang/SAM_evaluation_Paper/FakeSEReads75bpInclude0-60bpAdapters/GRCm38.100.genome.reads.75bpInclude_$adapterLength\_bp_Adapter.read1.bwa_aln.se.sam
echo "After bwa aln: $(date)"

done

# Evaluate mapping quality of each SAM
cd /data/guang/SAM_evaluation_Paper/FakeSEReads75bpInclude0-60bpAdapters
for sam in $(ls *.sam)
do
/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_with_adapter_count $sam AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG mem
done


