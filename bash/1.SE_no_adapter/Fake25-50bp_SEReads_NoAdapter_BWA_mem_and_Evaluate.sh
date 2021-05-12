# Run the script in this way
# cd /data/guang/SAM_evaluation_Paper/FakeSEReads25-50bpNoAdapter/
# bash Fake25-50bp_SEReads_NoAdapter_BWA_mem_and_Evaluate.sh > Fake25-50bp_SEReads_NoAdapter_BWA_mem_and_Evaluate\_$(date +"%m-%d-%Y").log 2>&1

cd /data/guang/SAM_evaluation_Paper/FakeSEReads25-50bpNoAdapter

for length in {25..50}
do

# Align with "aln" algorithm
echo "Before bwa mem: $(date)"
echo "Now processing: GRCm38.100.genome.reads.$length\_bpwithoutAdapter.fq"
bwa mem -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakeSEReads25-50bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.fq > /data/guang/SAM_evaluation_Paper/FakeSEReads25-50bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.Read1.bwa_mem.sam
echo "After bwa mem: $(date)"

# Evaluate SAM file
/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_no_adapter_count /data/guang/SAM_evaluation_Paper/FakeSEReads25-50bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.Read1.bwa_mem.sam

/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_no_adapter_sum /data/guang/SAM_evaluation_Paper/FakeSEReads25-50bpNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.Read1.bwa_mem.sam

done

sum_tsv_files=$(ls *.bwa_mem.sam.different_mapping_count_SumUp.tsv)
cat $sum_tsv_files > merged_bwa_mem_different_mapping_count_SumUp_20200918.tsv
