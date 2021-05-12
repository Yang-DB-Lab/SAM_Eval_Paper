# cd /data/guang/SAM_evaluation_Paper/FakeSEReads1kbto10kbNoAdapter
# bash Fake1kto10kbSEReads_NoAdapter_Align_and_Evaluate.sh > Fake1kto10kbSEReads_NoAdapter_Align_and_Evaluate\_$(date +"%m-%d-%Y").log 2>&1
# Get into right folder
cd /data/guang/SAM_evaluation_Paper/FakeSEReads1kbto10kbNoAdapter

# fake fastq files 1kb to 10kb reads without adapter
for length in $(seq 1000 1000 10000) 
do
# bwa mem
echo "Before bwa mem: $(date)"
echo "Now processing: GRCm38.100.genome.reads.$length\_bpwithoutAdapter.fq"
bwa mem -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakeSEReads1kbto10kbNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.fq > /data/guang/SAM_evaluation_Paper/FakeSEReads1kbto10kbNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.Read1.bwa_mem.sam
echo "After bwa mem: $(date)"

# count mapping of each read
/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_no_adapter_count /data/guang/SAM_evaluation_Paper/FakeSEReads1kbto10kbNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.Read1.bwa_mem.sam

/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_no_adapter_sum /data/guang/SAM_evaluation_Paper/FakeSEReads1kbto10kbNoAdapter/GRCm38.100.genome.reads.$length\_bpwithoutAdapter.Read1.bwa_mem.sam


done


# Merge SumUp tsv files together
cd /data/guang/SAM_evaluation_Paper/FakeSEReads1kbto10kbNoAdapter
sum_tsv_files=$(ls *.different_mapping_count_SumUp.tsv)
cat $sum_tsv_files > merged_different_mapping_count_SumUp_20200825.tsv

