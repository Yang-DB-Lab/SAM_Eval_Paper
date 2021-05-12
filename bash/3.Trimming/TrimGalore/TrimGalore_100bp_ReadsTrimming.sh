# Trim reads with adapters with trim_galore
# Run the code in this way:
# cd /data/guang/SAM_evaluation_Paper/Trimming/Trim_Galore
# bash TrimGalore_100bp_ReadsTrimming.sh > TrimGalore_100bp_ReadsTrimming\_$(date +"%m-%d-%Y").log 2>&1

# Trim 100bp Reads
cd /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-64bpAdapters
# list all fq files
for fq_file in $(ls *.fq)
do

echo "Before trim_galore trimming of 50bp reads: $(date)"
echo "Now processing: $fq_file"
trim_galore --phred33 --illumina --output_dir /data/guang/SAM_evaluation_Paper/Trimming/Trim_Galore/Fake100bpWithAdapters $fq_file
echo "After trim_galore trimming of 50bp reads: $(date)"
done



