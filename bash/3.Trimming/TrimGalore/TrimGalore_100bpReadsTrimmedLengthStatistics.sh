# Trim reads with adapters with trim_galore
# Run the code in this way:
# cd /data/guang/SAM_evaluation_Paper/Trimming/Trim_Galore
# bash TrimGalore_100bpReadsTrimmedLengthStatistics.sh > TrimGalore_100bpReadsTrimmedLengthStatistics\_$(date +"%m-%d-%Y").log 2>&1

# Trimmed 100bp Reads
cd /data/guang/SAM_evaluation_Paper/Trimming/Trim_Galore/Fake100bpWithAdapters
# list all fq files
for fq_file in $(ls *.fq)
do
/data/guang/SAM_evaluation_Paper/Trimming/fastq_read_length_calculate $fq_file
done

