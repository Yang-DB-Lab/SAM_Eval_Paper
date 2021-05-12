# Trim reads with adapters with trim_galore
# Run the code in this way:
# cd /data/guang/SAM_evaluation_Paper/Trimming/Trimmomatic
# bash trimmomatic_100bpReadsTrimmedLengthStatistics.sh > trimmomatic_100bpReadsTrimmedLengthStatistics\_$(date +"%m-%d-%Y").log 2>&1

# Trimmed 100bp Reads
cd /data/guang/SAM_evaluation_Paper/Trimming/Trimmomatic/Fake100bpWithAdapters
# list all fq files
for fq_file in $(ls *.fq)
do
/data/guang/SAM_evaluation_Paper/Trimming/fastq_read_length_calculate $fq_file
done

