# Trim reads with adapters with trim_galore
# Run the code in this way:
# cd /data/guang/SAM_evaluation_Paper/Trimming/Trimmomatic
# bash trimmomatic_50bpReadsTrimmingAndTrimmedLengthStatistics.sh > trimmomatic_50bpReadsTrimmingAndTrimmedLengthStatistics\_$(date +"%m-%d-%Y").log 2>&1

# Trim 50bp Reads
cd /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters
# list all fq files
for fq_file in $(ls *.fq)
do

echo "Before Trimmomatic trimming of 50bp reads: $(date)"
echo "Now processing: $fq_file"
java -jar /home/guang/bio_softwares/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 $fq_file \
/data/guang/SAM_evaluation_Paper/Trimming/Trimmomatic/Fake50bpWithAdapters/$fq_file.trimmomatic_trimmed.fq \
ILLUMINACLIP:/home/guang/bio_softwares/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
echo "After Trimmomatic trimming of 50bp reads: $(date)"
done

# Trimmed 50bp Reads
cd /data/guang/SAM_evaluation_Paper/Trimming/Trimmomatic/Fake50bpWithAdapters
# list all fq files
for fq_file in $(ls *.fq)
do
/data/guang/SAM_evaluation_Paper/Trimming/fastq_read_length_calculate $fq_file
done
