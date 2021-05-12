# Trim reads with adapters with trim_galore
# Run the code in this way:
# cd /data/guang/SAM_evaluation_Paper/Trimming/Trimmomatic
# bash trimmomatic_100bpReadsTrimming.sh > trimmomatic_100bpReadsTrimming\_$(date +"%m-%d-%Y").log 2>&1

# Trim 100bp Reads
cd /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-64bpAdapters
# list all fq files
for fq_file in $(ls *.fq)
do

echo "Before Trimmomatic trimming of 100bp reads: $(date)"
echo "Now processing: $fq_file"
java -jar /home/guang/bio_softwares/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 $fq_file \
/data/guang/SAM_evaluation_Paper/Trimming/Trimmomatic/Fake100bpWithAdapters/$fq_file.trimmomatic_trimmed.fq \
ILLUMINACLIP:/home/guang/bio_softwares/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
echo "After Trimmomatic trimming of 100bp reads: $(date)"
done

