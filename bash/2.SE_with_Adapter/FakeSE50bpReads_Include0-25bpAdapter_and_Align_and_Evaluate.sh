
# Run the code in this way:
#
# cd /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters
# bash FakeSE50bpReads_Include0-25bpAdapter_and_Align_and_Evaluate.sh > FakeSE50bpReads_Include0-25bpAdapter_and_Align_and_Evaluate\_$(date +"%m-%d-%Y").log 2>&1
#
# Get into the right folder
cd /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters

# fake fastq file of 100bp reads Include 0~40bp adapter
for adapterLength in {0..25}
do
# Generate fastq file Include adapters of different length
# generate fastq files of 100bp reads Include 1-40bp adapter
# Use "Index Adapter 5" from "TruSeq DNA and RNA Index Adapters"(TruSeq Single Indexes)
# NOTE: The first character of the adapter to be added to reads is 'A' !!!
/data/guang/SAM_evaluation_Paper/code/CreateFastqRead1File /data/guang/SAM_evaluation_Paper/MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters/GRCm38.100.genome.reads.50bpInclude_$adapterLength\_bp_Adapter.fq 50 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG $adapterLength
# do bwa mem alignment
echo "Before bwa mem: $(date)"
echo "Now processing: GRCm38.100.genome.reads.50bpInclude_$adapterLength\_bp_Adapter.fq"
bwa mem -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters/GRCm38.100.genome.reads.50bpInclude_$adapterLength\_bp_Adapter.fq > /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters/GRCm38.100.genome.reads.50bpInclude_$adapterLength\_bp_Adapter.read1.bwamem.se.sam
echo "After bwa mem: $(date)"

done



# bwa aln
# Align with "aln" algorithm
for adapterLength in {0..25}
do
echo "Before bwa aln: $(date)"
echo "Now processing: GRCm38.100.genome.reads.50bpInclude_$adapterLength\_bp_Adapter.fq"
bwa aln -t 4 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters/GRCm38.100.genome.reads.50bpInclude_$adapterLength\_bp_Adapter.fq > /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters/GRCm38.100.genome.reads.50bpInclude_$adapterLength\_bp_Adapter.bwa

# Make single-end alignment to SAM file
bwa samse /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters/GRCm38.100.genome.reads.50bpInclude_$adapterLength\_bp_Adapter.bwa /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters/GRCm38.100.genome.reads.50bpInclude_$adapterLength\_bp_Adapter.fq > /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters/GRCm38.100.genome.reads.50bpInclude_$adapterLength\_bp_Adapter.Read1.bwa_aln.se.sam
echo "After bwa aln: $(date)"
done


# Use full path!!! Otherwise the $sam file may not be found and through error(hopefully this is the reason).
# Evaluate mapping quality of bwa mem generated SAM files
cd /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters
bwamem_sams=$(ls *.sam | grep bwamem)
for sam in $bwamem_sams
do
echo "Evaluate $sam"
/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_with_adapter_count /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters/$sam AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG mem
done

# Evaluate SAM files from bwa aln
cd /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters
bwa_aln_sams=$(ls *.sam | grep bwa_aln)
for sam in $bwa_aln_sams
do
echo "Evaluate $sam"
/data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_with_adapter_count /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters/$sam AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG aln
done


