cd /data/guang/Aligner_Benchmark_Paper/
mkdir MouseChromosomes
cd /data/guang/Aligner_Benchmark_Paper/MouseChromosomes/
for i in {1..19} X Y MT
do
    wget ftp://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.$i.fa.gz -O GRCm38.100.chromosome.$i.fa.gz
done

gunzip *.gz
