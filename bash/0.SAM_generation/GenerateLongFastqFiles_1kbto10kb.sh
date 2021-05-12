
# Sampling chromosomes
cd /data/guang/SAM_evaluation_Paper/FakeSEReads1kbto10kbNoAdapter

## for i in {1..19}
## do
##     ./code/SamplingSequencesIn1KBpieces ./MouseChromosomes/GRCm38.100.chromosome.$i.fa ./MouseChromosome1kbPieces/GRCm38.100.chr$i.1k.pieces.tsv chr$i
## done

## for i in {X,Y,MT}
## do
##     ./code/SamplingSequencesIn1KBpieces ./MouseChromosomes/GRCm38.100.chromosome.$i.fa ./MouseChromosome1kbPieces/GRCm38.100.chr$i.1k.pieces.tsv chr$i
## done
for length in $(seq 1000 1000 10000) 
#1k to 10k sequential read lengths

do 
# Work on each chromosome to generate fastq files with $length long read
    for i in {1..19} X Y MT
    do
        /data/guang/SAM_evaluation_Paper/code/GenerateLongReadFastqFileFromOneChromosomeWithGivenStartCoordinateAndLength \
        /data/guang/SAM_evaluation_Paper/MouseChromosomes/GRCm38.100.chromosome.$i.fa \
        /data/guang/SAM_evaluation_Paper/MouseChromosome1kbPieces/GRCm38.100.chr$i.1k.pieces.tsv $length chr$i.$length\_bp.fastq

    done
    # Merge all pieces of from each chromosome into one final fastq file
    cat ./*.fastq > ./GRCm38.100.genome.reads.$length\_bpwithoutAdapter.fq
    rm ./*.fastq
done
