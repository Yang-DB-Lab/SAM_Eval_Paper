
# Sampling chromosomes
cd /data/guang/Aligner_Benchmark_Paper/

## for i in {1..19}
## do
##     ./code/SamplingSequencesIn1KBpieces ./MouseChromosomes/GRCm38.100.chromosome.$i.fa ./MouseChromosome1kbPieces/GRCm38.100.chr$i.1k.pieces.tsv chr$i
## done

## for i in {X,Y,MT}
## do
##     ./code/SamplingSequencesIn1KBpieces ./MouseChromosomes/GRCm38.100.chromosome.$i.fa ./MouseChromosome1kbPieces/GRCm38.100.chr$i.1k.pieces.tsv chr$i
## done

for i in {1..19} X Y MT
do
    ./code/SamplingSequencesIn1KBpieces ./MouseChromosomes/GRCm38.100.chromosome.$i.fa ./MouseChromosome1kbPieces/GRCm38.100.chr$i.1k.pieces.tsv chr$i
done
# Merge all pieces of all chromosomes into one file
cat ./MouseChromosome1kbPieces/*.tsv > ./MouseChromosome1kbPieces/GRCm38.100.1kb.Pieces.tsv
