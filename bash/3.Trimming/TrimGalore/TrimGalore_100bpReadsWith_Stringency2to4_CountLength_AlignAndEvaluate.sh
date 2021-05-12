# Trim reads with adapters with trim_galore
# Run the code in this way:
# cd /data/guang/SAM_evaluation_Paper/Trimming/Trim_Galore/
# bash TrimGalore_100bpReadsWith_Stringency2to4_CountLength_AlignAndEvaluate.sh > TrimGalore_100bpReadsWith_Stringency2to4_CountLength_AlignAndEvaluate\_$(date +"%m-%d-%Y").log 2>&1

# Trim 100bp Reads

# list all fq files
for stringency_value in {2..5}
do
    # Get into folder containing 100bp read fastq files with adapter sequences
    cd /data/guang/SAM_evaluation_Paper/FakeSEReads100bpInclude0-64bpAdapters
    for fq_file in $(ls *.fq)
    do

        echo "Before trim_galore trimming of 100bp reads: $(date)"
        echo "Now processing: $fq_file with stringency $stringency_value"
        trim_galore --phred33 --illumina --stringency $stringency_value \
                    --output_dir /data/guang/SAM_evaluation_Paper/Trimming/Trim_Galore/Fake100bpWithAdapters_Stringency$stringency_value $fq_file
        echo "After trim_galore trimming of 100bp reads: $(date)"
    done

    # Count, align trimmed fastq files and evaluate aligned SAM files
    # Get into folder containing trimmed 100bp Reads
    cd /data/guang/SAM_evaluation_Paper/Trimming/Trim_Galore/Fake100bpWithAdapters_Stringency$stringency_value
    # list all fq files
    for fq_file in $(ls *.fq)
    do
        # (1) count lengths of reads after trimming
        /data/guang/SAM_evaluation_Paper/Trimming/fastq_read_length_calculate $fq_file
        
        # (2) do bwa aln alignment
        echo "Before bwa aln: $(date)"
        echo "Now processing: $fq_file"
        
        bwa aln -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex $fq_file > $fq_file.bwa_aln.se.bwa
        bwa samse /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex $fq_file.bwa_aln.se.bwa $fq_file > $fq_file.bwa_aln.se.sam

        echo "After bwa aln: $(date)"
        
        # (3) do bwa mem alignment
        echo "Before bwa mem: $(date)"
        echo "Now processing: $fq_file"
        bwa mem -t 6 /data/guang/mouse_genome_index/bwa_index_chr1-19XYMT/GRCm38.100.chromosomes.bwaIndex $fq_file > $fq_file.bwamem.se.sam
        echo "After bwa mem: $(date)"
    done

    # After alignment evaluate bwa_aln sam files in the foloder
    for bwa_aln_sam in $(ls *.bwa_aln.se.sam)
    do
        /data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_with_adapter_count $bwa_aln_sam AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG aln
    done
    
    # After alignment evaluate bwamem sam files in the foloder    
    for bwa_mem_sam in $(ls *.bwamem.se.sam)
    do
        /data/guang/SAM_evaluation_Paper/code/SAM_Count_Statistics/SE_with_adapter_count $bwa_mem_sam AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG mem
    done
done
