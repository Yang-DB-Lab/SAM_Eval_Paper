# Trim reads with adapters with trim_galore
# Run the code in this way:
# cd /data/guang/SAM_evaluation_Paper/Trimming/Trim_Galore
# bash TrimGalore_50bpReadsTrimmingWithStringency2to4AndLengthStatistics.sh > TrimGalore_50bpReadsTrimmingWithStringency2to4AndLengthStatistics\_$(date +"%m-%d-%Y").log 2>&1

# Trim 50bp Reads

# list all fq files
for stringency_value in {2..4}

do
    cd /data/guang/SAM_evaluation_Paper/FakeSEReads50bpInclude0-25bpAdapters
    for fq_file in $(ls *.fq)
    do

    echo "Before trim_galore trimming of 50bp reads: $(date)"
    echo "Now processing: $fq_file with stringency $stringency_value"
    trim_galore --phred33 --illumina --stringency $stringency_value \
                --output_dir /data/guang/SAM_evaluation_Paper/Trimming/Trim_Galore/Fake50bpWithAdapters_Stringency$stringency_value $fq_file
    echo "After trim_galore trimming of 50bp reads: $(date)"
    done


    # Trimmed 50bp Reads
    cd /data/guang/SAM_evaluation_Paper/Trimming/Trim_Galore/Fake50bpWithAdapters_Stringency$stringency_value
    # list all fq files
    for fq_file in $(ls *.fq)
    do
    /data/guang/SAM_evaluation_Paper/Trimming/fastq_read_length_calculate $fq_file
    done
done
