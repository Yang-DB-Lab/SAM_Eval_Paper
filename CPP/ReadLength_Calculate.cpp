#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// example of fastq filename:
// GRCm38.100.genome.reads.100bpInclude_0_bp_Adapter_trimmed.fq
int main(int argc, char* argv[])
{
    std::string fastq_filename = argv[1];

    int read_length=0;
    std::size_t reads_loci_found = fastq_filename.find("reads.");
    std::size_t bpInclude_loci_found = fastq_filename.find("bpInclude");
    int read_length_stringLen = bpInclude_loci_found - reads_loci_found - 6; // "reads." length = 6
    std::string read_length_str = fastq_filename.substr(reads_loci_found + 6, read_length_stringLen);
    read_length = std::stoi(read_length_str);
    std::cout<<"read length is:"<<read_length<<std::endl;
    std::string fastq_line1, fastq_line2, fastq_line3, fastq_line4;
    std::string fastq_temp("abc");
    std::ifstream fastq(fastq_filename);
    std::ofstream fastq_read_length(fastq_filename+"_read_length.tsv");
    std::vector<int> read_length_count;
    for(int i=0;i!=read_length+1;++i)
    {
        read_length_count.push_back(0);
    }
    while(fastq>> fastq_line1 >> fastq_line2 >> fastq_line3 >> fastq_line4)
    {
        fastq_read_length<<fastq_line1<<'\t'<<fastq_line2.size()<<std::endl;
        ++read_length_count[fastq_line2.size()]; // Use read_length_count[50] to note how many 50bp reads
    }

    // std::cout<<"Count of reads with different lengths:"<<std::endl;
    std::ofstream fastq_read_length_count(fastq_filename+"_count_of_read_length.tsv");
    fastq_read_length_count<<"read_length"<<'\t'<<"total_number"<<std::endl;
    for(int i=0;i!=read_length_count.size();++i)
    {
        fastq_read_length_count<<i<<'\t'<<read_length_count[i]<<std::endl;
    }
    return 0;
}