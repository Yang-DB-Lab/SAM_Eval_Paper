#include <iostream>
#include <string>
#include <fstream>
#include <vector>

/***** Call the program using shell in this manner *********
*
CreateFastqRead1File chrom_info_file.tsv fastq_file.fq read_length adapter_seq adapter_length
*
************************************************************/
struct fastq_read_info //Each line of the sampled 1k chromosome tsv file!
{
    std::string chrom_num;
    std::string coordinate_start;
    std::string coordinate_end;
    std::string sequence;
};
// Overload operator >> for struct fastq_read_info
std::istream& operator>> (std::istream& is, fastq_read_info& dt)
{
    is>> dt.chrom_num>> dt.coordinate_start>> dt.coordinate_end>> dt.sequence;
    return is;
}

int main(int argc, char* argv[])
{
    struct fastq_read_info one_read_line;
    std::vector<fastq_read_info> all_read_info;

    // 1st command line parameter is the program name itself. Not useful.

    // 2nd command line parameter indicate file contains read_info
    std::ifstream chrom_info_file(argv[1]);
    // 3rd command line parameter for generated fastq_file location
    std::ofstream fastq_file(argv[2]);
    // 4th command line parameter is length of each read
    std::string read1_length(argv[3]);
    int read1_length_num = std::stoi(read1_length);
    // 5th command line parameter is adapter sequence
    std::string adapter_seq(argv[4]);
    // 6th command line parameter is adapter lenght in 3' end of reads
    std::string adapter_length=argv[5];
    int adapter_length_num = std::stoi(adapter_length);

    if(adapter_length_num > adapter_seq.size())
    {
        std::cout<<"ERROR: The adapter length parameter is longer than adapter sequence!"<<std::endl;
        return 0;
    }
    std::string read_name;
    std::string read_seq;
    std::string read_quality;
    while(chrom_info_file>>one_read_line)
    {
        read_name = '@' + one_read_line.chrom_num + '.' 
                    + one_read_line.coordinate_start + '-'
                    + std::to_string(std::stoi(one_read_line.coordinate_start) + read1_length_num -1) + '.'
                    + adapter_length +"bp.adapter";
        read_seq = one_read_line.sequence.substr(0, read1_length_num - adapter_length_num) + adapter_seq.substr(0, adapter_length_num);
        read_quality = std::string(read_seq.size(), 'J'); //faked highest score read quality
        fastq_file<<read_name<<std::endl
                  <<read_seq<<std::endl
                  <<'+'<<std::endl
                  <<read_quality<<std::endl;

    }
    
    chrom_info_file.close();
    fastq_file.close();

    return 0;
}
