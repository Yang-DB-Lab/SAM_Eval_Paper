#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

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
    // Note argv[0] is the name of the compiled program itself!!!
    std::string chrom_filename = argv[1]; // Where is the sequence file of one chromosome
    std::string sampled_tsv_filename = argv[2]; // The sampled tsv file of this chromosome which includes start coordinate
    int read_length = std::stoi(argv[3]); // How long you will need in the longer_read(>=1kb) fastq file.
    std::string chrom_fastq_filename = argv[4]; // Fastq file related to this chromosome 

    std::cout<<"chrom_file "<<chrom_filename<<std::endl;
    std::cout<<"Sampled_tsv_file with coordinates "<<sampled_tsv_filename<<std::endl;
    std::cout<<"read_length "<<read_length<<std::endl;
    std::cout<<"chrom_fastq_filname: "<<chrom_fastq_filename<<std::endl;

    std::ifstream chromosome(chrom_filename.c_str());
    std::ifstream sampled_tsv_file(sampled_tsv_filename.c_str());
    std::ofstream chrom_fastq_file(chrom_fastq_filename);

    std::string chrom_line;
    std::string whole_chrom_seq;
 //   int first_nonN_base_idx= -1;
 //   int last_nonN_base_idx = -1;
    
    std::getline(chromosome, chrom_line);//get the first line which is chromosome annotation; NOT sequence
    while(std::getline(chromosome, chrom_line))
    {
        whole_chrom_seq+=chrom_line;
    }
    chromosome.close(); // close the chromosome file containing sequence of the whole chromosome

    // One line for the sampled tsv file
    struct fastq_read_info one_read_line;
 
    std::string read_name;
    std::string read_seq;
    std::string read_quality;
    while(sampled_tsv_file>>one_read_line)
    {
        
        read_seq = whole_chrom_seq.substr(std::stoi(one_read_line.coordinate_start) -1, read_length); // If shorter than read_length, get as long as possible
        std::string::size_type N_pos = read_seq.find('N');
        if(N_pos != std::string::npos) // 'N' found in substring sequence from whole_chrom_seq with "read_length" length
        {
            read_seq = read_seq.substr(0,N_pos); //Cut out the "N" bases in the 3' end.
        }

        read_name = '@' + one_read_line.chrom_num + '.' 
                    + one_read_line.coordinate_start + '-'
                    + std::to_string(std::stoi(one_read_line.coordinate_start) + read_seq.size() -1) + '.'
                    + "0bp.adapter";

        read_quality = std::string(read_seq.size(), 'J'); //faked highest score read quality
        chrom_fastq_file<<read_name<<std::endl
                  <<read_seq<<std::endl
                  <<'+'<<std::endl
                  <<read_quality<<std::endl;

    }
    
    sampled_tsv_file.close();
    chrom_fastq_file.close();

    return 0;
}
