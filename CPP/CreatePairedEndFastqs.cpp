#include <iostream>
#include <string>
#include <fstream>
#include <vector>

/***** Call the program using shell in this manner *********
*
CreateFastqRead1File chrom_info_file.tsv fastq_file.fq read_length adapter_seq adapter_length
*
************************************************************/

// Get reverse complimentary sequence of a DNA
std::string get_RevCom(std::string DNA)
{
    std::string DNA_RevCom;
    for(int i=0;i!=DNA.size();++i)
    {
        switch (DNA[i])
        {
        case 'A':
            DNA_RevCom= 'T' + DNA_RevCom;
            break;
        case 'T':
            DNA_RevCom= 'A' + DNA_RevCom;
            break;
        case 'C':
            DNA_RevCom= 'G' + DNA_RevCom;
            break;
        case 'G':
            DNA_RevCom= 'C' + DNA_RevCom;
            break;    
        default:
            DNA_RevCom= 'N' + DNA_RevCom;
            break;
        }
    }
    return DNA_RevCom;
}


// One line of the 1kb pieces tsv file
struct fastq_read_info
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
    // 3rd command line parameter for fastq filename (and location)
    std::string fastq_filename(argv[2]);
    std::ofstream fastq_read1(fastq_filename + ".R1.fq");
    std::ofstream fastq_read2(fastq_filename + ".R2.fq");

    // 4th command line parameter is length of each read
    std::string read1_length(argv[3]);
    int read1_length_num = std::stoi(read1_length);
    // 5th command line parameter is adapter sequence
    std::string adapter_seq(argv[4]);
    // 6th command line parameter is adapter lenght in 3' end of reads
    std::string adapter_length=argv[5];
    // 7th command line parameter is the gap between forward and reverse read
    std::string for_rev_gap(argv[6]);
    int gap=std::stoi(for_rev_gap);

    int adapter_length_num = std::stoi(adapter_length);

    if(adapter_length_num > adapter_seq.size())
    {
        std::cout<<"ERROR: The adapter length parameter is longer than adapter sequence!"<<std::endl;
        return 0;
    }
    std::string read_name;
    std::string read1_seq, read2_seq;
    std::string read_quality;
    while(chrom_info_file>>one_read_line)
    {
        read_name = '@' + one_read_line.chrom_num + '.' 
                    + one_read_line.coordinate_start + '-'
                    + std::to_string(std::stoi(one_read_line.coordinate_start) + read1_length_num -1) + '.'
                    + adapter_length +"bp.adapter"+ '.'
                    + for_rev_gap + "bp.gap";
        read1_seq = one_read_line.sequence.substr(0, read1_length_num - adapter_length_num) + adapter_seq.substr(0, adapter_length_num);
        read_quality = std::string(read1_seq.size(), 'J'); //faked highest score read quality. Forward + reverse reads.
        fastq_read1<<read_name<<std::endl
                  <<read1_seq<<std::endl
                  <<'+'<<std::endl
                  <<read_quality<<std::endl;

        // Get string of read2 on genome(rev_com of the real read2)
        read2_seq = one_read_line.sequence.substr(read1_length_num + gap, read1_length_num - adapter_length_num) +
                    adapter_seq.substr(0, adapter_length_num);
        // Real read2 seq after RevCom();
        read2_seq = get_RevCom(read2_seq);
        fastq_read2<<read_name<<std::endl
                  <<read2_seq<<std::endl
                  <<'+'<<std::endl
                  <<read_quality<<std::endl;

    }
    
    chrom_info_file.close();
    fastq_read1.close();
    fastq_read2.close();

    return 0;
}
