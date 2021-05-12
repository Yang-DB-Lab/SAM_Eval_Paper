#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>


int main(int argc, char* argv[])
{
    // Note argv[0] is the name of the compiled program itself!!!
    std::string chrom_filename = argv[1];
    std::string sampled_chromfile = argv[2];
    std::string chromsome_number = argv[3];

    std::cout<<"chrom_file "<<chrom_filename<<std::endl;
    std::cout<<"Sampled_chromfile "<<sampled_chromfile<<std::endl;
    std::cout<<"chromosome_number "<<chromsome_number<<std::endl;

    std::ifstream chromosome(chrom_filename.c_str());
    std::ofstream sampled_file(sampled_chromfile.c_str());

    std::string chrom_line;
    std::string whole_chrom_seq;
    int first_nonN_base_idx= -1;
    int last_nonN_base_idx = -1;
    
    std::getline(chromosome, chrom_line);//get the first line which is chromosome annotation; NOT sequence
    while(std::getline(chromosome, chrom_line))
    {
        whole_chrom_seq+=chrom_line;
    }
    
    last_nonN_base_idx = whole_chrom_seq.size();

    for(int i=0;i!= whole_chrom_seq.size();++i)
    {
        if(whole_chrom_seq[i]=='N')
        {
            continue;
        }
        else
        {
            first_nonN_base_idx = i; //Get index of first real base(0 based index)
            break;
        }
        
    }
    std::cout<<"first_nonN_base_idx= "<<first_nonN_base_idx<<std::endl;
    std::cout<<"first nonN base is "<<whole_chrom_seq[first_nonN_base_idx]<<std::endl;

    for(int i= whole_chrom_seq.size() -1;i!= -1; --i)
    {
        if(whole_chrom_seq[i]=='N')
        {
            continue;
        }
        else
        {
            last_nonN_base_idx = i; //Get index of last real base(0 based index)
            break;
        }
        
    }
    std::cout<<"last_nonN_base_idx = "<<last_nonN_base_idx<<std::endl;
    std::cout<<"last_nonN_base is "<<whole_chrom_seq[last_nonN_base_idx]<<std::endl;
    // Generate random numbers which will be start of chromosome pieces.
    // Keep 0 based index
    std::vector<int> chrom_piece_start_indexes;
    std::srand(0); //set random number generator seed.
    for(int i=0;i!=whole_chrom_seq.size()/3000;++i)
    {
        int random_num; //  first_nonN_base_idx < random_num < last_nonN_base_idx - 1000
        random_num = std::rand()%(last_nonN_base_idx-1000 -first_nonN_base_idx) + first_nonN_base_idx;
        chrom_piece_start_indexes.push_back(random_num);
    }
    std::sort(chrom_piece_start_indexes.begin(), chrom_piece_start_indexes.end() );
    // remove duplicated coordinates, about 1/5000
    std::vector<int>::iterator vit;
    vit = std::unique(chrom_piece_start_indexes.begin(), chrom_piece_start_indexes.end() );
    chrom_piece_start_indexes.resize(std::distance(chrom_piece_start_indexes.begin(), vit) );

    // begin writing tsv file containing sampled 1kb genome pieces
    /* sampled_file<<"Chr"<<'\t'<<"Start"<<'\t'<<"End"<<'\t'<<"Sequence"<<std::endl; */
    // The coordinate is (string index + 1) (first index is 1)
    for(int i=0;i!= chrom_piece_start_indexes.size();++i)
    {
        /*
        sampled_file<<chromsome_number<<'\t'
                    <<chrom_piece_start_indexes[i]+1<<'\t'
                    <<chrom_piece_start_indexes[i]+1000<<'\t'
                    <<whole_chrom_seq.substr(chrom_piece_start_indexes[i], 1000)<<std::endl;
        */
        
        // Note some pieces may still have multiple Ns inside
        // Confirm no Ns inside before write to file
        std::string DNA_piece = whole_chrom_seq.substr(chrom_piece_start_indexes[i], 1000);
        std::string::size_type N_pos;
        N_pos = DNA_piece.find('N');
        if(N_pos == std::string::npos) // no 'N' in DNA_piece
        {
            sampled_file<<chromsome_number<<'\t'
                    <<chrom_piece_start_indexes[i]+1<<'\t'
                    <<chrom_piece_start_indexes[i]+1000<<'\t'
                    <<DNA_piece<<std::endl;
        }
    } 
    return 0;
}
