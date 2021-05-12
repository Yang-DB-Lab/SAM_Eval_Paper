#include <iostream>
#include <string>
#include <vector>
#include <fstream>

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



std::string read_in_chromosome(std::string chrom_filename)
{
    std::string chrom_seq;
    std::ifstream chrom_file(chrom_filename);
    if(!chrom_file)
    {
        std::cout<<"File open failed"<<std::endl;
        return 0;
    }

    std::string line;
    // First line is chromosome name, discard!
    std::getline(chrom_file, line);
    //std::cout<<"chrom_name: "<<line<<std::endl;
    while(std::getline(chrom_file, line) )
    {   /*
        if(line.back() <'A' || line.back() >'z' )
        {
            line.pop_back();
        }*/
        chrom_seq += line;
    }
    // std::cout<<"within read in chrom: chrom_seq.size()"<<chrom_seq.size()<<std::endl;
    chrom_file.close();
    return chrom_seq;
}

void read_in_fastq(std::string fastq_filename, std::vector<std::string>& fastq_read_names, std::vector<std::string>& fastq_read_seqs)
{
    std::string line1,line2,line3,line4;
    std::ifstream fastq_file(fastq_filename);
    while(std::getline(fastq_file, line1) &&
          std::getline(fastq_file, line2) &&
          std::getline(fastq_file, line3) &&
          std::getline(fastq_file, line4))
    {
        fastq_read_names.push_back(line1);
        fastq_read_seqs.push_back(line2);
    }
    fastq_file.close();

}

// get location of a fastq read(without adapter sequences) in the genome.
// If more than one location exists, get only two loci.
std::vector<std::string> get_matched_loci_of_read_in_genome(std::string& fastq_read, std::vector<std::string>& all_chromsomes)
{
    std::vector<std::string> read_loci;
    for(int i=0;i!= all_chromsomes.size();++i)
    {
        std::size_t found_pos=all_chromsomes[i].find(fastq_read);
        if(found_pos!=std::string::npos)
        {
            std::string found_location;
            if(i+1 <= 19) // Autochromosomes
            {
                found_location = "chr" + std::to_string(i+1) + ":" +
                                 std::to_string(found_pos+1) + "-" + std::to_string(found_pos+1+fastq_read.size() );
                read_loci.push_back(found_location);
            }
            if(i+1 == 20)
            {
                found_location = "chrX:" + 
                                 std::to_string(found_pos+1) + "-" + std::to_string(found_pos+1+fastq_read.size() );
                read_loci.push_back(found_location);
            }
            if(i+1 == 21)
            {
                found_location = "chrY:" + 
                                 std::to_string(found_pos+1) + "-" + std::to_string(found_pos+1+fastq_read.size() );
                read_loci.push_back(found_location);
            }
            if(i+1 == 22)
            {
                found_location = "chrMT:" + 
                                 std::to_string(found_pos+1) + "-" + std::to_string(found_pos+1+fastq_read.size() );
                read_loci.push_back(found_location);
            }
            if(read_loci.size() > 1)// get multiple loci already
            {
                break;
            }

            // 2nd trial to find the position of the read in this chromosome 
            found_pos = all_chromsomes[i].find(fastq_read, found_pos+1);
            if(found_pos != std::string::npos)
            {
                std::string found_location;
                if(i+1 <= 19) // Autochromosomes
                {
                    found_location = "chr" + std::to_string(i+1) + ":" +
                                    std::to_string(found_pos+1) + "-" + std::to_string(found_pos+1+fastq_read.size() );
                    read_loci.push_back(found_location);
                }
                if(i+1 == 20)
                {
                    found_location = "chrX:" + 
                                    std::to_string(found_pos+1) + "-" + std::to_string(found_pos+1+fastq_read.size() );
                    read_loci.push_back(found_location);
                }
                if(i+1 == 21)
                {
                    found_location = "chrY:" + 
                                    std::to_string(found_pos+1) + "-" + std::to_string(found_pos+1+fastq_read.size() );
                    read_loci.push_back(found_location);
                }
                if(i+1 == 22)
                {
                    found_location = "chrMT:" + 
                                    std::to_string(found_pos+1) + "-" + std::to_string(found_pos+1+fastq_read.size() );
                    read_loci.push_back(found_location);
                }
                if(read_loci.size() > 1)// get multiple loci already
                {
                    break;
                }
            }
        }
        
    }
    // check reverse complimentary sequence of the read in genome
    if(read_loci.size() <=1)
    {
        std::string read_rev_com = get_RevCom(fastq_read);
        for(int i=0;i!= all_chromsomes.size();++i)
        {
            std::size_t found_pos=all_chromsomes[i].find(read_rev_com);
            if(found_pos!=std::string::npos) //Found in genome
            {
                std::string found_location;
                if(i+1 <= 19) // Autochromosomes
                {
                    found_location = "RevCom_chr" + std::to_string(i+1) + ":" +
                                    std::to_string(found_pos+1) + "-" + std::to_string(found_pos+1+fastq_read.size() );
                    read_loci.push_back(found_location);
                }
                if(i+1 == 20)
                {
                    found_location = "RevCom_chrX:" + 
                                    std::to_string(found_pos+1) + "-" + std::to_string(found_pos+1+fastq_read.size() );
                    read_loci.push_back(found_location);
                }
                if(i+1 == 21)
                {
                    found_location = "RevCom_chrY:" + 
                                    std::to_string(found_pos+1) + "-" + std::to_string(found_pos+1+fastq_read.size() );
                    read_loci.push_back(found_location);
                }
                if(i+1 == 22)
                {
                    found_location = "RevCom_chrMT:" + 
                                    std::to_string(found_pos+1) + "-" + std::to_string(found_pos+1+fastq_read.size() );
                    read_loci.push_back(found_location);
                }
                break;
            }
        }
    }
    return read_loci;
}




int main(int argc, char* argv[])
{
    //Only one command line parameter, it is the fastq file
    std::string fastq_filename(argv[1]);

    std::vector<std::string> all_chromsomes;
    std::vector<std::string> chromosome_filenames;
    std::string chrom_file_dir("/data/guang/Aligner_Benchmark_Paper/MouseChromosomes/GRCm38.100.chromosome.");
    for(int i=1;i<=19;++i)
    {
        chromosome_filenames.push_back(chrom_file_dir + std::to_string(i) + std::string(".fa") );
    }
    chromosome_filenames.push_back(chrom_file_dir + "X.fa" );
    chromosome_filenames.push_back(chrom_file_dir + std::string("Y.fa") );
    chromosome_filenames.push_back(chrom_file_dir + std::string("MT.fa") );

    //Read in all chromosomes  
    for(int i=0;i!=chromosome_filenames.size();++i)
    {
        all_chromsomes.push_back(read_in_chromosome(chromosome_filenames[i]) );
        std::cout<<"all_chromosomes"<<i<<"size"<<all_chromsomes[i].size()<<std::endl;
    }

    // Import all sequences in fastq file
    std::vector<std::string> fastq_read_names;
    std::vector<std::string> fastq_read_seqs;
    read_in_fastq(fastq_filename, fastq_read_names, fastq_read_seqs);
    std::cout<<"fastq_read_seqs.size()"<<fastq_read_seqs.size()<<std::endl;

    std::ofstream count_read_num_in_genome(fastq_filename + ".read.counts.loci.tsv");
    for(int i=0;i!= fastq_read_seqs.size(); ++i)
    {
        std::vector<std::string> read_loci= get_matched_loci_of_read_in_genome(fastq_read_seqs[i], all_chromsomes);
        if(read_loci.size()==0)
        {
            count_read_num_in_genome<<fastq_read_names[i]<<'\t';
            count_read_num_in_genome<<read_loci.size()<<'\t';
            count_read_num_in_genome<<'\t'<<std::endl;
        }

        if(read_loci.size()==1)
        {
            count_read_num_in_genome<<fastq_read_names[i]<<'\t';
            count_read_num_in_genome<<read_loci.size()<<'\t';
            count_read_num_in_genome<<read_loci[0]<<std::endl;
        }
        if(read_loci.size()==2)
        {
            count_read_num_in_genome<<fastq_read_names[i]<<'\t';
            count_read_num_in_genome<<read_loci.size()<<'\t';
            count_read_num_in_genome<<(read_loci[0]+"_"+read_loci[1])<<std::endl;
        }
    }

    return 0;
}