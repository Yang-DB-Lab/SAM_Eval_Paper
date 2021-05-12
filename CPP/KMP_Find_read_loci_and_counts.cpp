#include <iostream>
#include <string>
#include <vector>
#include <fstream>


void computeLPSArray(std::string& pat, int M, int* lps); 
  
// Prints occurrences of txt[] in pat[] 
size_t KMPSearch(std::string& pat, std::string& txt) 
{ 
    int M = pat.size(); 
    int N = txt.size(); 
  
    // create lps[] that will hold the longest prefix suffix 
    // values for pattern 
    int lps[M]; 
  
    // Preprocess the pattern (calculate lps[] array) 
    computeLPSArray(pat, M, lps); 
  
    int i = 0; // index for txt[] 
    int j = 0; // index for pat[] 
    while (i < N) { 
        if (pat[j] == txt[i]) { 
            j++; 
            i++; 
        } 
  
        if (j == M) { 
            // printf("Found pattern at index %d ", i - j); 
            return (i-j); 
            j = lps[j - 1]; 
        } 
  
        // mismatch after j matches 
        else if (i < N && pat[j] != txt[i]) { 
            // Do not match lps[0..lps[j-1]] characters, 
            // they will match anyway 
            if (j != 0) 
                j = lps[j - 1]; 
            else
                i = i + 1; 
        } 
    } 
    if(i==N && j!=M)
    {
        return std::string::npos;
    }
} 


size_t KMPSearch(std::string& pat, std::string& txt, int start_pos_txt)
{
    std::string shortered_txt = txt.substr(start_pos_txt, txt.size() - start_pos_txt);
    int found_location= KMPSearch(pat, shortered_txt);
    if(found_location != std::string::npos)
    {
        found_location =  found_location + start_pos_txt;
    }

    return found_location;
}
// Fills lps[] for given patttern pat[0..M-1] 
void computeLPSArray(std::string& pat, int M, int* lps) 
{ 
    // length of the previous longest prefix suffix 
    int len = 0; 
  
    lps[0] = 0; // lps[0] is always 0 
  
    // the loop calculates lps[i] for i = 1 to M-1 
    int i = 1; 
    while (i < M) { 
        if (pat[i] == pat[len]) { 
            len++; 
            lps[i] = len; 
            i++; 
        } 
        else // (pat[i] != pat[len]) 
        { 
            // This is tricky. Consider the example. 
            // AAACAAAA and i = 7. The idea is similar 
            // to search step. 
            if (len != 0) { 
                len = lps[len - 1]; 
  
                // Also, note that we do not increment 
                // i here 
            } 
            else // if (len == 0) 
            { 
                lps[i] = 0; 
                i++; 
            } 
        } 
    } 
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
        std::size_t found_pos= KMPSearch(fastq_read, all_chromsomes[i]); //all_chromsomes[i].find(fastq_read);
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
            found_pos = KMPSearch(fastq_read, all_chromsomes[i], found_pos+1); // all_chromsomes[i].find(fastq_read, found_pos+1);
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