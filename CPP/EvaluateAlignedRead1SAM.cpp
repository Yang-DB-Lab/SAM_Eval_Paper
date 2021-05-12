#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream> //for stringstream
#include <algorithm> // for std::min() function

// Get sequence of a whole chromosome
std::string read_in_chromosome_sequence(std::string chr_filename) // Only one chromosome!
{
	std::string chromosome_seq;
	std::ifstream chr_file(chr_filename.c_str());
	if(!chr_file)
		return 0;
	std::string chr_line;
	//get rid of first line which is chromosome information
	std::getline(chr_file, chr_line);
	// get all sequence line by line, but ignore line breaks.
	while(std::getline(chr_file, chr_line) )
	{
		chromosome_seq+=chr_line;
	}
	return chromosome_seq;
}

// NOTE: the sequence of reads in SAM file always follow the direction of reference genome, 
//       no matter this region contains a gene has the same direction or opposite direction.

struct SAM_entry
{
	std::string qname;
	int FLAG;
	std::string ref_seq_name; // chromosome name
	unsigned long int position; // 1-based left-most mapping position
	int map_quality;
	std::string CIGAR; // Annotation of alignment of read
	std::string ref_seq_next; //RNEXT 	String 	Ref. name of the mate/next read 
	int pos_next; //PNEXT 	Int 	Position of the mate/next read 
	int tem_length; //TLEN 	Int 	observed Template LENgth
	std::string read_seq; // sequence of the read
	std::string base_quality; // Phred base quality

	std::string detailed_annotation; // Tab seprated columns annotate current alignment

	// overloading "<<" operator to print SNP lines from SAM_file ifstream
	friend std::ostream & operator << (std::ostream & SAM_output,  SAM_entry & one_SAM_line)
	{
		SAM_output<<one_SAM_line.qname<<'\t';
		SAM_output<<one_SAM_line.FLAG<<'\t';
		SAM_output<<one_SAM_line.ref_seq_name<<'\t';
		SAM_output<<one_SAM_line.position<<'\t';
		SAM_output<<one_SAM_line.map_quality<<'\t';
		SAM_output<<one_SAM_line.CIGAR<<'\t';
		SAM_output<<one_SAM_line.ref_seq_next<<'\t';
		SAM_output<<one_SAM_line.pos_next<<'\t';
		SAM_output<<one_SAM_line.tem_length<<'\t';
		SAM_output<<one_SAM_line.read_seq<<'\t';
		SAM_output<<one_SAM_line.base_quality<<'\t';
		SAM_output<<one_SAM_line.detailed_annotation<<'\n';

		return SAM_output;
	}
	
};

std::vector<SAM_entry> get_whole_SAM(std::string SAM_filename)
{
	std::ifstream SAM_input(SAM_filename);
	std::string one_string_line;
	SAM_entry one_SAM_line;
	std::vector<SAM_entry> whole_SAM;
	
	while( std::getline(SAM_input, one_string_line) )
	{
		if(one_string_line[0] =='@') // header line, omit
			continue;
		else
		{
			std::stringstream ss(one_string_line);
			ss>>one_SAM_line.qname;
			ss>>one_SAM_line.FLAG;
			ss>>one_SAM_line.ref_seq_name;
			ss>>one_SAM_line.position;
			ss>>one_SAM_line.map_quality;
			ss>>one_SAM_line.CIGAR;
			ss>>one_SAM_line.ref_seq_next;
			ss>>one_SAM_line.pos_next;
			ss>>one_SAM_line.tem_length;
			ss>>one_SAM_line.read_seq;
			ss>>one_SAM_line.base_quality;

			std::getline(ss, one_SAM_line.detailed_annotation);

			whole_SAM.push_back(one_SAM_line);
		}
		
	}
	return whole_SAM;

}

// Get reverse complementary sequence of a DNA
std::string get_rev_com(std::string DNA)
{
	std::string rev_com_DNA;
	for(int i=0;i!=DNA.size();++i)
	{
		switch (DNA[i])
		{
		case 'A':
			rev_com_DNA = 'T' + rev_com_DNA;
			break;
		case 'T':
			rev_com_DNA = 'A' + rev_com_DNA;
			break;
		case 'C':
			rev_com_DNA = 'G' + rev_com_DNA;
			break;
		case 'G':
			rev_com_DNA = 'C' + rev_com_DNA;
			break;
		default:
			rev_com_DNA = 'N' + rev_com_DNA;
			break;
		}
	}
	return rev_com_DNA;
}

// Get total number of aligned bases between two sequences
// This program expect same number of bases in ref_DNA and read_DNA
float max_num_aligned_bases(std::string ref_DNA, std::string read_DNA)
{
	float read_DNA_aligned_num=0;
	for(int i=0;i!=read_DNA.size();++i)
	{
		if(ref_DNA[i] == read_DNA[i])
			read_DNA_aligned_num = read_DNA_aligned_num + 1.0;
	}
	//Check aligned bases if use rev_com sequence of read_DNA
	std::string read_rev_com = get_rev_com(read_DNA);
	float read_rev_com_aligned_num=0;
	for(int i=0;i!= read_DNA.size(); ++i)
	{
		if(ref_DNA[i] == read_rev_com[i])
			read_rev_com_aligned_num = read_rev_com_aligned_num + 1.0;
	}
	if(read_DNA_aligned_num > read_rev_com_aligned_num)
	{
		return (read_DNA_aligned_num+0.1);// With final .1 means the read itself provides better alignment
	}
	else if (read_DNA_aligned_num == read_rev_com_aligned_num)
	{
		return (read_DNA_aligned_num+0.2);// With final .2 means the read itself provides same alignment bases as its reverse complementary
	}
	else
	{
		return (read_rev_com_aligned_num + 0.3);// With final .1 means the reverse complemetary sequence of read provides better alignment
	}
	
}

/******************** NOTE of extract_CIGAR ***************
 * CIGAR is one column of SAM file like "51M2S","45M861N7M1S"
 * CIGAR string indicate how the raw read could align to reference genome.
 * The extract_CIGAR function get 
 *        a CIGAR sting, reference to a vector<int> and reference to a vector<char>
 * and split the numbers and symbol annotation into the two vectors.
 * Input CIGAR string 45M861N7M1S, will change numbers into 45 861 7 1
 *                                        anno_symbols into M N M S.
*/

// numbers and anno_symbols are extracted from CIGAR string
// numbers and anno_symbols should be defined first and then used as parameters for extract CIGAR.
void extract_CIGAR(std::string& anno_line, std::vector<int>& numbers, std::vector<char>& anno_symbols)
{
    int num_start_idx=0, num_end_idx=0; //num_start_idx: position of the 1st number in CIGAR string.
    int num; // A number in the CIGAR anno_line
    char sym; // A symbol in the CIGAR anno_line

    for(int i=0;i!=anno_line.size(); ++i)
    {
        if(anno_line[i]>='A' && anno_line[i]<='Z') // A char symbol found
        {
            num_end_idx = i-1;
            std::stringstream num_ss(anno_line.substr(num_start_idx, num_end_idx - num_start_idx +1));
            num_ss>>num;
            numbers.push_back(num);
            std::stringstream sym_ss(anno_line.substr(i,1));
            sym_ss>>sym;
            anno_symbols.push_back(sym);
            num_start_idx = i+1;
        }
    }
}

/********** NOTES on CIGAR code in SAM file *************
 * In SAM file(at least for the SAM file generated by RNA-STAR),
 * the CIGAR column has only 5 different annotations for raw reads,
 * they are:
 * M : This postion is matched(including incorrected-match)
 * S : Soft clip(At beginning or end of reads, indicated a piece can not align to reference)
 * N : NOT matched.(Indicate how long of the sequence skipped in the reference before the read can align again) 
 * D : Deletion(The read has deletions where can not align to reference)
 * I : Insertion(The read has insertion)
 * "M" and "S" are the most frequent annotation symbols.
*********************************************************/

std::string modi_read_to_match_ref(std::string& raw_read, std::vector<int>& CIGAR_numbers, std::vector<char>& anno_symbols)
{
    std::string modified_read;
    int cur_loci_in_read=0; 
    // Note where the checking is happening on raw read so that to cut out the right piece of sequence from raw read
    for(int i=0;i!=CIGAR_numbers.size(); ++i)
    {
        switch (anno_symbols[i])
        {
        case 'S': // soft clip, delete this part from raw read 
        {
            cur_loci_in_read += CIGAR_numbers[i]; //Move the current pointer forward. modified_read get nothing
            break;
        }
        case 'M': // Match part. Get the sequence into modified_read, move the pointer in raw read forward.
        {
            modified_read += raw_read.substr(cur_loci_in_read, CIGAR_numbers[i]);
            cur_loci_in_read += CIGAR_numbers[i]; 
            break;
        }
        case 'N': // Large skipping in genome. Add 'N's into modified_read. Do not move pointer in raw read.
        {
            modified_read += std::string( CIGAR_numbers[i], 'N');
            break;
        }
        case 'D': // Deletion in raw read, add space to fill the gap
        {
            modified_read += std::string( CIGAR_numbers[i], ' ');
            break;
        }
        case 'I': // Insertion in raw read. Delete the insertion: just move point in raw read forward.
        {
            cur_loci_in_read += CIGAR_numbers[i];
            break;
        }

        default:
            break;
        }
    }
    return modified_read;
}

// Get corespondig sequence in the reference of the modified_read( of raw read)
std::string get_matched_ref(std::string& chromosome, unsigned long int chrom_start, std::vector<int>& CIGAR_numbers, std::vector<char>& anno_symbols)
{
    std::string matched_ref_seq;
    int matched_ref_seq_length = 0;
    for(int i=0;i!=CIGAR_numbers.size(); ++i)
    {
        switch (anno_symbols[i])
        {
        case 'S': // soft clip, no corelated sequence in reference.
        {
            break;
        }
        case 'M': // Match part. Need get corelated sequence in reference.
        {
            matched_ref_seq_length += CIGAR_numbers[i]; 
            break;
        }
        case 'N': // Large skipping in genome. Get reference sequence related to this part.
        {
            matched_ref_seq_length += CIGAR_numbers[i];
            break;
        }
        case 'D': // Deletion in raw read, add space in modified_read and also get reference sequence.
        {
            matched_ref_seq_length += CIGAR_numbers[i];
            break;
        }
        case 'I': // Insertion in raw read. Just delete the insertion in raw read. No reference sequence.
        {
            break;
        }

        default:
            break;
        }
    }
    matched_ref_seq = chromosome.substr(chrom_start - 1, matched_ref_seq_length);
    return matched_ref_seq;

}

// Return chromosome sequence index in std::vector<std::string> all_chromosomes
// based on chromosome name
int get_chrom_idx(std::string ref_chrom_name)
{
	int chrom_idx=-1; // index of chromosome in "all_chromosomes"
	if(ref_chrom_name.size() == 1) //One character for the chromosome name: 1-9 or X or Y (or 0 if no alignment) 
	{
		if(ref_chrom_name[0]=='X')
		{
			chrom_idx = 19;
		}
		else if(ref_chrom_name[0]=='Y') //Use "else if", or the last "else" do not kown which if to pair!
		{
			chrom_idx = 20;
		}
		else
		{
			chrom_idx = std::stoi(ref_chrom_name) -1; //chromosome 1 has index 0;// chromosome 0 return -1.
		}
		
	}
	if(ref_chrom_name.size() == 2) // chromosome 10-19 or chromsome MT
	{
		if(ref_chrom_name[0] == 'M') //chromosome MT
		{
			chrom_idx = 21;
		}
		else
		{
			chrom_idx = std::stoi(ref_chrom_name) - 1; //chromosome 10 has index 10-1 ...
		}
		
	}
	return chrom_idx;
}

// qname has this format:
//                        chr10.3103041-3103090.2bp.adapter
bool extract_chrom_name_positions_from_qname(std::string& qname, std::string& chrom_name, unsigned long int& start, unsigned long int& end )
{
	int first_dot_idx;
	int second_dot_idx;
	int first_slash_idx;
	int i;
	for(i=0;i!= qname.size();++i)
	{
		if(qname[i]=='.')
		{
			first_dot_idx = i;
		}
		if(qname[i] == '-')
		{
			first_slash_idx = i;
			break;
		}
	}
	for(;i!=qname.size();++i)
	{
		if(qname[i] == '.')
		{
			second_dot_idx = i;
			break;
		}
	}
	chrom_name = qname.substr(3, first_dot_idx - 3); // chromosome name
	start = std::stoi(qname.substr(first_dot_idx+1, first_slash_idx-first_slash_idx-1) );
	end = std::stoi(qname.substr(first_slash_idx +1, second_dot_idx - first_slash_idx -1));

	if(i== qname.size())
		return false; //failed to extract
	else
	{
		return true;
	}
	

}

struct evaluation_line
{
	std::string qname;
	std::string read_seq;

	std::string theory_chromosome; // Where the read is generated!
	std::string theory_seq_full;
	// What is the number of aligned bases if align read to theratical genome string 
	// as defined in qname
	float aligned_bases_to_original_genome;

	std::string ref_chrom_in_SAM; // Align to which chromosome in SAM file
	std::string modi_ref_seq_in_SAM; // Align to this sequence in SAM file; modified to match alignment
	std::string modi_read_in_SAM; // modi to show alignment to modi_ref_seq
	// Actual aligned bases by align read to genome location defined in SAM
	float aligned_bases_to_ref_in_SAM;

	// overloading "<<" operator to print evaluation lines to ofstream(or other ostream)
	friend std::ostream & operator << (std::ostream & SAM_eval_output,  evaluation_line & one_eval_line)
	{
		SAM_eval_output<<one_eval_line.qname<<'\t';
		SAM_eval_output<<one_eval_line.read_seq<<'\t';
		SAM_eval_output<<one_eval_line.theory_chromosome<<'\t';
		SAM_eval_output<<one_eval_line.theory_seq_full<<'\t';
		SAM_eval_output<<one_eval_line.aligned_bases_to_original_genome<<'\t';
		SAM_eval_output<<one_eval_line.ref_chrom_in_SAM<<'\t';
		SAM_eval_output<<one_eval_line.modi_ref_seq_in_SAM<<'\t';
		SAM_eval_output<<one_eval_line.modi_read_in_SAM<<'\t';
		SAM_eval_output<<one_eval_line.aligned_bases_to_ref_in_SAM<<'\n';

		return SAM_eval_output;
	}
};

// While generating the fake fastq files with M bp genome sequence + N bp adapter,
// two other fastq files were also generated:
// (1) theory_clean_seq_of_reads: stores only the sequence of M bp genome sequence
// (2) theory_full_seq_of_reads:  stores (M+N) bp genome sequence start with start position of faked read
std::vector<std::string> get_theory_sequences_of_reads(std::string theory_seq_of_reads_filename)
{
	std::ifstream theory_seq_of_reads(theory_seq_of_reads_filename);
	std::vector<std::string> whole_theory_reads;
	std::string one_read_name;
	std::string one_theory_seq;
	std::string one_plus_symbol;
	std::string one_quality_line;

	while(std::getline(theory_seq_of_reads, one_read_name)&&
		  std::getline(theory_seq_of_reads, one_theory_seq)&&
		  std::getline(theory_seq_of_reads, one_plus_symbol)&&
		  std::getline(theory_seq_of_reads, one_quality_line) )
	{
		whole_theory_reads.push_back(one_theory_seq);
	}

	return whole_theory_reads;
}

std::vector<evaluation_line> get_evaluation_from_SAM(std::vector<SAM_entry>& partial_SAM,
												 std::vector<std::string>& all_chromosomes)
{
	std::vector<evaluation_line> evaluation_table;
	evaluation_line one_eval_line;
	int i;
	for(i=0;i!=partial_SAM.size(); ++i)
	{
		// evaluation_line one_eval_line = evaluation_table[i];
		std::string read_qname = partial_SAM[i].qname;
		one_eval_line.qname = read_qname;
		one_eval_line.read_seq = partial_SAM[i].read_seq;

		std::string read_chrom_name;
		unsigned long int read_start;
		unsigned long int read_end;
		extract_chrom_name_positions_from_qname(read_qname, read_chrom_name, read_start, read_end);

		// Read is generated from this chromosome:
		one_eval_line.theory_chromosome = read_chrom_name;
		

		int read_chrom_idx = get_chrom_idx(read_chrom_name);
		std::string read_chromosome = all_chromosomes[read_chrom_idx]; // Read is generated from this chromosome
		// If the whole read came from the chromsome (no adapter), the sequence should be this:
		std::string raw_read_chrom_seq = read_chromosome.substr(read_start-1, read_end-read_start); 
		one_eval_line.theory_seq_full = raw_read_chrom_seq;


		std::string read_ref_seq_name = partial_SAM[i].ref_seq_name;
		unsigned long int read_ref_position = partial_SAM[i].position;

		// Read is aligned to this chromosome in SAM
		one_eval_line.ref_chrom_in_SAM = read_ref_seq_name;

		// align to genome sequence which has same length as the read and also has same start coordiante.
		float max_num_aligned_bases_to_read_chrom= max_num_aligned_bases(raw_read_chrom_seq, partial_SAM[i].read_seq);
		one_eval_line.aligned_bases_to_original_genome = max_num_aligned_bases_to_read_chrom;
		// Alignment position is the same as faked fastq read.
		if(read_ref_seq_name == read_chrom_name && read_start == read_ref_position)
		{	
			one_eval_line.modi_ref_seq_in_SAM = raw_read_chrom_seq;
			one_eval_line.modi_read_in_SAM = partial_SAM[i].read_seq;
			// If reference position in SAM for a read is the same as theory position. Make the SAM align value smaller!
			one_eval_line.aligned_bases_to_ref_in_SAM = max_num_aligned_bases_to_read_chrom - 0.4;
			//evaluation_table[i] = one_eval_line;
			evaluation_table.push_back(one_eval_line);
			continue;
		}
		else
		{
			std::string CIGAR_anno = partial_SAM[i].CIGAR;
			std::vector<int> CIGAR_numbers;
			std::vector<char> CIGAR_symbols;
			extract_CIGAR(CIGAR_anno, CIGAR_numbers, CIGAR_symbols);

			int ref_chrom_idx = get_chrom_idx(read_ref_seq_name);

			std::string ref_chromosome = all_chromosomes[ref_chrom_idx]; // Aligned to this chromosome
			std::string raw_read_seq = partial_SAM[i].read_seq;

			std::string modified_read_matched = modi_read_to_match_ref(raw_read_seq, CIGAR_numbers, CIGAR_symbols);
			std::string matched_ref_seq = get_matched_ref(ref_chromosome, read_ref_position, CIGAR_numbers, CIGAR_symbols);
			one_eval_line.modi_ref_seq_in_SAM = matched_ref_seq;
			one_eval_line.modi_read_in_SAM = modified_read_matched;

			float read_max_aligned_base_to_SAM_ref_chrom = max_num_aligned_bases(matched_ref_seq, raw_read_seq);
			one_eval_line.aligned_bases_to_ref_in_SAM = read_max_aligned_base_to_SAM_ref_chrom;
			evaluation_table.push_back(one_eval_line);
		}
	}
	return evaluation_table;
}


/***************************************************************
argv[1] get SAM filename
// discarded: argv[2] get theory sequence fastq filename of all reads
argv[2] get SAM file evaluation table filename(tsv file)
****************************************************************/

int main(int argc, char* argv[])
{
	std::string SAM_filename=argv[1];
	// std::string theory_seq_of_reads_filename= argv[2];
	std::string evaluation_filename = argv[2];
   
	// std::cout<<"before get fq"<<std::endl;
	// std::vector<std::string> full_theory_seq_of_reads = get_theory_sequences_of_reads(theory_seq_of_reads_filename);
	// std::cout<<"after get fq"<<std::endl;

	std::vector<std::string> all_chromosomes; // Store sequence of all chromosomes.
	std::string chrom_temp;
	std::vector<std::string> chrom_filenames; // Store filenames of chromosome fasta files.
	
	for(int i=1;i!=20;++i)
	{
		chrom_filenames.push_back("MouseChromosomes/GRCm38.100.chromosome."+std::to_string(i) + ".fa");
	}
	chrom_filenames.push_back("MouseChromosomes/GRCm38.100.chromosome.X.fa");
	chrom_filenames.push_back("MouseChromosomes/GRCm38.100.chromosome.Y.fa");
	chrom_filenames.push_back("MouseChromosomes/GRCm38.100.chromosome.MT.fa");

	std::cout<<"chromsome_filename 5 "<<chrom_filenames[4]<<std::endl;

	for(int i=0;i!=chrom_filenames.size();++i)
	{
		chrom_temp=read_in_chromosome_sequence(chrom_filenames[i]);
		all_chromosomes.push_back(chrom_temp);
	}

	std::string DNA_Seq="AGCGTTTAA";
	std::string rev_com_DNA = get_rev_com(DNA_Seq);

	std::cout<<"Rev_com of "<<DNA_Seq<<" is "<<rev_com_DNA<<std::endl;

	std::cout<<"The 4000 to 4020 bases of chromsome 10 "<<all_chromosomes[9].substr(3999, 4020-4000)<<std::endl;



	// get whole SAM file
	std::vector<SAM_entry> whole_SAM = get_whole_SAM(SAM_filename);
	std::ofstream evaluation_file(evaluation_filename);
	evaluation_file<<"Read_Name"<<'\t';
	evaluation_file<<"Read_Seq"<<'\t';
	evaluation_file<<"Read_Generated_From"<<'\t';
	evaluation_file<<"Full_Theory_Seq_Without_Adapter"<<'\t';
	evaluation_file<<"Align_Score_with_full_Theory_Seq"<<'\t';
	evaluation_file<<"Read_Aligned_To(in_SAM)"<<'\t';
	evaluation_file<<"Modi_Ref_Seq_inSAM"<<'\t';
	evaluation_file<<"Modi_Read_ToAlign_inSAM"<<'\t';
	evaluation_file<<"Align_Score_based_on_SAM"<<std::endl;
	
	std::vector<evaluation_line> evaluation_table_from_SAM = get_evaluation_from_SAM(whole_SAM, all_chromosomes);

	for(int i=0;i!= evaluation_table_from_SAM.size(); ++i)
	{
		evaluation_file<< evaluation_table_from_SAM[i];
	}

	return 0;
    
}




