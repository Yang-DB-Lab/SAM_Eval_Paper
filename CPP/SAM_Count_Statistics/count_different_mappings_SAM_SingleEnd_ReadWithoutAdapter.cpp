/******************************************************
 * Purpose of this code:
 * 
 * Since all reads are generated from genome, all reads could map to the genome.
 * But some may be able to map to mutiple locations(nonUniq_mappint)
 * 
 * Get counts of different types of reads:
 * nonUniq_wrong: the read mapped to the wrong location of genome in SAM file.
 * nonUniq_right: the read mapped to the right location but is not unique.
 * Uniq_wrong: the read mapped to the right location and this mapping is unique. This is what we really want!
 * 
 * There is no "Uniq_wrong" because all reads can map to the right postion 
 * since the reads have no adapters(exactly the same as part of the reference genome!!)
 * ****************************************************/


#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream> //for stringstream
#include <algorithm> // for std::min() function

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

// Count how read1 and read2 of a sequencing piece is aligned: Uinq or nonUniq
void count_different_mapping(std::vector<SAM_entry>& SAM_data, 
							 std::ofstream& SAM_statistics_file)
{
    float SAM_size = SAM_data.size();
	// Single-End read, only Read1
	SAM_statistics_file<<"Read_Name"<<'\t'<<"Read1_Uniqueness"<<std::endl;


    for(int i=0;i!=SAM_data.size(); ++i)
    {
        std::string chrom_name;
        unsigned long int start;
        unsigned long int end;

		
		SAM_statistics_file<<SAM_data[i].qname<<'\t';

		// get theory chromosome name and start coordinate
		extract_chrom_name_positions_from_qname(SAM_data[i].qname, chrom_name, start, end);

			
		if(chrom_name == SAM_data[i].ref_seq_name && start == SAM_data[i].position)
		{//map to right position in SAM
			if(SAM_data[i].map_quality == 0) //multiple mapping(nonUniq)
			{//map_quality==0, sometimes can still be unique mapping. Check "XT:A:U" for unique mapping
				if(SAM_data[i].detailed_annotation.find("XT:A:U") != std::string::npos) //With "XT:A:U"
				{
					SAM_statistics_file<<"Uniq_right_0_score";
				}
				else
				{
					SAM_statistics_file<<"nonUniq_right";
				}
					
					
			}
			else //map_quality>0, must be uniq!
			{
				SAM_statistics_file<<"Uniq_right";
			}
				
		}
		else // map to wrong position in SAM(Must be nonUniq_wrong!!)
		{
			if(SAM_data[i].map_quality == 0) //multiple mapping(nonUniq)
			{//map_quality==0, sometimes can still be unique mapping. Check "XT:A:U" for unique mapping
				if(SAM_data[i].detailed_annotation.find("XT:A:U") != std::string::npos) //With "XT:A:U"
				{
					SAM_statistics_file<<"Uniq_wrong_0_score";
				}
				else
				{
					SAM_statistics_file<<"nonUniq_wrong";
				}
					
					
			}
			else //map_quality>0, must be uniq!
			{
				SAM_statistics_file<<"Uniq_wrong";
			}
				
		}
			SAM_statistics_file<<std::endl; // Prepare for read2
	}
}

/*************************************
how to call the program:
count_error_mapping_rate SAM_filename
Output:
SAM_filename.different_mapping_count.tsv: check whether each read is uniq or not.
 * **********************************/
int main(int argc, char* argv[])
{
    std::string SAM_filename = argv[1];

    std::vector<SAM_entry> whole_SAM_data = get_whole_SAM(SAM_filename);
    
    std::string output_tsv_filename = SAM_filename+ std::string(".different_mapping_count.tsv");
	std::ofstream output_tsv(output_tsv_filename);
	count_different_mapping(whole_SAM_data, output_tsv);

	
    output_tsv.close();

    return 0;
}