/******************************************************
 * Purpose of this code:
 * 
 * Get counts of different types of reads:
 * Uniq_wrong: the read mapped to a wrong location, but uniquely.
 * nonUniq_nonRepairable_wrong: the read mapped to a wrong location in SAM. Meanwhile, although this mapping is not unique,
 * 							     the right location is not one of the optimal mapping locations.
 * nonUniq_Repairable_wrong: the read mapped to a wrong location in SAM. Meanwhile there are more than one optimal
 * 							 mapping locations(nonUniq) and the right location is one of the optimal locations.
 * nonUniq_right: the read mapped to the right location but is not unique.
 * Uniq_right: the read mapped to the right location and this mapping is unique. This is what we really want!
 * ****************************************************/


#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream> //for stringstream
#include <algorithm> // for std::min() function


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

// Local alignment function. Here only get possible max score, 
// do not note how to get the acutal alignment produce the score 
//Referrence
//http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Smith-Waterman

struct score_well //Scoring system for alignment
{
	int pre_x;
	int pre_y;
	int current_score;
	score_well()
	{
		pre_x=0;
		pre_y=0;
		current_score=0;
	}
	score_well(int a, int b, int c)
	{
		pre_x=a;
		pre_y=b;
		current_score=c;
	}
};

float max_score_Smith_Waterman_local_alignment(std::string DNA_Seq1, std::string DNA_Seq2, 
                                               std::string align_algorithm, 
											   std::vector<std::vector<score_well> >& matrix_scored) //DNA_Seq1 is reference
{
		
	DNA_Seq1 = '-' + DNA_Seq1;
	DNA_Seq2 = '-' + DNA_Seq2;
	//std::vector<std::vector<score_well> > matrix_scored;   //Matrix_scored: DNA_Seq1.size()+1 rows, DNA_Seq2.size()+1 columns. //DNA_Seq1 uses columns; DNA_Seq2 uses rows
	int match_award = 1; // ai == bj, get 1 point!
	int gap_open_penalty= -6; // value from bwa; for case of s(ai,-) or s(-,bj).   //Add '-' to head of DNA_Seq1 and DNA_Seq2 to marker only gap no symbol move case. The left up corner means (-,-);
	int gap_extending_penalty = -3; //value from bwa
	int mismatch_penalty = -3; //value from bwa

	int mem_match_award = 1;
	int mem_gap_open_penalty = -6;
	int mem_gap_extend_penalty = -1;
	int mem_mismatch_penalty = -4;
	int mem_clip_penalty = -5;

	int aln_match_award = 1;
	int aln_gap_open_penalty = -11;
	int aln_gap_extend_penalty = -4;
	int aln_mismatch_penalty = -3;

	if(align_algorithm=="aln") //bwa aln
	{
		match_award = aln_match_award;
		gap_open_penalty = aln_gap_open_penalty;
		gap_extending_penalty = aln_gap_extend_penalty;
		mismatch_penalty = aln_mismatch_penalty;
	}
	else
	{
		match_award = mem_match_award;
		gap_open_penalty = mem_gap_open_penalty;
		gap_extending_penalty = mem_gap_extend_penalty;
		mismatch_penalty = mem_mismatch_penalty;
	}
	
	score_well init_well(-1,-1,0); //pre_x = -1, pre_y=-1,score=0; The upleft corner have no previous coordinate. label as (-1,-1);
	
	for(int i=0; i!= DNA_Seq1.size();++i) 
	{
		std::vector<score_well> temp_sw_vec;
		for(int j=0;j!=DNA_Seq2.size();++j)
		{
			temp_sw_vec.push_back(init_well);
		}
		matrix_scored.push_back(temp_sw_vec);
	}
	//No need to set first row and first column, they all have score 0(Manually assigned, not from previous well)

	float max_score = 0; // Now start score is 0;

	for(int i=1;i!=matrix_scored.size();++i)
	{
		for(int j=1;j!=matrix_scored[0].size();++j) //now matrix_scored[i][j].current_score ==0 after the first initialization.
		{                                           //No assign 0 choice!!! Must choose from s(ai,bj),s(ai,-) or s(-,bj)
			
			int s_ai_bj = mismatch_penalty; //ai != bj; Mis-match
			if(DNA_Seq1[i]==DNA_Seq2[j]) //ai == bj ; Match
				s_ai_bj = match_award;
			//move right down from i-1 j-1 position: s(ai,bj); First choice.

			matrix_scored[i][j].current_score = matrix_scored[i-1][j-1].current_score + s_ai_bj;
			matrix_scored[i][j].pre_x = i-1;
			matrix_scored[i][j].pre_y = j-1;
			
			int gap_penalty_in_use = gap_open_penalty;// set gap_open_penalty as default
			if(matrix_scored[i-1][j].pre_y == j) // The next operation gap extention
			{
				gap_penalty_in_use = gap_extending_penalty;
			}
			if(matrix_scored[i-1][j].current_score + gap_penalty_in_use > matrix_scored[i][j].current_score)  //s(ai,-) // Move down
			{
				matrix_scored[i][j].current_score = matrix_scored[i-1][j].current_score + gap_penalty_in_use;
				matrix_scored[i][j].pre_x = i-1;
				matrix_scored[i][j].pre_y = j;
			}


			gap_penalty_in_use = gap_open_penalty;
			if(matrix_scored[i][j-1].pre_x == i) // The next operation gap extention
			{
				gap_penalty_in_use = gap_extending_penalty;
			}
			if(matrix_scored[i][j-1].current_score + gap_penalty_in_use > matrix_scored[i][j].current_score)  //s(-,bj)
			{
				matrix_scored[i][j].current_score = matrix_scored[i][j-1].current_score + gap_penalty_in_use;
				matrix_scored[i][j].pre_x = i;
				matrix_scored[i][j].pre_y = j-1;
			}
			if(matrix_scored[i][j].current_score < 0) // If current_score<0, assign to 0.
			{
				matrix_scored[i][j].current_score = 0;
				matrix_scored[i][j].pre_x = -1;
				matrix_scored[i][j].pre_y = -1;
				
			}
			if(matrix_scored[i][j].current_score > max_score)
				max_score = matrix_scored[i][j].current_score;
		}
	}
	return max_score;

}

score_well get_max_scored_well(std::vector<std::vector<score_well> > scored_matrix) //Get value of the highest score well; With the value: (pre_x+1, pre_y+1) is the indexes of this well.
{//for local alignment
	score_well max_scored_well;

	for(int i=0;i!=scored_matrix.size();++i)
	{
		for(int j=0;j!=scored_matrix[0].size();++j)
		{
			if(max_scored_well.current_score <=scored_matrix[i][j].current_score)
			{
				max_scored_well = scored_matrix[i][j];
			}
		}
	}
	return max_scored_well;
}
std::vector<std::string> get_longest_path_of_local_alignment(std::vector<std::vector<score_well> > scored_matrix, std::string DNA_Seq1, std::string DNA_Seq2)
{//local aligned_results     //Even the two sequence can not align at all. A path will also be returned. This path may not be the longest local fully alignment
	std::string Seq1_aligned, Seq2_aligned;
	std::string Seq1_left_piece_remained, Seq1_right_piece_remained, Seq2_left_piece_remained, Seq2_right_piece_remained; // Also return sequence pieces do not exist in the path.
	score_well max_scored_well = get_max_scored_well(scored_matrix);
	
	if(max_scored_well.current_score==0)//So no alignment between two DNA sequences!! Just align them together(align to the right edge) and no remain pieces left!
	{
		if(DNA_Seq1.size() == 0) //First sequence is empty
		{
			if(DNA_Seq2.size() >0)                      //return:
			{                                           //         - - ...-
				for(int i=0;i!=DNA_Seq2.size();++i)         //         b1b2...bn  
				{
					Seq1_aligned ='-' + Seq1_aligned;
				}
				Seq2_aligned = DNA_Seq2;
			} //else: both empty; do nothing.
		}
		else if(DNA_Seq2.size()==0)
		{
			if(DNA_Seq1.size() > 0) // NO need to do this if clause because of previous 'else'.
			{
				for(int i=0;i!=DNA_Seq1.size();++i)    //return:
				{                                      //        a1a2 ... an
					Seq2_aligned ='-' + Seq2_aligned;  //         - - ... -
				}
				Seq1_aligned = DNA_Seq1;
			}
		}
		else //Both DNA_Seq1 and DNA_Seq2 are not empty
		{
			if(DNA_Seq1.size() > DNA_Seq2.size())
			{                                                           //return in form:
				for(int i=0;i!=(DNA_Seq1.size() - DNA_Seq2.size());++i) //      a1a2a3a4a5
				{                                                       //      - - - b1b2
					Seq2_aligned = '-' + Seq2_aligned;                  //
				}
			}
			else //DNA_Seq2 is longer
			{                                                            //return in form:
				for(int i=0;i!=(DNA_Seq2.size() - DNA_Seq1.size() );++i) //      - - - a1a2
				{                                                        //      b1b2b3b4b5  
					Seq1_aligned = '-' + Seq1_aligned;
				}
			}
			Seq1_aligned = Seq1_aligned + DNA_Seq1;
			Seq2_aligned = Seq2_aligned + DNA_Seq2;
		}
		std::vector<std::string> aligned_seqs_and_remainings; //Return structure: Seq1_aligned, Seq2_aligned, Seq1_left_piece_remained, Seq1_right_piece_remained, Seq2_left_piece_remained, Seq2_right_piece_remained
		aligned_seqs_and_remainings.push_back(Seq1_aligned);
		aligned_seqs_and_remainings.push_back(Seq2_aligned);
		aligned_seqs_and_remainings.push_back(Seq1_left_piece_remained);
		aligned_seqs_and_remainings.push_back(Seq1_right_piece_remained);
		aligned_seqs_and_remainings.push_back(Seq2_left_piece_remained);
		aligned_seqs_and_remainings.push_back(Seq2_right_piece_remained);
	    
		return aligned_seqs_and_remainings;
	}
	//This is the main part:The score > 0. The else part!
	
	int x= max_scored_well.pre_x + 1, y= max_scored_well.pre_y +1; //Get indexes of biggest scored well.//This well must move right and down from previous well.
	
	DNA_Seq1 = '-' + DNA_Seq1;
	DNA_Seq2 = '-' + DNA_Seq2;
	
	Seq1_right_piece_remained = DNA_Seq1.substr(x+1, DNA_Seq1.size() - (x+1) );
	Seq2_right_piece_remained = DNA_Seq2.substr(y+1, DNA_Seq2.size() - (y+1) );
	while( scored_matrix[x][y].current_score > 0 ) 
	{
		if( scored_matrix[x][y].pre_x == x ) //Moved right to get here!!
		{
			Seq1_aligned = '-' + Seq1_aligned;
			Seq2_aligned = DNA_Seq2[y] + Seq2_aligned;
			//std::cout<<"-"<<std::endl<<DNA_Seq2[y]<<std::endl;
			
		}
		else if( scored_matrix[x][y].pre_y == y )//Moved down to get here!!
		{
			Seq1_aligned = DNA_Seq1[x] + Seq1_aligned;
			Seq2_aligned = '-' + Seq2_aligned;
			//std::cout<<DNA_Seq1[x]<<std::endl<<"-"<<std::endl;
		}
		else
		{
			Seq1_aligned = DNA_Seq1[x] + Seq1_aligned;
			Seq2_aligned = DNA_Seq2[y] + Seq2_aligned;
			//std::cout<<DNA_Seq1[x]<<std::endl<<DNA_Seq2[y]<<std::endl;
		}
		
		
		double current_x=x;
		double current_y=y;
		x= scored_matrix[current_x][current_y].pre_x;
		y= scored_matrix[current_x][current_y].pre_y;
		
		//AA-TGCATCG
        //AACATCAT-G
	}
	Seq1_left_piece_remained = DNA_Seq1.substr(1,x); //Note DNA_Seq1[0] now is '-'
	Seq2_left_piece_remained = DNA_Seq2.substr(1,y);
	std::vector<std::string> aligned_seqs_and_remainings;
	aligned_seqs_and_remainings.push_back(Seq1_aligned);
	aligned_seqs_and_remainings.push_back(Seq2_aligned);
	aligned_seqs_and_remainings.push_back(Seq1_left_piece_remained);
	aligned_seqs_and_remainings.push_back(Seq1_right_piece_remained);
	aligned_seqs_and_remainings.push_back(Seq2_left_piece_remained);
	aligned_seqs_and_remainings.push_back(Seq2_right_piece_remained);
	
	return aligned_seqs_and_remainings;
}

// The right alignment method to use

std::vector<std::vector<score_well> > Needleman_Wunsch_global_alignment(std::string DNA_Seq1, std::string DNA_Seq2, std::string align_algorithm)
{
	DNA_Seq1 = '-' + DNA_Seq1;
	DNA_Seq2 = '-' + DNA_Seq2;

	//std::vector<std::vector<score_well> > matrix_scored;   //Matrix_scored: DNA_Seq1.size()+1 rows, DNA_Seq2.size()+1 columns. //DNA_Seq1 uses columns; DNA_Seq2 uses rows
	int match_award = 1; // ai == bj, get 1 point!
	int gap_open_penalty= -6; // value from bwa; for case of s(ai,-) or s(-,bj).   //Add '-' to head of DNA_Seq1 and DNA_Seq2 to marker only gap no symbol move case. The left up corner means (-,-);
	int gap_extending_penalty = -3; //value from bwa
	int mismatch_penalty = -3; //value from bwa

	int mem_match_award = 1;
	int mem_gap_open_penalty = -6;
	int mem_gap_extend_penalty = -1;
	int mem_mismatch_penalty = -4;
	int mem_clip_penalty = -5;

	int aln_match_award = 1;
	int aln_gap_open_penalty = -11;
	int aln_gap_extend_penalty = -4;
	int aln_mismatch_penalty = -3;

	if(align_algorithm=="aln") //bwa aln
	{
		match_award = aln_match_award;
		gap_open_penalty = aln_gap_open_penalty;
		gap_extending_penalty = aln_gap_extend_penalty;
		mismatch_penalty = aln_mismatch_penalty;
	}
	else
	{
		match_award = mem_match_award;
		gap_open_penalty = mem_gap_open_penalty;
		gap_extending_penalty = mem_gap_extend_penalty;
		mismatch_penalty = mem_mismatch_penalty;
	}
	int gap_penalty = gap_open_penalty; //default vaule as open penalty

	std::vector<std::vector<score_well> > matrix_scored;   //Matrix_scored: DNA_Seq1.size()+1 rows, DNA_Seq2.size()+1 columns. //DNA_Seq1 uses columns; DNA_Seq2 uses rows
	
	score_well init_well(-1,-1,0); //pre_x = -1, pre_y=-1,score=0; The upleft corner have no previous coordinate. label as (-1,-1);
	
	for(double i=0; i!= DNA_Seq1.size();++i) 
	{
		std::vector<score_well> temp_sw_vec;
		for(double j=0;j!=DNA_Seq2.size();++j)
		{
			temp_sw_vec.push_back(init_well);
		}
		matrix_scored.push_back(temp_sw_vec);
	}
	
	for(double i=1;i!=matrix_scored.size();++i) //First column: can only move down!
	{
		matrix_scored[i][0].pre_x= i-1;
		matrix_scored[i][0].pre_y= 0;
		matrix_scored[i][0].current_score= gap_open_penalty + (i-1)*gap_extending_penalty;
	}
	for(double j=1;j!=matrix_scored[0].size();++j) //First row: can only move right!
	{
		matrix_scored[0][j].pre_x = 0;
		matrix_scored[0][j].pre_y = j-1;
		matrix_scored[0][j].current_score=gap_open_penalty + (j-1)*gap_extending_penalty;
	}
	
	for(double i=1;i!=matrix_scored.size();++i)
	{
		for(double j=1;j!=matrix_scored[0].size();++j) //now matrix_scored[i][j].current_score ==0 after the first initialization.
		{                                           //No assign 0 choice!!! Must choose from s(ai,bj),s(ai,-) or s(-,bj)
			
			double s_ai_bj= mismatch_penalty; //ai != bj; Mis-match
			if(DNA_Seq1[i]==DNA_Seq2[j]) //ai == bj ; Match
				s_ai_bj = match_award;
			//move right down from i-1 j-1 position: s(ai,bj); First choice.
			matrix_scored[i][j].current_score = matrix_scored[i-1][j-1].current_score + s_ai_bj;
			matrix_scored[i][j].pre_x = i-1;
			matrix_scored[i][j].pre_y = j-1;
			if(matrix_scored[i-1][j].pre_y == j) // The next operation gap extention
			{
				gap_penalty = gap_extending_penalty;
			}
			if(matrix_scored[i-1][j].current_score + gap_penalty > matrix_scored[i][j].current_score)  //s(ai,-) // Move down
			{
				matrix_scored[i][j].current_score = matrix_scored[i-1][j].current_score + gap_penalty;
				matrix_scored[i][j].pre_x = i-1;
				matrix_scored[i][j].pre_y = j;
			}

			gap_penalty = gap_open_penalty;
			if(matrix_scored[i][j-1].pre_x == i) // The next operation gap extention
			{
				gap_penalty = gap_extending_penalty;
			}
			if(matrix_scored[i][j-1].current_score + gap_penalty > matrix_scored[i][j].current_score)  //s(-,bj)
			{
				matrix_scored[i][j].current_score = matrix_scored[i][j-1].current_score + gap_penalty;
				matrix_scored[i][j].pre_x = i;
				matrix_scored[i][j].pre_y = j-1;
			}
		}
	}
	return matrix_scored;
}

int max_global_mapping_score(std::string DNA_Seq1, std::string DNA_Seq2, std::string align_algorithm)
{
	int soft_clipping_penalty=0; 
	
	//soft_clipping_penalty = mismatch_penalty-1;
	if(align_algorithm =="aln")
	{
		soft_clipping_penalty = -4;
	}
	else
	{
		soft_clipping_penalty = -5;
	}
	
	std::vector<std::vector<score_well> > scored_matrix = Needleman_Wunsch_global_alignment(DNA_Seq1, DNA_Seq2, align_algorithm);
	//score_well max_score_well = get_max_scored_well(scored_matrix);
	// Get max_score_well without function calling.
	score_well max_scored_well;
	// When the alignment rich last character of DNA_Seq2(the read), what is the max score;
	int max_last_j_score = -10000;
	int max_j, j;
	for(int i=0;i!=scored_matrix.size();++i)
	{
		for(j=0;j!=scored_matrix[0].size();++j)
		{
			if(max_scored_well.current_score <=scored_matrix[i][j].current_score)
			{
				max_scored_well = scored_matrix[i][j];
				max_j = j;
			}
		}
		if(scored_matrix[i][j-1].current_score > max_last_j_score)
		{
			max_last_j_score = scored_matrix[i][j-1].current_score;
		}
	}

	// set the final score as align until last character of DNA_seq2(the read)
	int max_final_score=max_last_j_score;

	// Consider soft-clipping penalty as DNA-Seq2(the read) does not reach end of its sequence.
	if(max_j < DNA_Seq2.size() )
	{
		int max_score_well_minus_clipping_penalty = max_scored_well.current_score + (DNA_Seq2.size() - max_j) * soft_clipping_penalty;
		if(max_score_well_minus_clipping_penalty > max_final_score)
		{
			max_final_score = max_score_well_minus_clipping_penalty;
		}
	}
	return max_final_score;
}

std::vector<std::string> get_global_alignment_result(std::vector<std::vector<score_well> > scored_matrix, std::string DNA_Seq1, std::string DNA_Seq2)
{//local aligned_results
	std::string Seq1_aligned, Seq2_aligned;
	DNA_Seq1 = '-' + DNA_Seq1;
	DNA_Seq2 = '-' + DNA_Seq2;
	double x=scored_matrix.size()-1, y=scored_matrix[0].size()-1; //largest index of scored_matrix
	
	
	while( !(x==0 && y==0)) //x==0, y==0 the left up corner
	{
		if( scored_matrix[x][y].pre_x == x ) //Moved right to get here!!
		{
			Seq1_aligned = '-' + Seq1_aligned;
			Seq2_aligned = DNA_Seq2[y] + Seq2_aligned;
			//std::cout<<"-"<<std::endl<<DNA_Seq2[y]<<std::endl;
			
		}
		else if( scored_matrix[x][y].pre_y == y )//Moved down to get here!!
		{
			Seq1_aligned = DNA_Seq1[x] + Seq1_aligned;
			Seq2_aligned = '-' + Seq2_aligned;
			//std::cout<<DNA_Seq1[x]<<std::endl<<"-"<<std::endl;
		}
		else
		{
			Seq1_aligned = DNA_Seq1[x] + Seq1_aligned;
			Seq2_aligned = DNA_Seq2[y] + Seq2_aligned;
			//std::cout<<DNA_Seq1[x]<<std::endl<<DNA_Seq2[y]<<std::endl;
		}
		
		
		double current_x=x;
		double current_y=y;
		x= scored_matrix[current_x][current_y].pre_x;
		y= scored_matrix[current_x][current_y].pre_y;
		
	}
	std::vector<std::string> aligned_seqs;
	aligned_seqs.push_back(Seq1_aligned);
	aligned_seqs.push_back(Seq2_aligned);
	
	return aligned_seqs;
	
}

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

			// for SAM line can not align to genome, there is no annotation
			one_SAM_line.detailed_annotation = "";
			if(ss)
			{
				std::getline(ss, one_SAM_line.detailed_annotation);
			}
			
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

// qname has this format:
//                        chr10.3103041-3103090.2bp.adapter_theorySequence
bool extract_read_adapter_length_from_qname(std::string& qname, 							  
								  			int& read_length,
								  			int& adapter_length) 
{
	// std::string& chrom_name, 
	unsigned long int start; 
	unsigned long int end;

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
	//chrom_name = qname.substr(3, first_dot_idx - 3); // chromosome name
	start = std::stoi(qname.substr(first_dot_idx+1, first_slash_idx-first_slash_idx-1) );
	end = std::stoi(qname.substr(first_slash_idx +1, second_dot_idx - first_slash_idx -1));
	
	if(start < end)
	{
		read_length = end - start + 1;
	}
	else
	{
		read_length = start - end + 1;
	}
	
	
	//get adapter length
	std::string adapter_length_str;
	int j=second_dot_idx+1;
	while( qname[j]!='b') //get everything before 'b' which is part of "bp.adapter"
	{
		adapter_length_str = adapter_length_str + qname[j];
		++j;
	}
	if(adapter_length_str.size()>0)
	{
		adapter_length = std::stoi(adapter_length_str);
	}

	if(i== qname.size())
		return false; //failed to extract
	else
	{
		return true;
	}
	
}
// chromesome_name:  1, 2, 3, ... , 17, 18, 19, X, Y, MT
// chromesome_index: 0, 1, 2, ... , 16, 17, 18, 19,20,21
int get_chromesome_index(std::string chrom_name)
{
	if(chrom_name == "X") //1, 2, 3, 4, 5, 6, 7, 8, 9, X, Y
	{
		return 19;
	}
	else if(chrom_name == "Y")
	{
		return 20;
	}
	else if(chrom_name == "MT")
	{
		return 21;
	}
	else if(chrom_name[0]>='0' && chrom_name[0]<='9')
	{
		int chrom_name_integar = std::stoi(chrom_name); //NOTE the no mapping reads!!!!
		return chrom_name_integar-1;
	}
	else
	{
		return -1;
	}
	
	
}

// Use chromesome_index, start_coordinate and length to get corresponding sequence from genome
std::string extract_chrom_seq_using_coordinate_length(std::vector<std::string>& all_chromesomes, 
													  int chrom_index, unsigned int start, unsigned int length) 
{
	std::string seq_wanted = all_chromesomes[chrom_index].substr(start-1, length);
	return seq_wanted;
}

int extract_optimal_mapping_score_from_sam_annotation(std::string& SAM_line_detailed_annotation)
{
	int mapping_score= 0;
	// mapping score has this format: MD:Z:30. This is part of SAM_line last part: detailed_annotation

	// Find "Z" from the end(since "MD:Z:" is at the relative end part of annotation)
	// NOTE: sometimes one read may have "XA:Z" column(62300th row of GRCm38.100.genome.reads.50bpWith_3_bp_Adapter.bwa_aln.sam)
	// so only check "Z" is not safe, must also check "MD" exists before 'Z'!!
	int Z_inMDZ_index=0; 
	for(int i= SAM_line_detailed_annotation.size()-1; i>=0; --i)
	{
		if(SAM_line_detailed_annotation[i]=='Z') // found 'Z'
		{
			if(SAM_line_detailed_annotation[i-2] =='D' && SAM_line_detailed_annotation[i-3]=='M')//MD:Z found!!
			{
				Z_inMDZ_index = i;
			}
		}

	}
	int score_idx= Z_inMDZ_index + 2;
	std::string score_str;
	while(SAM_line_detailed_annotation[score_idx]>='0' && SAM_line_detailed_annotation[score_idx] <= '9' 
		  && score_idx < SAM_line_detailed_annotation.size() )
	{
		score_str = score_str + SAM_line_detailed_annotation[score_idx];
		++ score_idx;
	}

	mapping_score = std::stoi(score_str);

	return mapping_score;
}

float get_error_mapping_rate(std::vector<SAM_entry>& SAM_data)
{
    float SAM_size = SAM_data.size();

    float error_mapping_count = 0;
    for(int i=0;i!=SAM_data.size(); ++i)
    {
        std::string chrom_name;
        unsigned long int start;
        unsigned long int end;;

        // get theory chromosome name and start coordinate
        extract_chrom_name_positions_from_qname(SAM_data[i].qname, chrom_name, start, end);

        if(chrom_name == SAM_data[i].ref_seq_name && start == SAM_data[i].position)
        {
            continue;
        }
        else
        {
            ++ error_mapping_count;
        }
        
    }
    return (error_mapping_count/SAM_size);
}


void count_different_mapping(std::vector<SAM_entry>& SAM_data, 
							 std::vector<std::string>& all_chromosomes,//Contains seqs of all chromosomes
							 std::string full_adapter_seq, // Long sequence used to generate adater_sequence
							 std::string align_algorithm, //What bwa algorithm was used to get SAM file("aln" or "mem")
							 std::ofstream& SAM_statistics_file) //tsv file to store results of each read
{
	int read_length, adapter_length; // get read_length and adapter_length of this SAM file
	                                 // by checking the 0st SAM read qname.
	extract_read_adapter_length_from_qname(SAM_data[0].qname, read_length, adapter_length);
	// Get adapter sequence
    std::string adapter_seq = full_adapter_seq.substr(0, adapter_length);

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
		if(SAM_data[i].ref_seq_name == "*") // Align to chromosome "*" in SAM file, which means no alignment.
		{
			SAM_statistics_file<<"No_Mapping";
		}
			
		else // Map to someplace
		{
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
					else //nonUniq_wrong
					{
						std::string seq_used_for_mapping_in_sam = SAM_data[i].read_seq;//This sequence is used to count mapping quality
						std::string chrom_name_in_sam = SAM_data[i].ref_seq_name;//The read map to this chromosome
						int chrom_index_in_sam = get_chromesome_index(chrom_name_in_sam);
						unsigned int mapping_position_in_sam = SAM_data[i].position;
						//Use a longer ref_seq to calculate mapping score for the read in SAM.
						//Its total length=read_length+adapter_length=read_without_adapter + 2*adapter_length
						std::string longer_ref_seq_in_sam = extract_chrom_seq_using_coordinate_length(all_chromosomes,
																									chrom_index_in_sam,
																									mapping_position_in_sam,
																									read_length+adapter_length);
						// For debug
						//std::cout<<"longer_ref_seq_in_sam: "<<longer_ref_seq_in_sam<<std::endl;
						//std::cout<<"seq_used_for_mapping_in_sam: "<<seq_used_for_mapping_in_sam<<std::endl;

						int max_mapping_score_of_read_in_sam = max_global_mapping_score(longer_ref_seq_in_sam,  
																						seq_used_for_mapping_in_sam,
																						align_algorithm);

						int index_theory_chrom = get_chromesome_index(chrom_name); // The sequence is generated here
						std::string read_seq_in_fastq= extract_chrom_seq_using_coordinate_length(all_chromosomes,
																								index_theory_chrom,
																								start,
																								read_length - adapter_length) + adapter_seq;
						std::string longer_theory_ref_seq = extract_chrom_seq_using_coordinate_length(all_chromosomes,
																									index_theory_chrom,
																									start, read_length + adapter_length);
						// For debug
						//std::cout<<"read_seq_in_fastq: "<<read_seq_in_fastq<<std::endl;
						//std::cout<<"longer_theory_ref_seq: "<<longer_theory_ref_seq<<std::endl;

						int max_mapping_score_theory = max_global_mapping_score(longer_theory_ref_seq, read_seq_in_fastq, align_algorithm);
						
						// For debug
						//std::cout<<"max sam score: "<<max_mapping_score_of_read_in_sam <<"max theory score: "<<max_mapping_score_theory<<std::endl;

						if(max_mapping_score_theory < max_mapping_score_of_read_in_sam)
						{ // The sequence align to wrong place better!!
							SAM_statistics_file<<"nonUniq_wrong_nonRepairable";
						} 
						else // Also check whether read mapped into reverse complimentary strand.
						{
							/*int max_mapping_score_of_read_in_sam_RevCom = max_global_mapping_score(longer_ref_seq_in_sam,  
																				get_RevCom(seq_used_for_mapping_in_sam),
																						align_algorithm);
							if(max_mapping_score_of_read_in_sam < max_mapping_score_of_read_in_sam_RevCom)
							{
								max_mapping_score_of_read_in_sam = max_mapping_score_of_read_in_sam_RevCom;
							}
							if(max_mapping_score_theory < max_mapping_score_of_read_in_sam)
							{ // The sequence align to wrong place better!!
								SAM_statistics_file<<"nonUniq_wrong_nonRepairable";
							} 
							else*/ if(max_mapping_score_theory == max_mapping_score_of_read_in_sam)
							{
								SAM_statistics_file<<"nonUniq_wrong_Repairable";
							}
							else
							{
								SAM_statistics_file<<"RightPlaceBetterScoreButWrongMapped";
							}
						}
						
						
						
						
					}
						
						
				}
				else //map_quality>0, must be uniq!
				{
					SAM_statistics_file<<"Uniq_wrong";
				}
					
			}
		}
		SAM_statistics_file<<std::endl; // Finshed one read, change to next line for next read.
	}
}

/*
void SumUp_different_mapping(std::vector<SAM_entry>& SAM_data, 
							 std::vector<std::string>& all_chromosomes,//Contains seqs of all chromosomes
							 //int read_length,                  // Length of read(including adapter)
							 //int adapter_length,               // Length of 3' adapter in the read
							 int& Uniq_right,                  // The read mapped to the right location and this mapping is unique. 
							                                   // This is what we really want!
							 int& Uniq_right_0_score,          // For "bwa-aln" sam only. Some reads have "XT:A:U" label but 0 mapq.
							 int& nonUniq_right,               // The read mapped to the right location but is not unique.
							 int& Uniq_wrong,                  // The read mapped to a wrong location, but uniquely.
							 int& Uniq_wrong_0_score,          // Also determinated wrong!
							 int& nonUniq_nonRepairable_wrong, // The read mapped to a wrong location in SAM. 
							                                   // Meanwhile, although this mapping is not unique,
															   // the right location is not one of the optimal mapping locations.
							 int& nonUniq_Repairable_wrong)    // The read mapped to a wrong location in SAM. 
							                                   // Meanwhile there are more than one optimal mapping locations(nonUniq)
															   // and the right location is one of the optimal locations.
							 
							 
{
    float SAM_size = SAM_data.size();

    Uniq_wrong = 0;
	nonUniq_nonRepairable_wrong = 0; 
	nonUniq_Repairable_wrong = 0;
	nonUniq_right = 0;
	Uniq_right = 0;

	int read_length, adapter_length; // get read_length and adapter_length of this SAM file
	                                 // by checking the 0st SAM read qname.
	extract_read_adapter_length_from_qname(SAM_data[0].qname, read_length, adapter_length);

    for(int i=0;i!=SAM_data.size(); ++i)
    {
        std::string chrom_name;
        unsigned long int start;
        unsigned long int end;;

        // get theory chromosome name and start coordinate
        extract_chrom_name_positions_from_qname(SAM_data[i].qname, chrom_name, start, end);

        if(chrom_name == SAM_data[i].ref_seq_name && start == SAM_data[i].position)
        {//map to right position in SAM
            if(SAM_data[i].map_quality == 0) //multiple mapping(nonUniq)
			{
				++ nonUniq_right;
			}
			else
			{
				++ Uniq_right;
			}
			
        }
        else // map to wrong position in SAM(Must be nonUniq_wrong!!)
        {
            if(SAM_data[i].map_quality == 0) // nonUniq wrong
			{
				int mapping_score = extract_optimal_mapping_score_from_sam_annotation(SAM_data[i].detailed_annotation);
				if(mapping_score <= read_length - adapter_length) // Adapter sequence can not map to reference
				{                                                 // The same read can also map to its right position
					++ nonUniq_Repairable_wrong; 
				}
				else
				{
					++ nonUniq_nonRepairable_wrong;
				}
				

			}
			else //Wrong and unique position
			{
				++ Uniq_wrong;
			}
			
        }
        
    }
}
*/

/*************************************
how to call the program:
count_error_mapping_rate SAM_filename full_adapter_sequence align_algorithm
Output:
SAM_filename.tsv(which has only one line: SAM_filename  error_mapping_rate)
 * **********************************/


int main(int argc, char* argv[])
{
    std::string SAM_filename = argv[1];
	std::string full_adapter_seq = argv[2];
	std::string align_algorithm = argv[3];

    std::vector<SAM_entry> whole_SAM_data = get_whole_SAM(SAM_filename);

	std::vector<std::string> chrom_filenames;
	for(int i=1;i!=20;++i)
	{
		chrom_filenames.push_back("/data/guang/SAM_evaluation_Paper/MouseChromosomes/GRCm38.100.chromosome."+std::to_string(i) + ".fa");
	}
	chrom_filenames.push_back("/data/guang/SAM_evaluation_Paper/MouseChromosomes/GRCm38.100.chromosome.X.fa");
	chrom_filenames.push_back("/data/guang/SAM_evaluation_Paper/MouseChromosomes/GRCm38.100.chromosome.Y.fa");
	chrom_filenames.push_back("/data/guang/SAM_evaluation_Paper/MouseChromosomes/GRCm38.100.chromosome.MT.fa");

	std::cout<<"chromsome_filename 5 "<<chrom_filenames[4]<<std::endl;

    std::vector<std::string> all_chromosomes;
	for(int i=0;i!=chrom_filenames.size();++i)
	{
		std::string chrom_temp=read_in_chromosome_sequence(chrom_filenames[i]);
		all_chromosomes.push_back(chrom_temp);
	}
    
	//std::string output_tsv_filename = SAM_filename+ std::string(".count_different_mapping.tsv");
    //std::ofstream output_tsv(output_tsv_filename);
	std::ofstream SAM_statistics_file(SAM_filename + ".different_mapping_count.tsv");
	
	count_different_mapping(whole_SAM_data, 
							all_chromosomes,//Contains seqs of all chromosomes
							full_adapter_seq, // Long sequence used to generate adater_sequence
							align_algorithm, //What bwa algorithm was used to get SAM file("aln" or "mem")
							SAM_statistics_file);

	SAM_statistics_file.close();
    return 0;
}

/*
int main()
{
	std::string DNA1="ATCGAGTGCGAGGACGGTTTCAATTTT";
	std::string DNA2="ATGGAGCTAGTCGACTTTTTGCAAA";
	std::string align_method="mem";
	std::vector<std::vector<score_well> > scored_matrix;
	float score= max_score_Smith_Waterman_local_alignment(DNA1, DNA2, align_method, scored_matrix);

	std::cout<<"Align score: "<<score<<std::endl;


	std::vector<std::string> local_alignment_results = get_longest_path_of_local_alignment(scored_matrix, DNA1,DNA2);
	std::cout<<"Local Aligned Results:"<<std::endl;
	for(double i=0;i!=local_alignment_results.size();++i)
	{
		std::cout<<local_alignment_results[i]<<std::endl;
	}

	std::vector<std::vector<score_well> > global_scored_matrix = Needleman_Wunsch_global_alignment(DNA1,DNA2,align_method);
	std::vector<std::string> global_align_result = get_global_alignment_result(global_scored_matrix, DNA1, DNA2);
	std::cout<<"Global Aligned Results:"<<std::endl;
	for(double i=0;i!=global_align_result.size();++i)
	{
		std::cout<<global_align_result[i]<<std::endl;
	}
	int i= global_scored_matrix.size();
	int global_align_score = global_scored_matrix[i-1].back().current_score;
	std::cout<<"Global alignment final score= "<<global_align_score<<std::endl;

	score_well max_global_align_score = get_max_scored_well(global_scored_matrix);
	std::cout<<"Max global alignment score of along the alignment: "<<max_global_align_score.current_score<<std::endl;
	return 0;
} 
*/