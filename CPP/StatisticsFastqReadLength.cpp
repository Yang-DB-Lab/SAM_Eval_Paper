#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility> //for std::pair
#include <algorithm> //for std::sort

/******************************************
Call program in this way:
StatisticsFastqReadLength fastq_filename statistic_tsv_filename
 * ****************************************/
int main(int argc, char* argv[])
{
    std::string fastq_filename(argv[1]);
    std::ifstream fastq(fastq_filename);
    std::string statistics_filename(argv[2]);
    std::ofstream statistics(statistics_filename);

    std::string one_read_name;
	std::string one_read_seq;
	std::string one_plus_symbol;
	std::string one_quality_line;

    std::vector<int> all_read_length;
    std::pair<int,int> one_readLength_count;
    std::vector<std::pair<int,int> > all_readLength_counts;

	while(std::getline(fastq, one_read_name)&&
		  std::getline(fastq, one_read_seq)&&
		  std::getline(fastq, one_plus_symbol)&&
		  std::getline(fastq, one_quality_line) )
    {
        while(one_read_seq.back() < 'A') // get rid of special characters in the end of read if exist.
        {
            one_read_seq.pop_back();
        }
        int one_read_length=one_read_seq.size();
        all_read_length.push_back(one_read_length);
    }
    std::sort(all_read_length.begin(), all_read_length.end() );
    std::cout<<"How many reads get? "<<all_read_length.size()<<std::endl;
    
    one_readLength_count = std::make_pair(all_read_length[0],1);//shortest read length: count 1 time.
    for(int i=1; i < all_read_length.size();++i)
    {
        if(all_read_length[i]!=one_readLength_count.first)// new read length found
        {
            all_readLength_counts.push_back(one_readLength_count); // get old read length count
            one_readLength_count = std::make_pair(all_read_length[i], 1); // create 1 count for a new read length
        }
        else
        {
            one_readLength_count.second = one_readLength_count.second +1; // Same read length met, make it count + 1
        }
        
    }
    // get last readlengthCount
    all_readLength_counts.push_back(one_readLength_count);

    for(int i=0;i!= all_readLength_counts.size();++i)
    {
        statistics<<all_readLength_counts[i].first<<'\t'<<all_readLength_counts[i].second<<std::endl;
    }

    return 0;

}
