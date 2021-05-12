#include <iostream>
#include <vector>
#include <algorithm>

int main()
{
    std::vector<int> chrom_piece_start_indexes={1,2,3,3,4,2,5,6,9,7,5,0,9};
    std::cout<<"content of my chrom_piece_start_indexes vector: "<<std::endl;
    for(int i=0;i!=chrom_piece_start_indexes.size();++i)
    {
        std::cout<<chrom_piece_start_indexes[i]<<" ";
    }
    std::cout<<std::endl;

    std::sort(chrom_piece_start_indexes.begin(), chrom_piece_start_indexes.end() );
    // remove duplicated coordinates, about 1/5000
    std::vector<int>::iterator vit;
    vit = std::unique(chrom_piece_start_indexes.begin(), chrom_piece_start_indexes.end() );
    chrom_piece_start_indexes.resize(std::distance(chrom_piece_start_indexes.begin(), vit) );

    std::cout<<"content of my chrom_piece_start_indexes vector after sort and remove duplicates: "<<std::endl;
    for(int i=0;i!=chrom_piece_start_indexes.size();++i)
    {
        std::cout<<chrom_piece_start_indexes[i]<<" ";
    }
    std::cout<<std::endl;
}