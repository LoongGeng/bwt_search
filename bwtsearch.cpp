#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <map>
#include <algorithm>
#include <sstream>

#define block_size 2048 // max records read at once
#define max_check_point 126 // max ASCII code
#define check_point_offset 9 // min ASCII code

void build_idxfile(char *bwtfile, char *idxfile, int checkpoint[]); // build index file
void calculate_check_point(char *block, int checkpoint[]); // calculate the time of every symbol in a block read from bwtfile 
int Occ(char c, int num, std::ifstream &idx_in,std::ifstream &bwt_in); // calculate the time of an explicit symbol in bwtfile some where
int count_symbol(char c, int checkpoint[]); // count the number of symbols which has small ASCII code than an explicit symbol
void count_parttern(char *partten, int checkpoint[], char *idxfile, char *bwtfile, int pos[]); // count the number of pattern which we want to search
void remove_duplicate_record(char *idxfile, char *bwtfile, int first, int last, std::map<int, int> &count_out_du, int checkpoint[]); // remove patterns which have the same offset
int find_pre_in_bwtfile(std::ifstream &idx_in, std::ifstream &bwt_in, char *c, int pos_bwt, int checkpoint[]); // find the previous symbol of an explicit symbol
void find_offset(char *idxfile, char *bwtfile, std::map<int, int> &count_out_du, std::vector<unsigned int> &offset, int checkpoint[]); // find all offset of pattern which we want to search
int main (int argc, char* argv[]) {
    char *cmd_para = argv[1];
    char *bwtfile = argv[2];
    char *idxfile = argv[3];
    char *parttern = argv[4];
    int pos[2];
    int checkpoint[118] = {0};
    std::ifstream in;
    in.open(idxfile, std::ios_base::binary);
    //if index file does not exist, build the file
    //else calculate time of every symbol in the whole bwtfile
    if (!in.is_open()) {
        build_idxfile(bwtfile,idxfile,checkpoint);
    } else {
        int block_num;
        in.seekg(-(int)sizeof(int) ,std::ios_base::end);
        in.read((char *) &block_num, sizeof(int));
        //if block number is not zero
        //read the last record of index file
        if(block_num != 0) {
            in.seekg(-(int)((max_check_point - check_point_offset + 2)*sizeof(int)) ,std::ios_base::end);
            for (int i = 0; i <= (max_check_point - check_point_offset); i++)
                in.read((char *) &checkpoint[i], sizeof(int));
        }
        in.close();
        //read the last block of the bwt file
        in.open(bwtfile, std::ios_base::in);
        in.seekg(block_size * block_num);
        char *block = (char*) calloc(block_size, sizeof(char));
        in.read(block, block_size);
        if (strlen(block) < block_size && strlen(block) > 0)
            calculate_check_point(block, checkpoint);
        free(block);
        in.close();
    }
    //count the time of pattern
    count_parttern(parttern, checkpoint, idxfile, bwtfile, pos);
    if (strcmp(cmd_para, "-n") == 0) {
        std::cout<< pos[1] - pos[0] + 1<<std::endl;
        return 0;
    }else {
        if(pos[1] < pos[0]) {
            if (strcmp(cmd_para, "-r") == 0)
                std::cout << 0 << std::endl;
            return 0;
        }else {
        		//remove duplicate patterns which have the same offset in bwt file
                std::map<int, int> count_out_du;
                remove_duplicate_record(idxfile, bwtfile, pos[0], pos[1], count_out_du, checkpoint);
                if (strcmp(cmd_para, "-r") == 0) {
                    std::cout << count_out_du.size() << std::endl;
                    return 0;
                } else {
                    std::vector<unsigned int> offset;
                    //count offset for each pattern
                    find_offset(idxfile, bwtfile, count_out_du, offset, checkpoint);
                    std::sort(offset.begin(), offset.end());
                    std::vector<unsigned int>::iterator it;
                    for (it = offset.begin(); it != offset.end(); it++)
                        std::cout << "[" << *it << "]" << std::endl;
                    return 0;
                }
            }
    }

}

void build_idxfile(char *bwtfile, char *idxfile, int checkpoint[]) {
	//read a block from bwt file
	//count the time of each symbol in the block
	//write the check point into index file
    char *block = (char*) calloc(block_size, sizeof(char));
    std::ifstream in;
    std::ofstream out;
    int block_num = 0;
    in.open(bwtfile, std::ios_base::in);
    out.open(idxfile, std::ios_base::binary);
    in.read(block, block_size);
    while (!in.eof()) {
        block_num += 1;
        calculate_check_point(block, checkpoint);
        for (int i = 0; i <= (max_check_point - check_point_offset); i++)
            out.write((char *) &checkpoint[i], sizeof(int));
        free(block);
        block = (char *) calloc(block_size, sizeof(char));
        in.read(block, block_size);
    }
    out.write((char*) &block_num, sizeof(int));
    if (strlen(block) < block_size && strlen(block) > 0)
        calculate_check_point(block, checkpoint);
    free(block);
    in.close();
    out.close();
}

void calculate_check_point(char *block, int checkpoint[]) {
    for(int i = 0; i < strlen(block); i++)
        checkpoint[(int) block[i] - 9] += 1;
}

int count_symbol(char c, int checkpoint[]) {
	//count the number of symbols which has smaller ASCII code than an explicit symbol
    int pos = ((int) c)- check_point_offset;
    int count = 0;
    for (int i = 0; i < pos; i++)
        count += checkpoint[i];
    return count;
}

void count_parttern(char *partten, int checkpoint[], char *idxfile, char *bwtfile, int pos[]) {
    //count the time of pattern we want to search in bwt file
    //use FL transition and backward search 
    int i = strlen(partten) - 1;
    char c = partten[i];
    int first = count_symbol(c, checkpoint) + 1;
    int last = first -1 +checkpoint[(int) c - check_point_offset];
    std::ifstream idx_in, bwt_in;
    idx_in.open(idxfile, std::ios_base::binary);
    bwt_in.open(bwtfile, std::ios_base::in);
    //last > first means that the pattern does not exist
    while (first <= last && i > 0) {
        c = partten[--i];
        first = count_symbol(c, checkpoint) + Occ(c, first - 2, idx_in, bwt_in) + 1;
        last = count_symbol(c, checkpoint) + Occ(c, last - 1, idx_in, bwt_in);
    }
    pos[0] = first;
    pos[1] = last;
    idx_in.close();
    bwt_in.close();
}

int Occ(char c, int num, std::ifstream &idx_in,std::ifstream &bwt_in) {
	//count the time of an explicit symbol in bwt file from an explicit position
	//find the nearst index record of the explicit position
	//read the rest part from the recorded position to the explicit position in bwt file
	//count the time of the symbol in the resr part
	//return the sum of two counts
    int line = (num + 1) / block_size;
    int offset = (num + 1) % block_size;
    int sym_count = 0;
    int offest_count = 0;
    if (line != 0) {
        idx_in.seekg((118 * (line - 1) + (int) c - check_point_offset) * sizeof(int));
        idx_in.read((char *) &sym_count, sizeof(int));
    }
    char *buffer = (char*) calloc (offset, sizeof(char));
    bwt_in.seekg(line * block_size);
    bwt_in.read(buffer, offset);
    for (int i = 0; i < offset; i++) {
        if (buffer[i] == c)
            offest_count += 1;
    }
    free(buffer);
    return sym_count + offest_count;
}

int find_pre_in_bwtfile(std::ifstream &idx_in, std::ifstream &bwt_in, char *c, int pos_bwt, int checkpoint[]) {
	// FL transition
	// find previous symbol in bwt file
	// return the position of previous symbol in the first row
    bwt_in.seekg(pos_bwt);
    bwt_in.read(c, 1);
    return Occ(*c, pos_bwt,idx_in, bwt_in) + count_symbol(*c, checkpoint) - 1;
}

void remove_duplicate_record(char *idxfile, char *bwtfile, int first, int last, std::map<int, int> &count_out_du, int checkpoint[]) {
	// remove records which have the same offset in original file
	// find their nearst "]"
	// record the positions of "]"
	// remove duplicate records by using map
    std::ifstream idx_in, bwt_in;
    idx_in.open(idxfile, std::ios_base::in);
    bwt_in.open(bwtfile, std::ios_base::in);
    first = first -1;
    last = last -1;
    int pattern [last - first + 1] = {0};
    char *pre_sym = (char *) calloc(1, sizeof(char));
    std::string record_idx;
    for(int i = (last -first); i >= 0; i--) {
        if(pattern[i] == 0 ) {
            int pos_bwt = i + first;
            pos_bwt = find_pre_in_bwtfile(idx_in, bwt_in, pre_sym, pos_bwt, checkpoint);
            while ((*pre_sym) != ']') {
            	// when we seach the previous symbol, we may find the previous symbol is in the range where we find the pattern
            	// mark the positions in the range, they will not be visited in the following iteration
                if (pos_bwt >= first && pos_bwt <= last)
                    pattern[pos_bwt - first] = 1;
                pos_bwt = find_pre_in_bwtfile(idx_in, bwt_in, pre_sym, pos_bwt, checkpoint);
            }
            count_out_du.insert(std::map<int,int>::value_type(pos_bwt, 1));
        } else
            continue;
    }
    free(pre_sym);
    idx_in.close();
    bwt_in.close();
}

void find_offset(char *idxfile, char *bwtfile, std::map<int, int> &count_out_du, std::vector<unsigned int> &offset, int checkpoint[]) {
	// start from the position of "]"
	// find previous symbol by FL transition
	// if previous symbol is "[" then stop
	// the symbols between "[" and "]" is offset
    std::ifstream idx_in, bwt_in;
    idx_in.open(idxfile, std::ios_base::in);
    bwt_in.open(bwtfile, std::ios_base::in);
    std::map<int ,int>::iterator it;
    int pos_bwt;
    unsigned int offset_record;
    char *pre_sym = (char *) calloc(1, sizeof(char));
    std::stringstream ss;
    for (it = count_out_du.begin(); it!=count_out_du.end(); ++it) {
        std::string record = "";
        pos_bwt = find_pre_in_bwtfile(idx_in, bwt_in, pre_sym, it -> first, checkpoint);
        int j = 0;
        while((*pre_sym) != '[') {
            record = *pre_sym + record;
            pos_bwt = find_pre_in_bwtfile(idx_in, bwt_in, pre_sym, pos_bwt, checkpoint);
        }
        ss.clear();
        ss << record;
        ss >> offset_record;
        offset.push_back(offset_record);
    }
    free(pre_sym);
    idx_in.close();
    bwt_in.close();
}