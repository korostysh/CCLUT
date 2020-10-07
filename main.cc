#include "CCLUT.h"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <bitset>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <vector>

using namespace cw;

int main(int argc, char * argv[]) {
/*
	//unsigned int uint = 120;
	//unsigned char a , b = 255;

	gemPacket pack;

	//std::cout << std::hex << uint << std::endl;
	std::cout << std::hex << pack << std::endl;

	std::fstream ofs("test_file.pat", std::ios_base::out);
	for (int i = 0; i < 512; i++) {
		ofs << pack;
	}
	ofs.close();

	return 0;

 UNCOMMENT TO TEST PATTERN GEN
*/
  std::srand(std::time(0));
	if(argc <= 1) {
        	std::cout << "Error! No filename specified!" << std::endl;
        	return -1;
	}

	std::string file_name(argv[1]);
	std::vector<CCLUT> ccluts;
	std::vector<Cluster> gem;

	std::cout << "Reading from File:  " << file_name << std::endl;
	ReadTxt(file_name, ccluts, gem);
	std::cout << "Successfully read " << ccluts.size() << " Muons" << std::endl;
	std::cout << "Successfully read " << gem.size()   << " GEM clusters" << std::endl;

	for(int i=0;i<ccluts.size();i++){
		PrintCCLUT(ccluts.at(i));
	}

	std::cout << std::endl;

	std::cout << "Generating " << file_name << "(.pat) files" << std::endl;

	if(WritePat(file_name, ccluts)){
		return 0;
	}
	else{
		return -1;
	}

/*

	std::string pref = "SamplePattern";
	std::string title = "TestRPattern";
	std::vector<CLCT> pat;
	std::vector<Hit> raw;

	int N = ReadTxt(pref, pat);
	ExtractHits(pat, raw);

	std::cout << pat.size() << std::endl;
	std::cout << raw.size() << std::endl;

	std::cout << std::endl;

	for(int i=0; i < pat.size(); i++){
		PrintCLCT(pat.at(i));
	}

	WriteTxt(title, pat);

	bool ek = WritePat(pat, title);

	char ch = 'a';
	while (ch == 'a')
	{
		std::cout << "Press Any Key to Quit:	";
		std::cin >> ch;
	}

*/


	return 0;
}
