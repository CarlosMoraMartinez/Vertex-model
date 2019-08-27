#include "VertexSystem.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <random>

using namespace std;



int DEFAULT_STEPS = 10e7;
int DEFAULT_WRITE = 10e4;
int main(int argc, char *argv[]){

	int moves, print_step;
	if(argc < 4){
		moves = DEFAULT_STEPS;
		print_step = DEFAULT_WRITE;		
	}else{
		moves = std::stoi(argv[2]);
		print_step = std::stoi(argv[3]);
	}
	cout << "Reading from file\n\n"<<endl;
	std::string inputfile = argv[1];
	Tissue t = Tissue(inputfile, moves, print_step);
	try{
		t.simulate();
	}catch(const char* msg){
		t.produceOutputs("crashed");
		cout << msg << endl;
		exit(1);
	}
	cout << t.getStats() << endl;
	exit(0);
	
}
