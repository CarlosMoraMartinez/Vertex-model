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
	if(argc < 5){
		moves = DEFAULT_STEPS;
		print_step = DEFAULT_WRITE;		
	}else{
		moves = std::stoi(argv[3]);
		print_step = std::stoi(argv[4]);
	}
	cout << "Reading arguments"<<endl;
	std::string inputfile = argv[1];
	std::string paramfile = argv[2];
	cout << "Reading from file\n\n"<<endl;
	Tissue t = Tissue(inputfile, paramfile, moves, print_step, inputfile);

	srand(RANDOM_SEED);
	std::default_random_engine generator;
	std::uniform_real_distribution<double> unif;

	try{
		t.simulate(generator, unif);
		cout << t.getStats() << endl; 
		//t.addAcceptedMovements(10);
		//t.simulate(generator, unif);
		//cout << t.getStats() << endl; 
	}catch(const char* msg){
		t.produceOutputs("crashed");
		cout << msg << endl;
		exit(1);
	}
	cout << t.getStats() << endl;
	exit(0);
	
}
