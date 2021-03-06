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
	std::string second_paramfile, simulname;

	std::string inputfile = argv[1];
	std::string paramfile = argv[2];

	if(argc < 5){
		moves = DEFAULT_STEPS;
		print_step = DEFAULT_WRITE;		
	}else{
		moves = std::stoi(argv[3]);
		print_step = std::stoi(argv[4]);
	}
	if(argc < 6){
		second_paramfile = "";
	}else{
		second_paramfile = argv[5];
	}
	if(argc < 7){
		simulname = "";
	}else{
		simulname = argv[6];
	}
	if(REPORT_OUT) cout << "Reading from file\n\n"<<endl;
	if(REPORT_OUT) cout << "Arguments: input condition: " << inputfile << ", params: " << paramfile << ", moves: " << moves << ", print step: " << print_step << ", sim name: " << simulname << endl;
	Tissue t = Tissue(inputfile, paramfile, moves, print_step, simulname);

	srand(RANDOM_SEED);
	std::default_random_engine generator;
	std::uniform_real_distribution<double> unif;

	try{
		t.simulate(generator, unif);
		if(REPORT_OUT) cout << t.getStats() << endl; 
		if(second_paramfile != ""){
			if(REPORT_OUT) cout << "**** Adding springs ****"<<endl;
			t.addSpringsAutomatically();
			t.produceOutputs("movedspr");
			t.addAcceptedMovements(moves);
			t.readNewParameters(second_paramfile);
			t.simulate(generator, unif);
			if(REPORT_OUT) cout << t.getStats() << endl; 
		}
	}catch(const char* msg){
		t.produceOutputs("crashed");
		cout << msg << endl;
		exit(1);
	}

	if(REPORT_OUT) cout << t.getStats() << endl;
	exit(0);
	
}
