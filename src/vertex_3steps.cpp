#include "VertexSystem.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <random>

using namespace std;

int main(int argc, char *argv[]){

	std::string inputfile = argv[1];
	std::string paramfile = argv[2];
    std::string second_paramfile = argv[3];
    std::string third_paramfile = argv[4];
	std::string simulname;
	if(argc < 6){
		simulname = "";
	}else{
		simulname = argv[5];
	}

	if(REPORT_OUT) cout << "Arguments: input condition: " << inputfile << "\n params: " << paramfile << ", " << second_paramfile << ", " << third_paramfile << "\n sim name: " << simulname << endl;
	if(REPORT_OUT) cout << "Declaring Tissue Object. Reading from files\n\n" << endl;
    Tissue t = Tissue(inputfile, paramfile, 0, 0, simulname);

	srand(RANDOM_SEED);
	std::default_random_engine generator;
	std::uniform_real_distribution<double> unif;

	try{
		t.simulate(generator, unif);
		if(REPORT_OUT) cout << t.getStats() << endl; 
		if(second_paramfile != ""){
			if(REPORT_OUT) cout << "**** Reading Second parameter file ****"<<endl;
			//t.addSpringsAutomatically();
			//t.produceOutputs("movedspr");
			//t.addAcceptedMovements(moves);
			t.readNewParameters(second_paramfile);
			t.simulate(generator, unif);
			if(REPORT_OUT) cout << t.getStats() << endl; 
		}
		if(third_paramfile != ""){
			if(REPORT_OUT) cout << "**** Adding springs and Restoring shape ****"<<endl;
			t.addSpringsAutomatically();
			t.restoreShape();
			t.produceOutputs("movedspr");
			if(REPORT_OUT) cout << "**** Reading Third parameter file ****"<<endl;
			//t.addAcceptedMovements(moves);
			t.readNewParameters(third_paramfile);
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