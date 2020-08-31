#include "VertexSystem.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <random>
#include <chrono> 

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
		auto start = chrono::high_resolution_clock::now(); 
		t.simulate(generator, unif);
		auto stop = chrono::high_resolution_clock::now(); 
		auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
		cout << ">>Step 1 lasted " << duration.count() << " seconds." << endl;
		if(REPORT_OUT) cout << t.getStats() << endl; 
		if(second_paramfile != ""){
			if(REPORT_OUT) cout << "**** Reading Second parameter file ****"<<endl;
			//t.addSpringsAutomatically();
			//t.produceOutputs("movedspr");
			//t.addAcceptedMovements(moves);
			t.readNewParameters(second_paramfile);
			auto start2 = chrono::high_resolution_clock::now(); 
			t.simulate(generator, unif);
			auto stop2 = chrono::high_resolution_clock::now(); 
			auto duration2 = chrono::duration_cast<chrono::seconds>(stop2 - start2);
			cout << ">>Step 2 lasted " << duration2.count() << " seconds." << endl;
			if(REPORT_OUT) cout << t.getStats() << endl; 
		}
		if(third_paramfile != ""){
			if(REPORT_OUT) cout << "**** Adding springs and Restoring shape ****"<<endl;
			t.addSpringsAutomatically();
			if(REPORT_OUT) cout << "Added springs"<<endl;
			t.makeVeinsThinner(1);
			if(REPORT_OUT) cout << "Veins made thinner"<< endl;
			//t.produceOutputs("movedspr");
			if(REPORT_OUT) cout << "**** Reading Third parameter file ****"<<endl;
			//t.addAcceptedMovements(moves);
			t.readNewParameters(third_paramfile);
			cout << "New Parameters read"<<endl;
			auto start3 = chrono::high_resolution_clock::now(); 
			t.simulate(generator, unif);
			auto stop3 = chrono::high_resolution_clock::now(); 
			auto duration3 = chrono::duration_cast<chrono::seconds>(stop3 - start3);
			cout << ">>Step 3 lasted " << duration3.count() << " seconds." << endl;
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
