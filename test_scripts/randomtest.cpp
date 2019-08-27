
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <random>
using namespace std;

int main(int argc, char **argv){

	std::srand(1234); //move this to main loop
	std::default_random_engine generator;
	std::uniform_real_distribution<double> unif(0.0,1.0);

	double angle, radius, x, y, move;	
	for(int i = 0; i< atoi(argv[1]); i++){

		angle = 2*M_PI*unif(generator);
		radius = unif(generator) * 0.02;

		x = cos(angle)*radius;
		y =  sin(angle)*radius;
	
		move = unif(generator);

		cout << x << "\t" << y << "\t" << move << std::endl;
	}
}
