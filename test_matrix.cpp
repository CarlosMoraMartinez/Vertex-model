#include "GXMatrix.h"
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

	GXMatrix<int> g(10,10, 5);
	vector<int> v(10, -2);
	cout << "check element: " << g(1, 1) << endl;
	g(1, 1) = 3;
	cout  << "check changed element: "<< g(1, 1) << endl;
	g.add_row(v);
	cout  << "check added row: " << g(10, 1) << endl;
	
	GXMatrix<int> g2(g);
	g2(10, 1)++;
	cout << "check initialize with another"  << g(10, 1) << " plus 1 = "<< g2(10, 1) << endl;


	cout << "TEST SUM: \n";
	cout << g.toString() << "\n" << "\n" << g2.toString();
	cout << "g + g2: \n";
	GXMatrix<int> g3(g + g2);
	cout << g.toString() << "\n" << "\n" << g2.toString() << "\n" << "result:\n" << g3.toString();


	cout << "TEST MULT: \n";
	cout << g.toString() << "\n" << "\n" << g2.toString();
	cout << "g * g2: \n";
	GXMatrix<int> g4 = g * g2;
	cout << g.toString() << "\n" << "\n" << g2.toString() << "\n" << "result:\n" << g4.toString();


	cout << "TEST DIV: \n";
	cout << g4.toString() << "\n" << "\n" << g2.toString();
	cout << "g4 / g2: \n";
	GXMatrix<int> g5 = g4 / g2;
	cout << g4.toString() << "\n" << "\n" << g2.toString() << "\n" << "result:\n" << g5.toString();
	exit(0);
}
