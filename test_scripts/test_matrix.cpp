#include "GXMatrix.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <random>

using namespace std;


int main(int argc, char *argv[]){

	GXMatrix<int> g(10,10, 5);
	vector<int> v(10, -2);
	cout << "check element: " << g(1, 1) << endl;
	g(1, 1) = 3;
	cout  << "check changed element: "<< g(1, 1) << endl;
        cout << g.toString() << "<BEFORE ADDING ROW" << endl;
	g.add_row(g.getRow(1));
        cout << g.toString() << "<AFTER ADDING ROW" << endl;
        g(10, 1)++;
	cout  << "check added row element (should be 4): " << g(10, 1)<< ", element in original row (should be 3): " <<  g(1, 1) << endl;
	

	g.add_row(2, 19);
        cout << g.toString() << "<AFTER ADDING 2 ROWS WITH VALUE = 19" << endl;


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


	cout << "ASSIGNMENT: \n";
	cout << g4.toString() << "\n" << "\n" << g2.toString();
	g4 = g2;
	g4(0,0) = 12345;
	cout << "after: \n";
	cout << g4.toString() << "\n" << "\n" << g2.toString();


	cout << "g5 /= g*g*g: \n";
	cout << g5.toString() << "\n" << "\n" << g.toString();
	g5 /= g*g*g;
	g5(0,0) = 12345;
	cout << "after: \n";
	cout << g5.toString() << "\n" << "\n" << g.toString();


	cout << "SCALAR SUM: \n";
	cout << g.toString();
	g2 = g + 100;
	g2(0,0) = 12345;
	cout << "after: \n";
	cout << g.toString() << "\n" << "\n" << g2.toString();

	cout << "TRANSPOSE: \n";
	GXMatrix<double> g6(5, 7, 4.2);
	GXMatrix<double> g7;
	for(int i = 0; i < 5; i++) for(int j = 0; j < 7; j++) g6(i,j)+= double(i*j + j);
	cout << g6.toString();
	g7 = g6.transpose();
	cout << "after: \n";
	cout << g6.toString() << "\n" << "\n" << g7.toString();
	exit(0);
}
