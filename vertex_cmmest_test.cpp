#include "VertexSystem.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <random>


using namespace std;

int main(){
	//Test that enum types work
	CellType cell = CellType::hinge;
	EdgeType edge = EdgeType::tissue_boundary;
	cout << int(cell) << "\t" << int(edge) << endl;	
	
		cout << "\nTEST EDGE VECTOR"<<endl;
	//Test edge_v
	Edge e2;
	edge_v edges;
	for(int i = 0; i < 10; i++){
		e2.ind=i;
		e2.length = i/100;
		e2.tension = i;
		e2.type = EdgeType::hinge;
		e2.vertices[0] = i;
		e2.vertices[1] = i*2;
		e2.cells[0] = i;
		e2.cells[1] = i*2;
		edges.push_back(e2);
	}
	for(Edge e:edges){
		cout << e << endl;
	}
	
	// test cell_v
	cout << "\nTEST CELL VECTOR"<<endl;
	Cell c;
	cell_v cells;
	for(int i = 0; i < 10; i++){
		c.ind=i;
		c.area = i/100;
		c.preferred_area = PREFERRED_AREA_INITIAL;
		c.type = CellType::hinge;
		c.num_vertices = 5;
		c.vertices[0] = i;
		c.vertices[1] = i*2;
		c.vertices[2] = i*3;
		c.vertices[3] = i*4;
		c.vertices[4] = i*5;
		c.edges[0] = i;
		c.edges[1] = i*2;
		c.edges[2] = i*3;
		c.edges[3] = i*4;
		c.edges[4] = i*5;
		cells.push_back(c);
	}
	for(Cell c:cells){
		cout << c << endl;
	}
	
	// Test vertex_v
	cout << "\nTEST VERTEX VECTOR"<<endl;
	vertex_v vertices;
	Vertex v;
	for(int i=0; i< 10; i++){
		v.ind = i;
		v.x = 5.3;
		v.y = 6.2*2;
		v.cells[0] = i;
		v.cells[1] = i + 1;
		v.cells[2] = i + 2;
		v.edges[0] = i;
		v.edges[1] = i*2;
		v.edges[2] = i*3;
		v.neighbour_vertices[0] = i;
		v.neighbour_vertices[1] = i*2;
		v.neighbour_vertices[2] = i*3;
		vertices.push_back(v);
	}
	for(Vertex v:vertices){
		cout << v << endl;
	}
	cout << "Remove element 5"<<endl;
	vertices.erase(vertices.begin() + 5);
	for(Vertex v:vertices){
		cout << v << endl;
	}
	//Testing Tissue initialization
	cout << "Reading from file\n\n"<<endl;
	std::string inputfile = "hex2_s2.5_r20_c20_n0.2_0";//"hex_s2.5_r5_c20_n0.4_0";//"test1016_11";
	Tissue t = Tissue(inputfile);
	cout << "Read\n";
	cout << t;
	cout << "Perimeter in cell 2 should be: " << 0.936 + 0.84 + 0.217 + 1.256 + 0.27 + 1.459 + 0.1911 << endl;
	cout << "Energy of Vertex 0 should be: 2.5948" << endl;
	cout << "Energy of Vertex 26 should be: 1.285702" << endl;
	//Test cell area
	Cell cellx;
	Vertex a, b, c2, d, e;
	a.x = 3; a.y = 4;
	b.x = 5; b.y = 11;
	c2.x = 12; c2.y = 8;
	d.x = 9; d.y = 5;
	e.x = 5; e.y = 6;
	a.ind = 0; b.ind = 1; c2.ind = 2; d.ind = 4; e.ind = 5;
	cellx.num_vertices = 5;
	cellx.ind = 0;
	for(int i = 0; i < cellx.num_vertices; ++i) cellx.vertices[i] = i;
	Tissue t2 = Tissue();
	t2.addVertex(a);
	t2.addVertex(b);
	t2.addVertex(c2);
	t2.addVertex(d);
	t2.addVertex(e);
	t2.addCell(cellx);
	double area = t2.calculateCellArea(cellx);
	cout << t2 <<endl;
	cout << "area test cell should be 30. Area calculated: " << area << endl;

	t.simulate(1e8, "test_sim_", 1e5, true);
	/*
	cout << "Tests move vertices " << endl;

	
	srand(RANDOM_SEED);
	std::default_random_engine generator;
	std::uniform_real_distribution<double> unif;

	std::string fname;
	int fprint = 0;
	for(int i = 0; i< 100000; i++){
		if(i % 1000 == -1){
		fname = inputfile + "_moved_" + std::to_string(fprint);
		t.writeCellsFile(fname);
		t.writePointsFile(fname);
		fprint++;
		}
		t.tryMoveVertex(generator, unif);

	}*/
	cout << t.getStats() << endl;
	exit(0);
	
}
