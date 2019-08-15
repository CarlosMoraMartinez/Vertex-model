

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <random>
#include <limits>
#include "VertexSystem.h"


using namespace std;

//Empty constructor useful to debug
Tissue::Tissue() : num_cells(0), num_vertices(0), num_edges(0), counter_move_trials(0), counter_moves_accepted(0), counter_t1(0), counter_t1_abortions(0), counter_divisions(0), counter_t2(0), counter_t1_outwards(0), counter_t1_inwards(0){
	//Do nothing
}

//Constructor that reads vertex positions and cells from two different files, and initializes all variables using constants defined in VertexSystem.h
Tissue::Tissue(std::string starting_tissue_file, int max_accepted_movements,  int write_every_N_moves) : num_cells(0), num_vertices(0), num_edges(0),  counter_move_trials(0), counter_moves_accepted(0), counter_t1(0), counter_t1_abortions(0), counter_divisions(0), counter_t2(0), counter_t1_outwards(0), counter_t1_inwards(0), simname(starting_tissue_file), max_accepted_movements(max_accepted_movements), write_every_N_moves(write_every_N_moves){
	//Read file of vertices (indicates coordinates for each vertex)
	string vertexfile = starting_tissue_file + VERTEX_FILE_EXTENSION;
	std::ifstream fin_vertex;
	fin_vertex.open(vertexfile);
	initialize_vertices(fin_vertex);
	fin_vertex.close();
	//Read file of cells (indicates vertices for each cell)
	cout << ".points file read...\n";
	try{
		string cellfile = starting_tissue_file + CELLS_FILE_EXTENSION;
		ifstream fin_cells(cellfile);
		initialize_cells(fin_cells);
		fin_cells.close();
		cout << ".cells file read...\n";
	}catch(const char* msg){
		cout << msg << endl;
		exit(1);
	}
	//Initialize edges from vertices and cells
	initialize_edges();	
	//No file with parameters, therefore use constants defined in VertexSystem.h.
	//Also initializes cell area and vertex energy
	set_default_params();
}



/*
Initializes vertices from file
Input: ifstream pointing to a file defining vertex coordinates
*/
void Tissue::initialize_vertices(std::ifstream& inp){
	string s;
	int split_ind;
	getline(inp, s);
	this->num_vertices = stoi(s);
	Vertex v;
	v.dead=false;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		v.cells[i] = EMPTY_CONNECTION;
		v.edges[i] = EMPTY_CONNECTION;
		v.neighbour_vertices[i] = EMPTY_CONNECTION;
	}
	int i = 0;
	while(getline(inp, s)){
		v.ind = i;
		split_ind = s.find("\t");
		v.x = stod(s.substr(0, split_ind));
		v.y = stod(s.substr(split_ind));	
		this->vertices.push_back(v);
		i++;
	}	
}

/*
Initializes cells from file; adds cells to vertices
Input: ifstream pointing to a file defining vertex identifiers (starting from 0) for each cell
*/
void Tissue::initialize_cells(std::ifstream& inp){
	string s, s2;
	int split_ind, split_ind2;
	getline(inp, s);
	split_ind = s.find("\t");
	//cout << "reading first line...\n" << s << split_ind << endl;
	this->num_cells = stoi(s.substr(0, split_ind));
	//cout << "first line read \n" << split_ind << "\t" << this->num_cells <<endl;
	if(stoi(s.substr(split_ind)) > MAX_SIDES_PER_CELL){
		throw "Cell file has more vertices by cels than allowed by MAX_SIDES_PER_CELL constant";
	}
	int vertex, vert_count;

	Cell c;
	c.dead=false;
	int cell_count = 0;
	while(getline(inp, s)){ //Each row in file represents a cell
		vert_count = 0; //Store number of vertices read so far for this cell
		c.ind = cell_count;
		c.num_vertices = 0;
		
		for(int i=0; i < MAX_SIDES_PER_CELL; i++){ //Reset all connections of cell c to EMPTY_CONNECTIONS
			c.vertices[i] = EMPTY_CONNECTION;
			c.edges[i] = EMPTY_CONNECTION; //edges are initialized here but will be added later, in method initialize_edges
		}

		split_ind = 0;	//Position in line where last number was found
		split_ind2 = s.find("\t", split_ind);
		vertex = stod(s.substr(split_ind, split_ind2)); //Number is between split_ind and split_ind2
		while(vertex != EMPTY_CONNECTION){		//Each number in a row is the index of a vertex. Consecutive vertices in a row are adjacent to each other.
			c.vertices[vert_count] = vertex;
			addCellToVertex(vertex, cell_count);
			vert_count++;
			split_ind = split_ind2 + 1;
			split_ind2 = s.find("\t", split_ind);
			vertex = stod(s.substr(split_ind, split_ind2)); //Number is between split_ind and split_ind2		
		}
		c.num_vertices = vert_count;
		this->cells.push_back(c);
		cell_count++;
	}//while rows in file

} 


/*
Initializes edges from vertex and cell data already read from file
Assumes that vertices in cell are in order (i.e., consecutive vertices are bound by an edge, and last vertex is bound to first)
Calls: Tissue::addEdge
*/

void Tissue::initialize_edges(){
	int previous_vertex;
	for(int cell = 0; cell < this->num_cells; cell++){ //for each cell
		previous_vertex = this->cells[cell].vertices[this->cells[cell].num_vertices - 1];
		for(int v = 0; v < this->cells[cell].num_vertices; v++){ //For each vertex in a cell
			addEdge(previous_vertex, this->cells[cell].vertices[v], cell);
			previous_vertex = this->cells[cell].vertices[v];
		}//for vertex in cell
	}//for each cell
}
/*
Initializes a single edge. If edge between two vertices already exists, adds the second cell. 
If the edge does not exist, it creates one and adds the two vertices and the first cell.
Inputs:
index of vertex1
index ofvertex2
index of cell to which v1 and v2 belong
*/
void Tissue::addEdge(int v1, int v2, int cell){
	bool edge_found = false;
	//If edge between v1 and v2 already exists, add Cell to it
	for(int i = 0; i < this->num_edges; i++){
		if((this->edges[i].vertices[0] == v1 && this->edges[i].vertices[1] == v2) ||
			(this->edges[i].vertices[1] == v1 && this->edges[i].vertices[0] == v2)){
				this->edges[i].cells[1] = cell;
				addEdgeToCell(cell, i);
				edge_found = true;
				break;
		}		
	}
	//If edge does not exist, create edge and add the two vertices and the cell passed as argument to it. The other cell will be added later
	if(!edge_found){
		Edge e;
		e.dead=false;
		e.ind = this->num_edges;
		e.vertices[0] = v1;
		e.vertices[1] = v2;
		e.cells[0] = cell;
		e.cells[1] = EMPTY_CONNECTION;
		e.length = sqrt(pow(this->vertices[v1].x - this->vertices[v2].x, 2) + pow(this->vertices[v1].y - this->vertices[v2].y, 2));
		this->edges.push_back(e);
		addEdgeToVertex(v1, e.ind);
		addEdgeToVertex(v2, e.ind);
		addNeighbourVertex(v1, v2);
			
		addEdgeToCell(cell, e.ind);
		this->num_edges++;	
	}
}


void Tissue::set_default_params(){
	for(int c = 0; c < cells.size(); c++){
		cells[c].type = CellType::blade;
		cells[c].preferred_area = PREFERRED_AREA_INITIAL;
		cells[c].area = calculateCellArea(cells[c]);
		cells[c].perimeter_contractility = PERIMETER_CONTRACT_BLADE;
		cells[c].perimeter = calculateCellPerimeter(cells[c]);
	}
	for(int e = 0; e < edges.size(); e++){
		if(contains(EMPTY_CONNECTION, edges[e].cells, VERTEX_PER_EDGE)){
			edges[e].type = EdgeType::tissue_boundary;
			edges[e].tension = LINE_TENSION_TISSUE_BOUNDARY;
		}else{
			edges[e].type = EdgeType::blade;
			edges[e].tension = LINE_TENSION_BLADE;
		}
	}
	for(int v = 0; v < vertices.size(); v++){
		vertices[v].energy = calculateEnergy(vertices[v]);
	}

}

/*
Adds a cell index to the vector of cells in a vertex structure
Input: 
- vertex: the index of vertex to which cell must be added
- cell: the index of the cell to add
*/
void Tissue::addCellToVertex(int vertex, int cell){
		int i = 0;		
		while(this->vertices[vertex].cells[i] != EMPTY_CONNECTION) i++;
		this->vertices[vertex].cells[i] = cell;
}


/*
Adds an edge index to the vector of edges in a vertex structure
Input: 
- vertex: the index of vertex to which cell must be added
- edge: the index of the edge to add
*/
void Tissue::addEdgeToVertex(int vertex, int edge){
		int i = 0;
		while(this->vertices[vertex].edges[i] != EMPTY_CONNECTION) i++;
		this->vertices[vertex].edges[i] = edge;
}

/*
Adds an edge index to the vector of edges in a cell structure
Input: 
- cell: the index of edge to which cell must be added
- edge: the index of the edge to add
*/
void Tissue::addEdgeToCell(int cell, int edge){
	int i = 0;
	while(this->cells[cell].edges[i] != EMPTY_CONNECTION) i++;
	this->cells[cell].edges[i] = edge;
}


void Tissue::addNeighbourVertex(int vertex1, int vertex2){
	int i = 0;
	while(this->vertices[vertex1].neighbour_vertices[i] != EMPTY_CONNECTION) i++;
	this->vertices[vertex1].neighbour_vertices[i] = vertex2;
	i = 0;
	while(this->vertices[vertex2].neighbour_vertices[i] != EMPTY_CONNECTION) i++;
	this->vertices[vertex2].neighbour_vertices[i] = vertex1;
}

//dead_vertices, dead_cells, dead_edges

int Tissue::newVertex(){
	if (!dead_vertices.empty()){
		Vertex* v = &this->vertices[dead_vertices.front()];
		dead_vertices.pop();
		for(int i = 0; i < CELLS_PER_VERTEX; i++){
			v->cells[i] = EMPTY_CONNECTION;
			v->edges[i] = EMPTY_CONNECTION;
			v->neighbour_vertices[i] = EMPTY_CONNECTION;
		}
		v->dead = false;
		this->num_vertices++;
		return v->ind;
	}else{
		Vertex v;
		v.ind = static_cast<int>(this->vertices.size());
		for(int i = 0; i < CELLS_PER_VERTEX; i++){
			v.cells[i] = EMPTY_CONNECTION;
			v.edges[i] = EMPTY_CONNECTION;
			v.neighbour_vertices[i] = EMPTY_CONNECTION;
		}
		v.dead = false;
		this->vertices.push_back(v);
		this->num_vertices++;
		return v.ind;
	}
		

}




int Tissue::newVertex(double x, double y){
	int ind = newVertex();
	vertices[ind].x = x;
	vertices[ind].y = y;
	return ind;
}
int Tissue::newCell(){
	if(!dead_cells.empty()){
		Cell * c = &this->cells[dead_cells.front()];
		dead_cells.pop();
		for(int i = 0; i < MAX_SIDES_PER_CELL; i++){
			c->vertices[i] =  EMPTY_CONNECTION;
			c->edges[i] = EMPTY_CONNECTION;
		}
		c->dead = false;
		c->num_vertices = 0;
		this->num_cells++;
		return c->ind;
	}else{
		Cell c;
		c.ind = static_cast<int>(this->cells.size());
		for(int i = 0; i < MAX_SIDES_PER_CELL; i++){
			c.vertices[i] =  EMPTY_CONNECTION;
			c.edges[i] = EMPTY_CONNECTION;
		}
		c.dead = false;
		c.num_vertices = 0;
		this->num_cells++;
		this->cells.push_back(c);
		return c.ind;		
	}

}
int Tissue::newEdge(){
	if(!dead_edges.empty()){
		//cout << "DEAD EDGES NOT EMPTY\n";
		Edge* e = &this->edges[dead_edges.front()];
		dead_edges.pop();	
		e->vertices[0] = EMPTY_CONNECTION;
		e->vertices[1] = EMPTY_CONNECTION;
		e->cells[0] = EMPTY_CONNECTION;
		e->cells[1] = EMPTY_CONNECTION;
		e->dead=false;
		this->num_edges++;
		return e->ind;
	}else{
		//cout << "here something weird happens\n";
		Edge e;	
		e.ind = static_cast<int>(this->edges.size());
		e.vertices[0] = EMPTY_CONNECTION;
		e.vertices[1] = EMPTY_CONNECTION;
		e.cells[0] = EMPTY_CONNECTION;
		e.cells[1] = EMPTY_CONNECTION;
		e.dead=false;
		this->num_edges++;
		edges.push_back(e);
		return e.ind;
	}



}
///////////////////////////////////////////////////////////////////////////////////////////////
//public functions only useful for testing
int Tissue::addVertex(Vertex v){
	vertices.push_back(v);
	num_vertices++;
	return num_vertices -1;
}
int Tissue::addCell(Cell c){
	cells.push_back(c);
	num_cells++;
	return num_cells - 1;
}
int Tissue::addEdge(Edge e){
	edges.push_back(e);
	num_edges ++;
	return num_edges;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//Area of a polygon using Shoelace formula 
double Tissue::calculateCellArea(const Cell& c){
	double area = 0.0;
	int previous = c.num_vertices - 1;
	for(int i = 0; i < c.num_vertices; i++){
		area += (vertices[c.vertices[previous]].x * vertices[c.vertices[i]].y - vertices[c.vertices[previous]].y * vertices[c.vertices[i]].x);
		previous = i;
	}
	return abs(area/2.0);
}

double Tissue::calculateCellPerimeter(const Cell& c){
	double perim = 0.0;
	for(int i = 0; i < c.num_vertices; i++) perim += edges[c.edges[i]].length;
	return perim;
}

double Tissue::distance(int v1, int v2){
	return sqrt(pow(this->vertices[v1].x - this->vertices[v2].x, 2) + pow(this->vertices[v1].y - this->vertices[v2].y, 2));
}

double Tissue::calculateEnergy(Vertex& v){
	double term1 = 0, term2 = 0, term3 = 0;
	
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(v.cells[i] != EMPTY_CONNECTION) term1 += pow(cells[v.cells[i]].area/cells[v.cells[i]].preferred_area - 1, 2);
	}
	//Second term is edge length divided by preferred area of cell.
	//Most of the time both cells will be of the same type, but sometimes one will be hinge and the other will be blade
	//Therefore, I take the mean between the preferred areas of both cells,
	//except for Edges in the border (obviously)
	double pref_area;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(v.edges[i] != EMPTY_CONNECTION){
			pref_area = edges[v.edges[i]].type == EdgeType::tissue_boundary ? 
										edges[v.edges[i]].cells[0]  == EMPTY_CONNECTION ? 
																cells[ edges[v.edges[i]].cells[1] ].preferred_area
																: cells[ edges[v.edges[i]].cells[0] ].preferred_area	
										: (cells[ edges[v.edges[i]].cells[0] ].preferred_area +  cells[ edges[v.edges[i]].cells[1] ].preferred_area)*0.5;
			term2 += edges[v.edges[i]].tension * edges[v.edges[i]].length/sqrt(pref_area);
		}
		
	}
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(v.cells[i] != EMPTY_CONNECTION) term3 += 0.5*cells[v.cells[i]].perimeter_contractility * pow(cells[v.cells[i]].perimeter, 2)/cells[v.cells[i]].preferred_area;
	}
	return 0.5*term1 + term2 + term3;
}

//
void Tissue::moveVertex(Vertex& v, float x, float y){
	//Set new vertex coordinates. Original object is modified because v is passed by reference
	v.x = x; 
	v.y = y;

	//Loop to re-calculate Edge lengths
	//Length of edge in index v.edges[i] is updated
	//Distance function is called with indexes of vertices touching edge (edge.vertices[0] and edge.vertices[1], where edge is edges[ v.edges[i] ])
	for(int i = 0; i < CELLS_PER_VERTEX; i++){ 
		if(v.edges[i] != EMPTY_CONNECTION) this->edges[v.edges[i]].length = distance(this->edges[v.edges[i]].vertices[0], this->edges[v.edges[i]].vertices[1]);
	}
	
	// Loop to re-calculate Cell areas
	// calculateCellArea is called with the Cell structure located in vector this->cells, in position v.cells[i], as only argument
	// same for calculateCellPerimeter
	for(int i = 0; i < CELLS_PER_VERTEX; i++){ // re-calculate cell areas
		if(v.cells[i] != EMPTY_CONNECTION){
			this->cells[v.cells[i]].area = calculateCellArea(this->cells[v.cells[i]]);
			this->cells[v.cells[i]].perimeter = calculateCellPerimeter(this->cells[v.cells[i]]);
		}
	}
	v.energy = calculateEnergy(v);
}
bool Tissue::tryMoveVertex(std::default_random_engine& generator, std::uniform_real_distribution<double>& unif){

	int vertex_to_move;
	do{
		vertex_to_move = std::rand() % static_cast<int>(vertices.size());
	}while(vertices[vertex_to_move].dead); // A dead vertex cannot be selected
	double old_x = vertices[vertex_to_move].x;
	double old_y = vertices[vertex_to_move].y;
	double old_energy = calculateEnergy(vertices[vertex_to_move]); //Calculate again because it is not updated every time a cell area changes etc
	double angle, radius, new_x, new_y;
	double move_prob;
	bool cell_borders_cross;
	int cross_trials = 0;
	do{
		angle = 2*M_PI*unif(generator);
		radius = unif(generator) * MAX_RANGE_VERTEX_MOVEMENT;
		new_x = old_x + cos(angle)*radius;
		new_y = old_y + sin(angle)*radius;
		moveVertex(vertices[vertex_to_move], new_x, new_y);
		move_prob = vertices[vertex_to_move].energy <= old_energy? 
							exp((vertices[vertex_to_move].energy - old_energy)/TEMPERATURE_NEGATIVE_ENERGY) : 
							exp(-(vertices[vertex_to_move].energy - old_energy)/TEMPERATURE_POSITIVE_ENERGY) ;
		cell_borders_cross = check_if_edges_cross(vertex_to_move);
		if(cell_borders_cross) moveVertex(vertices[vertex_to_move], old_x, old_y);
		cross_trials++;
		if(cross_trials > 10) return false;
	}while(cell_borders_cross);
	double move = unif(generator);

	if(move < move_prob){ 
		detectChangesAfterMove(vertex_to_move);
		return true;
	}else{
		moveVertex(vertices[vertex_to_move], old_x, old_y);
		return false;
	}
	
}

void Tissue::produceOutputs(std::string add_to_name="moved"){
	std:string fname = simname + "_" + add_to_name + "_" + std::to_string(int(counter_moves_accepted / write_every_N_moves));
	writeCellsFile(fname);
	writePointsFile(fname);
	writeAllData(fname);
}

void Tissue::simulate(){
	srand(RANDOM_SEED);
	std::default_random_engine generator;
	std::uniform_real_distribution<double> unif;

	produceOutputs();
	//MAIN SIMULATION LOOP
	do{
		if(tryMoveVertex(generator, unif)){ //Tries a random movement; if accepted:
			counter_moves_accepted++;	
			performRearrangements();	// Performs any transition/rearrangement needed
			if(counter_moves_accepted % write_every_N_moves == 0) produceOutputs();
			//else if(counter_moves_accepted >= 4730400) produceOutputs("moved" + std::to_string(counter_moves_accepted));
		}//End if move accepted
		counter_move_trials++;
	}while(counter_moves_accepted < max_accepted_movements);

}

void Tissue::detectChangesAfterMove(int vertex_moved){

	//Check if T1 transitions are needed
	Edge* ee;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(vertices[vertex_moved].edges[i] == EMPTY_CONNECTION) continue;
		ee = &edges[ vertices[vertex_moved].edges[i] ];
		Vertex* v1 = &vertices[ ee->vertices[0] ];
		Vertex* v2 = &vertices[ ee->vertices[1] ];
		if(ee->length <= T1_TRANSITION_CRITICAL_DISTANCE){ //if length < critical distance
			if(ee->cells[0] != EMPTY_CONNECTION && ee->cells[1] != EMPTY_CONNECTION){ //if the edge touches two cells, and there are 4 cells involved, then T1
				//int empty_connections = 0;
				//for(int j = 0; j < CELLS_PER_VERTEX; j++) if(v1->cells[j] == EMPTY_CONNECTION) empty_connections++;
				//for(int j = 0; j < CELLS_PER_VERTEX; j++) if(v2->cells[j] == EMPTY_CONNECTION) empty_connections++;
				int cellnum = 0;
				for(int j = 0; j < CELLS_PER_VERTEX; j++) if(v1->cells[j] != EMPTY_CONNECTION) cellnum++;
				for(int j = 0; j < CELLS_PER_VERTEX; j++) if(v2->cells[j] != EMPTY_CONNECTION && !contains(v2->cells[j], v1->cells, CELLS_PER_VERTEX)) cellnum++;
				if(cellnum == 4 && cells[ee->cells[0]].num_vertices > 3 && cells[ee->cells[1]].num_vertices > 3 ){ //empty_connections == 0 
					rearrangements_needed.push(Rearrangement{ee->ind, RearrangementType::t1});
				//If less than 4 cells involved, t1 transition outwards //empty_connections == 1
				}else if(cellnum == 3 && cells[ee->cells[0]].num_vertices > 3 && cells[ee->cells[1]].num_vertices > 3){  
					rearrangements_needed.push(Rearrangement{vertices[vertex_moved].edges[i], RearrangementType::t1_at_border_outwards});
				}
			}else{  //Otherwise, length of union of edges touching two vertices touching the edge
				int edgenum = 0;
				for(int j = 0; j < CELLS_PER_VERTEX; j++) if(v2->edges[j] != EMPTY_CONNECTION) edgenum++;
				for(int j = 0; j < CELLS_PER_VERTEX; j++) if(v1->edges[j] != EMPTY_CONNECTION && !contains(v1->edges[j], v2->edges, CELLS_PER_VERTEX)) edgenum++;
				if(edgenum == 5)  rearrangements_needed.push(Rearrangement{vertices[vertex_moved].edges[i], RearrangementType::t1_at_border_inwards});
			}

		} 
	}

	//Check if T2 transitions (cell death) or cell divisions are needed
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(vertices[vertex_moved].cells[i] == EMPTY_CONNECTION) continue;
		if(cells[ vertices[vertex_moved].cells[i] ].area <= T2_TRANSITION_CRITICAL_AREA && cells[ vertices[vertex_moved].cells[i] ].num_vertices == 3){
			rearrangements_needed.push(Rearrangement{vertices[vertex_moved].cells[i], RearrangementType::t2});
		} else if (cells[ vertices[vertex_moved].cells[i] ].area >= MAX_CELL_AREA){
			rearrangements_needed.push(Rearrangement{vertices[vertex_moved].cells[i], RearrangementType::divide_cell});
		}
	}

}//detect rearrangements

void Tissue::performRearrangements(){
	
	Rearrangement r;
	while(! rearrangements_needed.empty()){
		r = rearrangements_needed.front();
		rearrangements_needed.pop();
		switch (r.type){
			case RearrangementType::t1:
				cout << "Entering t1 \n";
				if(T1_ACTIVE) make_t1(r);
				break;
			case RearrangementType::t2:
				cout << "Entering t2 \n";
				if(T2_ACTIVE) make_t2(r);
				break;
			case RearrangementType::divide_cell:
				cout << "Entering division \n";
				if(DIVISION_ACTIVE) make_divide_cell(r);
				break;
			case RearrangementType::divide_edge:
				break;
			case RearrangementType::join_limit_edges:
				//if() make_join_limit_edges(r);
				break;
			case RearrangementType::t1_at_border_outwards:
				cout << "Entering t1 outwards\n";
				if(T1_BORDER_OUTWARDS_ACTIVE) make_t1_at_border_outwards(r);
				break;
			case RearrangementType::t1_at_border_inwards:
				cout << "Entering t1 inwards " << counter_move_trials << "\n";
				if(T1_BORDER_INWARDS_ACTIVE) make_t1_at_border_inwards(r);
				break;
			default:
				break;
		}

	}

}
		

void Tissue::make_t1_at_border_inwards(Rearrangement& r){
	int edge = r.element_index;
	if(edges[edge].length > T1_TRANSITION_CRITICAL_DISTANCE) return; //Check that condition is still true (other rearrangements could have taken place since detection)

	Vertex* v1 = &vertices[ edges[edge].vertices[0] ];
	Vertex* v2 = &vertices[ edges[edge].vertices[1] ];
	//cout << "A, v1: " << v1->ind << " v2: " << v2->ind << endl;
	if(!(edges[edge].cells[0] == EMPTY_CONNECTION || edges[edge].cells[1] == EMPTY_CONNECTION )) return;

	double old_x1 = v1->x, old_x2 = v2->x, old_y1 = v1->y, old_y2 = v2->y, old_length = edges[edge].length;

	//Get the 2 cells that are common to both vertices and cell that is exclusive to v1
	int common_cell1 =  -1, common_cell2 = -1, only_v1 = -1, only_v2 = -1;
	t1_getCellRelationships(v1, v2, common_cell1, common_cell2, only_v1, only_v2); 
	if(common_cell1 == EMPTY_CONNECTION){
		common_cell1 = common_cell2;
		common_cell2 = EMPTY_CONNECTION;
	}
	//cout << "B " << common_cell1 << " " << only_v1 << " " << only_v2 << "\n";
	Cell* cc1 = &this->cells[common_cell1]; //These pointers to Cell structures will be used later
	//Cell* cc2 = &this->cells[common_cell2];
	Cell* sp1 = &this->cells[only_v1];
	Cell* sp2 = &this->cells[only_v2];
	//cout << "Cn";
	if(cc1->num_vertices < 4) return; //One of the cells can't lose any other vertex
	//t1_rotate_edge(v1, v2, edge, cc1, cc2, sp1, sp2);
	t1_rotate_edge90degrees(v1, v2, edge);
	//cout << "D\n";
	//Determine which cell goes with which point. Compare new point locations of v1 and v2 with all the neighbours of v1 and v2. 

	double dist1 = t1_inwards_get_dist_sum(v1, v2, cc1, sp1, sp2);
	double dist2 = t1_inwards_get_dist_sum(v2, v1, cc1, sp2, sp1);

	//cout << "E\n";
	if(dist2 > 0 && (dist2 < dist1 || dist1 < 0)){ //common_cell1 goes with v1, and v2 goes with EMPTY_CONNECTION

		v1 = &vertices[ edges[edge].vertices[1] ];
		v2 = &vertices[ edges[edge].vertices[0] ];
		int aux = only_v1;
		only_v1 = only_v2;
		only_v2 = aux;
		sp1 = &this->cells[only_v1];
		sp2 = &this->cells[only_v2];
	}else if(dist2 < 0 && dist1 < 0){//edges always cross
		v1->x = old_x1;
		v2->x = old_x2;
		v1->y = old_y1;
		v2->y = old_y2;
		edges[edge].length = old_length;
		counter_t1_abortions++;
		return;		
	}
	ofstream ff;
	if(REPORT_T1){
		ff.open("T1_REPORT_inwards_" + std::to_string(counter_moves_accepted) +  "_" + std::to_string(v1->ind) + "_" + std::to_string(v2->ind) + ".txt");

		ff << "\nBEFORE\n\n";
		ff << "edge: " << edge << "\n\n";
		ff << VERTEX_HEADER << *v1 << "\n" << *v2 << "\n\n";
		ff << CELL_HEADER << *cc1 << "\tcommon_cell1\n" << *sp1 << "\tsp1\n"<< *sp2 << "\tsp2\n\n";
		ff << getStats();

		ff << "\n\nDist when cell " << cc1->ind << " goes with vertex " << v1->ind << ": " << dist1 << "\n";
		ff << "Dist when cell " << cc1->ind << " goes with vertex " << v2->ind << ": " << dist2 << "\n";
	}
	//cout << "F\n";
	//Now we know that v1 is closest to common_cell1, and v2 should be closer to common_cell2. Time to update everything
	//1) Cell exclussive of v1 is going to touch v2, and cell exclussive of v2 is going to touch v1
	v1->cells[ which(EMPTY_CONNECTION, v1->cells, CELLS_PER_VERTEX) ] = only_v2;
	v2->cells[ which(common_cell1, v2->cells, CELLS_PER_VERTEX) ] = only_v1;
	//cout << "G\n";
	//2)Remove edge from cell that was common to both vertices
	for(int i = which(edge, cc1->edges, MAX_SIDES_PER_CELL); i < cc1->num_vertices; i++) cc1->edges[i] = i == MAX_SIDES_PER_CELL-1 ? EMPTY_CONNECTION : cc1->edges[i+1];

	//3) Remove one of the vertices from each of the cells that were common to both vertices
	for(int i = which(v2->ind, cc1->vertices, cc1->num_vertices); i < cc1->num_vertices; i++) cc1->vertices[i] = i == MAX_SIDES_PER_CELL-1 ? EMPTY_CONNECTION : cc1->vertices[i+1];
	cc1->num_vertices--;
	//cout << "H\n";
	//4) Add edge to cells that were touched only by one of the vertices
	sp1->edges[sp1->num_vertices] = edge;
	sp2->edges[sp2->num_vertices] = edge;

	//5) Update neighbours 
	int remove_from_v1, remove_from_v2;
	t1_inwards_update_neighbours(v1, v2, edge, common_cell1, only_v1, only_v2, remove_from_v1, remove_from_v2);
	//cout << "I\n";
	//6) Update edges
	int e_remove_from_v1 = -1, e_remove_from_v2 = -1;
	t1_update_edges(v1, v2, edge, remove_from_v1, remove_from_v2, e_remove_from_v1, e_remove_from_v2); //Careful! after update edges, number of vertices still not updated
	//cout << "J\n";	 
	//7) Add v1 to cell that was specific to v2, and v2 to cell that was specific to v1. 
	t1_add_vertices_to_cells(v1, v2, sp2, remove_from_v2);
	t1_add_vertices_to_cells(v2, v1, sp1, remove_from_v1);

	//8)Update edge lengths and areas.
	t1_update_sizes(v1, v2, edge);

	//9) Update cells in edge
	edges[edge].cells[0] = sp1->ind;
	edges[edge].cells[1] = sp2->ind;
	//cout << "K\n";
	//10) Change edge tension and type
	if(sp1->type == CellType::vein){
		if(sp2->type == CellType::blade){
			edges[edge].type = EdgeType::vein_blade;
			edges[edge].tension = LINE_TENSION_VEIN_BLADE;
		}else{
			edges[edge].type = EdgeType::vein_hinge;
			edges[edge].tension = LINE_TENSION_VEIN_HINGE;
		}
	}else if(sp2->type == CellType::vein ){
		if(sp1->type == CellType::blade){
			edges[edge].type = EdgeType::vein_blade;
			edges[edge].tension = LINE_TENSION_VEIN_BLADE;
		}else{
			edges[edge].type = EdgeType::vein_hinge;
			edges[edge].tension = LINE_TENSION_VEIN_HINGE;
		}
	}else if(sp1->type == CellType::hinge || sp2->type == CellType::hinge){ //No special tension for hinge-blade separation; hinge tension used
		edges[edge].type = EdgeType::hinge;
		edges[edge].tension = LINE_TENSION_HINGE;		
	}else{
		edges[edge].type = EdgeType::blade;
		edges[edge].tension = LINE_TENSION_BLADE;
	}

	//cout << "L\n";

	this->counter_t1_inwards++;

	cout << "T1 transition inwards: v1=" << v1->ind << ", v2=" << v2->ind << "; mov. accepted: " << this->counter_moves_accepted << "; T1: " << counter_t1 << "; T1 abortions: " << counter_t1_abortions << ", T1 in: " << counter_t1_inwards << endl;

	if(REPORT_T1){
		ff << "\n\nAFTER\n\n";
		ff << VERTEX_HEADER << *v1 << "\n" << *v2 << "\n\n";
		ff << CELL_HEADER << *cc1 << "\tcommon_cell1\n"<< *sp1 << "\tsp1\n"<< *sp2 << "\tsp2\n\n";
		ff << getStats();
		ff.close();
	}

}// end t1 at border inwards

void Tissue::make_t1_at_border_outwards(Rearrangement& r){
	
	int edge = r.element_index;
	if(edges[edge].length > T1_TRANSITION_CRITICAL_DISTANCE) return; //Check that condition is still true (other rearrangements could have taken place since detection)

	Vertex* v1 = &vertices[ edges[edge].vertices[0] ];	//INSIDE CELL
	Vertex* v2 = &vertices[ edges[edge].vertices[1] ];	//IN BORDER OF CELL
	//cout << "A" << endl;
	int num_cells = 0;
	for(int  i = 0; i < CELLS_PER_VERTEX; i++) if(v1->cells[i] != EMPTY_CONNECTION) num_cells++;
	for(int  i = 0; i < CELLS_PER_VERTEX; i++) if(v2->cells[i] != EMPTY_CONNECTION && !contains(v2->cells[i], v1->cells, CELLS_PER_VERTEX)) num_cells++;
	if(num_cells != 3) return;

	//Check which of the two vertices is facing outwards. The one facing outwards will be v2. Assumes that all vertices inside tissue contact with three cells
	if(! contains(EMPTY_CONNECTION, v2->cells, CELLS_PER_VERTEX)){
		if(contains(EMPTY_CONNECTION, v1->cells, CELLS_PER_VERTEX)){
			v1 = &vertices[ edges[edge].vertices[1] ];
			v2 = &vertices[ edges[edge].vertices[0] ];
		}
	}
	//cout << "B" << endl;
	double old_x1 = v1->x, old_x2 = v2->x, old_y1 = v1->y, old_y2 = v2->y, old_length = edges[edge].length;

	//Get the 2 cells that are common to both vertices and cell that is exclusive to v1
	int common_cell1 =  -1, common_cell2 = -1, only_v1 = -1, only_v2 = -1; //, only_v1_indv1 = -1, only_v2_indv2 = -1;
	t1_getCellRelationships(v1, v2, common_cell1, common_cell2, only_v1, only_v2); //now only_v2 will be EMPTY_CONNECTION 
	Cell* cc1 = &this->cells[common_cell1]; //These pointers to Cell structures will be used later
	Cell* cc2 = &this->cells[common_cell2];
	Cell* sp1 = &this->cells[only_v1];
	//Cell* sp2 = &this->cells[only_v2];
	//cout << "C" << endl;
	if(cc1->num_vertices < 4 || cc2->num_vertices < 4) return; //One of the cells can't lose any other vertex
	//t1_rotate_edge(v1, v2, edge, cc1, cc2, sp1, sp2);
	t1_rotate_edge90degrees(v1, v2, edge);
	//cout << "D" << endl;
	//Determine which cell goes with which point. Compare new point locations of v1 and v2 with all the neighbours of v1 and v2. 

	double dist1 = t1_get_dist_sum(v1, v2, cc1, cc2, sp1, sp1);	//Not very efficient: to check for edge collisions will examine twice sp1
	double dist2 = t1_get_dist_sum(v1, v2, cc2, cc1, sp1, sp1);

	//cout << "E" << endl;
	if(dist2 > 0 && (dist2 < dist1 || dist1 < 0)){
		//cout << "E2" << endl;
		int aux = common_cell1;
		common_cell1 = common_cell2;
		common_cell2 = aux;
		cc1 = &this->cells[common_cell1]; 
		cc2 = &this->cells[common_cell2];
	}else if(dist2 < 0 && dist1 < 0){//edges always cross
		//cout << "E3" << endl;
		v1->x = old_x1;
		v2->x = old_x2;
		v1->y = old_y1;
		v2->y = old_y2;
		edges[edge].length = old_length;
		counter_t1_abortions++;
		return;		
	}
	ofstream ff;
	//cout << "F" << endl;	
	if(REPORT_T1){
		ff.open("T1_REPORT_outwards_" + std::to_string(counter_moves_accepted) +  "_" + std::to_string(v1->ind) + "_" + std::to_string(v2->ind) + ".txt");

		ff << "\nBEFORE\n\n";
		ff << "edge: " << edge << "\n\n";
		ff << VERTEX_HEADER << *v1 << "\n" << *v2 << "\n\n";
		ff << CELL_HEADER << *cc1 << "\tcommon_cell1\n"<< *cc2 << "\tcommon_cell2\n"<< *sp1 << "\tsp1\n\n";
		ff << getStats();

		ff << "\n\nDist when cell " << cc1->ind << " goes with vertex " << v1->ind << ": " << dist1 << "\n";
		ff << "Dist when cell " << cc1->ind << " goes with vertex " << v2->ind << ": " << dist2 << "\n";
	}

	//Now we know that v1 is closest to common_cell1, and v2 should be closer to common_cell2. Time to update everything
	//1) Cell exclussive of v1 is going to touch v2, and cell exclussive of v2 is going to touch v1
	v1->cells[ which(common_cell2, v1->cells, CELLS_PER_VERTEX) ] = EMPTY_CONNECTION;
	v2->cells[ which(common_cell1, v2->cells, CELLS_PER_VERTEX) ] = only_v1;
	//cout << "G" << endl;
	//2)Remove edge from cells that were common to both vertices
	for(int i = which(edge, cc1->edges, MAX_SIDES_PER_CELL); i < cc1->num_vertices; i++) cc1->edges[i] = i == MAX_SIDES_PER_CELL-1 ? EMPTY_CONNECTION : cc1->edges[i+1];
	for(int i = which(edge, cc2->edges, MAX_SIDES_PER_CELL); i < cc2->num_vertices; i++) cc2->edges[i] = i == MAX_SIDES_PER_CELL-1 ? EMPTY_CONNECTION : cc2->edges[i+1];
	//cout << "H" << endl;
	//3) Remove one of the vertices from each of the cells that were common to both vertices
	for(int i = which(v2->ind, cc1->vertices, cc1->num_vertices); i < cc1->num_vertices; i++) cc1->vertices[i] = i == MAX_SIDES_PER_CELL-1 ? EMPTY_CONNECTION : cc1->vertices[i+1];
	cc1->num_vertices--;
	for(int i = which(v1->ind, cc2->vertices, cc2->num_vertices); i < cc2->num_vertices; i++) cc2->vertices[i] = i == MAX_SIDES_PER_CELL-1 ? EMPTY_CONNECTION : cc2->vertices[i+1];
	cc2->num_vertices--;
	//cout << "I" << endl;
	//4) Add edge to cells that were touched only by one of the vertices
	sp1->edges[sp1->num_vertices] = edge;
	//sp2->edges[sp2->num_vertices] = edge;
	//cout << "J" << endl;
	//5) Update neighbours 
	int remove_from_v1, remove_from_v2;
	t1_update_neighbours(v1, v2, edge, common_cell1, common_cell2, only_v1, only_v2, remove_from_v1, remove_from_v2);
	//cout << "K" << endl;
	//6) Update edges
	int e_remove_from_v1 = -1, e_remove_from_v2 = -1;
	t1_update_edges(v1, v2, edge, remove_from_v1, remove_from_v2, e_remove_from_v1, e_remove_from_v2); //Careful! after update edges, number of vertices still not updated
	//cout << "L" << endl;	

	//7) Add v1 to cell that was specific to v2, and v2 to cell that was specific to v1. 
	//t1_add_vertices_to_cells(v1, v2, sp2, remove_from_v2);
	t1_add_vertices_to_cells(v2, v1, sp1, remove_from_v1);

	//8)Update edge lengths and areas.
	t1_update_sizes(v1, v2, edge);
	//cout << "M" << endl;
	//9) Update cells in edge
	edges[edge].cells[0] = sp1->ind;
	edges[edge].cells[1] = EMPTY_CONNECTION;
	edges[edge].tension = LINE_TENSION_TISSUE_BOUNDARY;
	edges[edge].type = EdgeType::tissue_boundary;
	this->counter_t1_outwards++;

	cout << "T1 transition outwards: v1=" << v1->ind << ", v2=" << v2->ind << "; mov. accepted: " << this->counter_moves_accepted << "; T1: " << counter_t1 << "; T1 abortions: " << counter_t1_abortions << ", T1 out: " << counter_t1_outwards <<endl;

	if(REPORT_T1){
		ff << "\n\nAFTER\n\n";
		ff << VERTEX_HEADER << *v1 << "\n" << *v2 << "\n\n";
		ff << CELL_HEADER << *cc1 << "\tcommon_cell1\n"<< *cc2 << "\tcommon_cell2\n"<< *sp1 << "\tsp1\n\n";
		ff << getStats();
		ff.close();
	}

}/// End make_t1_outwards

void Tissue::make_divide_cell(Rearrangement& r){
	int cell = r.element_index;
	if(cells[cell].area < MAX_CELL_AREA || cells[cell].num_vertices < 3 || cells[cell].dead) return; 

	int mv1, mv2;
	double dist, max_dist = -1;
	//get longest ditance
	for(int i = 0; i < cells[cell].num_vertices - 1; i++){
		for(int j  = i + 1; j < cells[cell].num_vertices; j++){
			dist = distance(cells[cell].vertices[i], cells[cell].vertices[j]);
			if(dist > max_dist || max_dist < 0){
				max_dist = dist;
				mv1 = cells[cell].vertices[i];
				mv2 = cells[cell].vertices[j];
			} 
		}

	}

	double center_x = (vertices[mv1].x + vertices[mv2].x)*0.5;
	double center_y = (vertices[mv1].y + vertices[mv2].y)*0.5;
	double angle = atan2(vertices[mv1].y - vertices[mv2].y, vertices[mv1].x - vertices[mv2].x);
	double angle_rotate = DIVISION_ANGLE_RANDOM_NOISE > 0 ? DIVISION_ANGLE_RANDOM_NOISE*(rand()%360)*M_PI/180 + 0.5*M_PI : 0.5*M_PI ; //ADD SOME NOISE TO THE ANGLE
	//Rotate angle 90ºC and calculate new positions (new edge will have a length of 1.5*T1_TRANSITION_CRITICAL_DISTANCE)
	double x1 = center_x + cos(angle - angle_rotate)*max_dist*20; 
	double y1 = center_y - sin(angle - angle_rotate)*max_dist*20; // calculate an orthogonal line (big in excess to be sure that it cuts the polygon in 2 pieces) and look which edges it cuts
	double x2 = center_x - cos(angle - angle_rotate)*max_dist*20;
	double y2 = center_y + sin(angle - angle_rotate)*max_dist*20;

	//cout << "C\n";
	//Find edges to cut
	int e1 = -1, e2 = -1;
	StraightLine l1, l2, l3, l4;
	l1 = getLineFromEdge(x1, x2, y1, y2);
	for(int i = 0; i < cells[cell].num_vertices; i++){
		l2 = getLineFromEdge(&this->edges[cells[cell].edges[i]]);
		if(lines_cross(l1, l2)){
			if(e1 < 0){
				e1 = cells[cell].edges[i];
				l3 = l2;
			}else{
				e2 = cells[cell].edges[i];
				l4 = l2;
				break;
			}

		}

	}//end for find edges that intersect with division edge

	if(e1 < 0  || e2 < 0){
		cout << "Error: new edge does not cut other cell edges. Cell: " << cell <<endl;
		return;
	}
	x1 = (l3.intercept - l1.intercept)/(l1.slope - l3.slope); //Now calculate the points where they actually cross
	x2 = (l4.intercept - l1.intercept)/(l1.slope - l4.slope);
	y1 = l1.intercept + l1.slope*x1;
	y2 = l1.intercept + l1.slope*x2;

	int newvind1 = newVertex(x1, y1);
	int newvind2 = newVertex(x2, y2);  //create new vertices that are going to be positioned at (x1, y1) and (x2, y2) insideedges e1 and e2
	int newcind = newCell();
	int newe = newEdge();

	splitEdgeWithVertex(e1, cell, newvind1); //Changes old cell but does nothing to new (only splits edges, not cell)
	splitEdgeWithVertex(e2, cell, newvind2);
	

	//connect new edge
	edges[newe].cells[0] = cell;
	edges[newe].cells[1] = newcind;
	edges[newe].vertices[0] = newvind1;
	edges[newe].vertices[1] = newvind2;
	edges[newe].length = distance(newvind1, newvind2);
	edges[newe].tension = edges[e1].type == EdgeType::tissue_boundary? edges[e2].tension : edges[e1].tension;//make this more sophisticated to take into account veins etc.
	edges[newe].type = edges[e1].type == EdgeType::tissue_boundary? edges[e2].type : edges[e1].type; //make this more sophisticated to take into account veins etc.

	cells[cell].edges[which(EMPTY_CONNECTION, cells[cell].edges, MAX_SIDES_PER_CELL)] = newe;
	cells[newcind].edges[which(EMPTY_CONNECTION, cells[newcind].edges, MAX_SIDES_PER_CELL)] = newe;

	vertices[newvind2].cells[which(EMPTY_CONNECTION, vertices[newvind2].cells, CELLS_PER_VERTEX)] = newcind;
	vertices[newvind1].cells[which(EMPTY_CONNECTION, vertices[newvind1].cells, CELLS_PER_VERTEX)] = newcind;
	vertices[newvind2].edges[which(EMPTY_CONNECTION, vertices[newvind2].edges, CELLS_PER_VERTEX)] = newe;
	vertices[newvind1].edges[which(EMPTY_CONNECTION, vertices[newvind1].edges, CELLS_PER_VERTEX)] = newe;
	addNeighbourVertex(newvind1, newvind2);

	int newv1_in_c1 = which(newvind1, cells[cell].vertices,  cells[cell].num_vertices);
	int newv2_in_c1 = which(newvind2,  cells[cell].vertices,  cells[cell].num_vertices);

	//Move half of the vertices of one cell to the other, and update edges accordingly
	int i = (newv1_in_c1 + 1) % cells[cell].num_vertices;
	int changeCounter = 0;
	while(i != newv2_in_c1){
		for(int j = 0; j < cells[cell].num_vertices; j++){
			newe =  cells[cell].edges[j];
			if(edges[newe].vertices[0] == cells[cell].vertices[i] || edges[newe].vertices[1] == cells[cell].vertices[i]){
				cells[newcind].edges[which(EMPTY_CONNECTION, cells[newcind].edges, MAX_SIDES_PER_CELL)] = newe;
				cells[cell].edges[j] = EMPTY_CONNECTION;
				edges[newe].cells[which(cell, edges[newe].cells, 2)] = newcind;
			}
		}
		cells[newcind].vertices[cells[newcind].num_vertices] = cells[cell].vertices[i];
		cells[newcind].num_vertices++;
		vertices[cells[cell].vertices[i]].cells[which(cell, vertices[cells[cell].vertices[i]].cells, CELLS_PER_VERTEX)] = newcind;
		cells[cell].vertices[i] = EMPTY_CONNECTION;
		changeCounter ++;
		i = (i+1)%cells[cell].num_vertices;  //Do not decrease cell->num_vertices !
	}
	cells[newcind].vertices[cells[newcind].num_vertices] = newvind2;
	cells[newcind].vertices[cells[newcind].num_vertices + 1] = newvind1;
	cells[newcind].num_vertices+=2;

	//Reagrupate edges and vertices in cell
	for(int j = 0; j < MAX_SIDES_PER_CELL - 1; j++){
		if(cells[cell].edges[j] == EMPTY_CONNECTION){
			for(int k = j + 1; k < MAX_SIDES_PER_CELL; k++){
				if(cells[cell].edges[k] != EMPTY_CONNECTION){
					cells[cell].edges[j] = cells[cell].edges[k];
					cells[cell].edges[k] = EMPTY_CONNECTION;
					break;
				}
			}
		}//end if reagrupate vertices
		if(cells[cell].vertices[j] == EMPTY_CONNECTION){
			for(int k = j + 1; k < MAX_SIDES_PER_CELL; k++){
				if(cells[cell].vertices[k] != EMPTY_CONNECTION){
					cells[cell].vertices[j] = cells[cell].vertices[k];
					cells[cell].vertices[k] = EMPTY_CONNECTION;
					break;
				}
			}//end if reagrupate edges
		}
	}	
	cells[cell].num_vertices -=  changeCounter;//newv1_in_c1 < newv2_in_c1 ? newv2_in_c1 - newv1_in_c1 - 1 : newv1_in_c1 - newv2_in_c1 + cells[cell].num_vertices - 1;

	//Now update cell areas and perimeters
	cells[newcind].type =  cells[cell].type;
	cells[newcind].perimeter_contractility = cells[cell].perimeter_contractility;
	cells[newcind].preferred_area = cells[cell].preferred_area;
	cells[newcind].area = calculateCellArea(cells[newcind]);
	cells[newcind].perimeter = calculateCellPerimeter(cells[newcind]);
	cells[cell].area = calculateCellArea(cells[cell]);
	cells[cell].perimeter = calculateCellPerimeter(cells[cell]);

	counter_divisions++;

	if(REPORT_DIV) writeAllData(simname + "_div_2" + to_string(counter_divisions));
	cout << "DIVISION AT " << counter_moves_accepted << " divi. accepted: " << counter_divisions << " Cell: " << cell << ", New cell: " << newcind << " cut in vertices: " << newvind1 << ", " << newvind2 << endl;
	
	
}// END make_divide_cell

void Tissue::splitEdgeWithVertex(int e, int cell, int  v){

	int i = 0;
	int n1 = edges[e].vertices[0];
	int n2 = edges[e].vertices[1];

	//Insert vertex in cell
	while(cells[cell].vertices[i] != n1) i++;
	if(cells[cell].vertices[(i + 1)%cells[cell].num_vertices] == n2){
		i = (i + 1)%cells[cell].num_vertices;
	}
	for(int j = cells[cell].num_vertices; j > i; j--){
		cells[cell].vertices[j] = cells[cell].vertices[(j - 1 + cells[cell].num_vertices)%cells[cell].num_vertices];
	}

	cells[cell].vertices[i] = v;
	cells[cell].num_vertices++;
	vertices[v].cells[0] = edges[e].cells[0];
	vertices[v].cells[1] = edges[e].cells[1];
	vertices[v].neighbour_vertices[0] = n1;
	vertices[v].neighbour_vertices[1] = n2;
	vertices[n1].neighbour_vertices[which(n2, vertices[n1].neighbour_vertices, CELLS_PER_VERTEX)] = v;
	vertices[n2].neighbour_vertices[which(n1, vertices[n2].neighbour_vertices, CELLS_PER_VERTEX)] = v;

	int e2ind = newEdge();
	edges[e].vertices[which(n2, edges[e].vertices, 2)] = v;
	edges[e2ind].vertices[0] = v;
	edges[e2ind].vertices[1] = n2;
	edges[e2ind].cells[0] = edges[e].cells[0];
	edges[e2ind].cells[1] = edges[e].cells[1];
	edges[e2ind].tension = edges[e].tension;
	edges[e2ind].type = edges[e].type;	
	edges[e].length = distance(v, n1);
	edges[e2ind].length = distance(v, n2);

	vertices[v].edges[0] = e;
	vertices[v].edges[1] = e2ind;
	vertices[n2].edges[which(e, vertices[n2].edges, CELLS_PER_VERTEX)] = e2ind; 
	cells[cell].edges[which(EMPTY_CONNECTION, cells[cell].edges, MAX_SIDES_PER_CELL)] = e2ind; //add edge to dividing cell

	int other_cell = edges[e].cells[0] == cell? edges[e].cells[1] : edges[e].cells[0]; //neighbour cell
	if(other_cell != EMPTY_CONNECTION){
		cells[other_cell].edges[cells[other_cell].num_vertices] = e2ind; //add edge to the neighbour cell
		i = which(n1, cells[other_cell].vertices, cells[other_cell].num_vertices); //neighbour cell also needs the new vertex
		if(cells[other_cell].vertices[(i + 1)%cells[other_cell].num_vertices] == n2){
			i = (i + 1)%cells[other_cell].num_vertices;
		}
		for(int j = cells[other_cell].num_vertices; j > i; j--){
			cells[other_cell].vertices[j] = cells[other_cell].vertices[(j - 1 + cells[other_cell].num_vertices)%cells[other_cell].num_vertices];
		}
		cells[other_cell].vertices[i] = v;
		cells[other_cell].num_vertices++;

	}
}//End splitEdge


void Tissue::make_t1(Rearrangement& r){


	int edge = r.element_index;
	if(edges[edge].length > T1_TRANSITION_CRITICAL_DISTANCE) return; //Check that condition is still true (other rearrangements could have taken place since detection)

	Vertex* v1 = &vertices[ edges[edge].vertices[0] ];
	Vertex* v2 = &vertices[ edges[edge].vertices[1] ];

	if(contains(EMPTY_CONNECTION, v1->cells, CELLS_PER_VERTEX) || contains(EMPTY_CONNECTION, v2->cells, CELLS_PER_VERTEX) || contains(EMPTY_CONNECTION, v1->edges, CELLS_PER_VERTEX) || contains(EMPTY_CONNECTION, v2->edges, CELLS_PER_VERTEX)) return;

	double old_x1 = v1->x, old_x2 = v2->x, old_y1 = v1->y, old_y2 = v2->y, old_length = edges[edge].length;

	//Get the 2 cells that are common to both vertices and cell that is exclusive to v1
	int common_cell1 =  -1, common_cell2 = -1, only_v1 = -1, only_v2 = -1; //, only_v1_indv1 = -1, only_v2_indv2 = -1;
	t1_getCellRelationships(v1, v2, common_cell1, common_cell2, only_v1, only_v2);
	Cell* cc1 = &this->cells[common_cell1]; //These pointers to Cell structures will be used later
	Cell* cc2 = &this->cells[common_cell2];
	Cell* sp1 = &this->cells[only_v1];
	Cell* sp2 = &this->cells[only_v2];

	if(cc1->num_vertices < 4 || cc2->num_vertices < 4) return; //One of the cells can't lose any other vertex
	//t1_rotate_edge(v1, v2, edge, cc1, cc2, sp1, sp2);
	t1_rotate_edge90degrees(v1, v2, edge);

	//Determine which cell goes with which point. Compare new point locations of v1 and v2 with all the neighbours of v1 and v2. 

	double dist1 = t1_get_dist_sum(v1, v2, cc1, cc2, sp1, sp2);
	double dist2 = t1_get_dist_sum(v1, v2, cc2, cc1, sp1, sp2);


	if(dist2 > 0 && (dist2 < dist1 || dist1 < 0)){
		int aux = common_cell1;
		common_cell1 = common_cell2;
		common_cell2 = aux;
		cc1 = &this->cells[common_cell1]; 
		cc2 = &this->cells[common_cell2];
	}else if(dist2 < 0 && dist1 < 0){//edges always cross
		v1->x = old_x1;
		v2->x = old_x2;
		v1->y = old_y1;
		v2->y = old_y2;
		edges[edge].length = old_length;
		counter_t1_abortions++;
		return;		
	}
	ofstream ff;
	
	if(REPORT_T1){
		ff.open("T1_REPORT_" + std::to_string(counter_moves_accepted) +  "_" + std::to_string(v1->ind) + "_" + std::to_string(v2->ind) + ".txt");

		ff << "\nBEFORE\n\n";
		ff << "edge: " << edge << "\n\n";
		ff << VERTEX_HEADER << *v1 << "\n" << *v2 << "\n\n";
		ff << CELL_HEADER << *cc1 << "\tcommon_cell1\n"<< *cc2 << "\tcommon_cell2\n"<< *sp1 << "\tsp1\n"<< *sp2 << "\tsp2\n\n";
		ff << getStats();

		ff << "\n\nDist when cell " << cc1->ind << " goes with vertex " << v1->ind << ": " << dist1 << "\n";
		ff << "Dist when cell " << cc1->ind << " goes with vertex " << v2->ind << ": " << dist2 << "\n";
	}

	//Now we know that v1 is closest to common_cell1, and v2 should be closer to common_cell2. Time to update everything
	//1) Cell exclussive of v1 is going to touch v2, and cell exclussive of v2 is going to touch v1
	v1->cells[ which(common_cell2, v1->cells, CELLS_PER_VERTEX) ] = only_v2;
	v2->cells[ which(common_cell1, v2->cells, CELLS_PER_VERTEX) ] = only_v1;

	//2)Remove edge from cells that were common to both vertices
	for(int i = which(edge, cc1->edges, MAX_SIDES_PER_CELL); i < cc1->num_vertices; i++) cc1->edges[i] = i == MAX_SIDES_PER_CELL-1 ? EMPTY_CONNECTION : cc1->edges[i+1];
	for(int i = which(edge, cc2->edges, MAX_SIDES_PER_CELL); i < cc2->num_vertices; i++) cc2->edges[i] = i == MAX_SIDES_PER_CELL-1 ? EMPTY_CONNECTION : cc2->edges[i+1];

	//3) Remove one of the vertices from each of the cells that were common to both vertices
	for(int i = which(v2->ind, cc1->vertices, cc1->num_vertices); i < cc1->num_vertices; i++) cc1->vertices[i] = i == MAX_SIDES_PER_CELL-1 ? EMPTY_CONNECTION : cc1->vertices[i+1];
	cc1->num_vertices--;
	for(int i = which(v1->ind, cc2->vertices, cc2->num_vertices); i < cc2->num_vertices; i++) cc2->vertices[i] = i == MAX_SIDES_PER_CELL-1 ? EMPTY_CONNECTION : cc2->vertices[i+1];
	cc2->num_vertices--;

	//4) Add edge to cells that were touched only by one of the vertices
	sp1->edges[sp1->num_vertices] = edge;
	sp2->edges[sp2->num_vertices] = edge;

	//5) Update neighbours 
	int remove_from_v1, remove_from_v2;
	t1_update_neighbours(v1, v2, edge, common_cell1, common_cell2, only_v1, only_v2, remove_from_v1, remove_from_v2);

	//6) Update edges
	int e_remove_from_v1 = -1, e_remove_from_v2 = -1;
	t1_update_edges(v1, v2, edge, remove_from_v1, remove_from_v2, e_remove_from_v1, e_remove_from_v2); //Careful! after update edges, number of vertices still not updated
	

	//7) Add v1 to cell that was specific to v2, and v2 to cell that was specific to v1. 
	t1_add_vertices_to_cells(v1, v2, sp2, remove_from_v2);
	t1_add_vertices_to_cells(v2, v1, sp1, remove_from_v1);

	//8)Update edge lengths and areas.
	t1_update_sizes(v1, v2, edge);

	//9) Update cells in edge
	edges[edge].cells[0] = sp1->ind;
	edges[edge].cells[1] = sp2->ind;
	this->counter_t1++;

	cout << "T1 transition: v1=" << v1->ind << ", v2=" << v2->ind << "; mov. accepted: " << this->counter_moves_accepted << "; T1: " << counter_t1 << "; T1 abortions: " << counter_t1_abortions <<endl;
	if(REPORT_T1){
		ff << "\n\nAFTER\n\n";
		ff << VERTEX_HEADER << *v1 << "\n" << *v2 << "\n\n";
		ff << CELL_HEADER << *cc1 << "\tcommon_cell1\n"<< *cc2 << "\tcommon_cell2\n"<< *sp1 << "\tsp1\n"<< *sp2 << "\tsp2\n\n";
		ff << getStats();
		ff.close();
	}
}//End make_t1


void Tissue::t1_update_sizes(Vertex* v1, Vertex* v2, int edge){

	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(v1->edges[i] != edge) this->edges[v1->edges[i]].length = distance(this->edges[v1->edges[i]].vertices[0], this->edges[v1->edges[i]].vertices[1]);
		if(v2->edges[i] != edge) this->edges[v2->edges[i]].length = distance(this->edges[v2->edges[i]].vertices[0], this->edges[v2->edges[i]].vertices[1]);	
	}

	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(v1->cells[i] != EMPTY_CONNECTION){
			this->cells[v1->cells[i]].area = calculateCellArea(this->cells[v1->cells[i]]);
			this->cells[v1->cells[i]].perimeter = calculateCellPerimeter(this->cells[v1->cells[i]]);
		}
	if(v2->cells[i] != EMPTY_CONNECTION){
			this->cells[v2->cells[i]].area = calculateCellArea(this->cells[v2->cells[i]]);
			this->cells[v2->cells[i]].perimeter =calculateCellPerimeter(this->cells[v2->cells[i]]);
		}
	}
}// End update_sizes


//7) Add v1 to cell that was specific to v2, and v2 to cell that was specific to v1. 
//Insertion must keep order of vertices in each cell
//Search rm_from_v2 vertex in sp2 cell (this vertex now is neighbour v1)
//Since vertex rm_from_v2 is touching common_cell1 (which now is specific to vertex v1), v1 should be added between this vertex and v2.
void Tissue::t1_add_vertices_to_cells(Vertex* v1, Vertex* v2, Cell* sp2, int remove_from_v2){

	int move_from;
	for(int i = 0; i < sp2->num_vertices; i++){
		if(sp2->vertices[i] == remove_from_v2){
			if(sp2->vertices[(i + 1) % sp2->num_vertices] == v2->ind){
				move_from = (i+1) % sp2->num_vertices;
			}else if(sp2->vertices[(i - 1 + sp2->num_vertices) % sp2->num_vertices] == v2->ind){
				move_from = i;
			}else{
				cout << "Error in t1 transition: unable to insert vertex in cell: v1=" << v1->ind << " v2=" << v2->ind << " cell:" << sp2->ind << endl;
			}
			break;
		}

	}

	for(int i = sp2->num_vertices; i > move_from; i--){
		sp2->vertices[i] = sp2->vertices[i - 1];
	}

	sp2->vertices[move_from] = v1->ind;
	sp2->num_vertices++;

}// end t1_add_vertices_to_cells

void Tissue::t1_update_edges(Vertex* v1, Vertex* v2, int edge, int remove_from_v1, int remove_from_v2, int& e_remove_from_v1, int& e_remove_from_v2){
	Edge* ee;
	int e_ind_remove_from_v1 = -1, e_ind_remove_from_v2;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		ee = &this->edges[v1->edges[i]];
		if(v1->edges[i] != edge && contains(remove_from_v1, ee->vertices, 2)){
			ee->vertices[which(v1->ind, ee->vertices, 2)] = v2->ind;
			e_ind_remove_from_v1 = i;
			e_remove_from_v1 = ee->ind;
			break;
		}	

	}
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		ee = &this->edges[v2->edges[i]];
		if(v2->edges[i] != edge && contains(remove_from_v2, ee->vertices, 2)){
			ee->vertices[which(v2->ind, ee->vertices, 2)] = v1->ind;
			e_ind_remove_from_v2 = i;
			e_remove_from_v2 = ee->ind;
			break;
		}	
	}
	v1->edges[e_ind_remove_from_v1] = e_remove_from_v2;
	v2->edges[e_ind_remove_from_v2] = e_remove_from_v1;
}// end update edges


void Tissue::t1_inwards_update_neighbours(Vertex* v1, Vertex* v2, int edge, int common_cell1, int only_v1, int only_v2, int& remove_from_v1, int& remove_from_v2){
	//cout << "h2" << endl;
	int ind_remove_from_v1, ind_remove_from_v2;
	Vertex* neighbour;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		neighbour = &this->vertices[v1->neighbour_vertices[i]];
		if(contains(EMPTY_CONNECTION, neighbour->cells, CELLS_PER_VERTEX) && contains(only_v1, neighbour->cells, CELLS_PER_VERTEX) && !contains(common_cell1, neighbour->cells, CELLS_PER_VERTEX) && neighbour->ind != v2->ind){ //LAST CONDITION IS NEEDED FOR ELEMENTS IN CORNERS, SINCE THE WRONG NEIGHBOUR CAN ALSO BE TOUCHING EMPTY_CONNECTION
			ind_remove_from_v1 = i;
			remove_from_v1 = neighbour->ind;
			neighbour->neighbour_vertices[which(v1->ind, neighbour->neighbour_vertices, CELLS_PER_VERTEX)] = v2->ind; //change neighbour's neighbours also
			break;
		}

	}
	//cout << "h2" << endl;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		neighbour = &this->vertices[v2->neighbour_vertices[i]];
		//cout << "nei: " << neighbour->ind << " cc1 " << common_cell1 << " sp2 " << only_v2 << endl;
		if(contains(common_cell1, neighbour->cells, CELLS_PER_VERTEX) && contains(only_v2, neighbour->cells, CELLS_PER_VERTEX) && neighbour->ind != v1->ind){
			ind_remove_from_v2 = i;
			remove_from_v2 = neighbour->ind;
			neighbour->neighbour_vertices[which(v2->ind ,neighbour->neighbour_vertices, CELLS_PER_VERTEX)] = v1->ind; //change neighbour's neighbours also
			break;
		}
	}
	//cout << "h3 " << remove_from_v1 << " " << ind_remove_from_v1 << " " << remove_from_v2 << " " << ind_remove_from_v2 << endl;
	v1->neighbour_vertices[ind_remove_from_v1] = remove_from_v2;
	v2->neighbour_vertices[ind_remove_from_v2] = remove_from_v1;
	//cout << "h4" << endl;

}// End update neighbours t1 inwards


void Tissue::t1_update_neighbours(Vertex* v1, Vertex* v2, int edge, int common_cell1, int common_cell2, int only_v1, int only_v2, int& remove_from_v1, int& remove_from_v2){

	//5) Update neighbour vertices
	//V1 has to remove a neighbour touching common_cell2 and only_v1 cells
	//V2 has to remove a neighbour touching common_cell1 and only_v2 cells
	//V1 has to add a neighbour of v2 that is touching common_cell1 and only_v2 cells (the one removed from v2)
	//V2 has to add a neighbour of v1 that is touching common_cell2 and only_v1 cells (the one removed from v1)
	//cout << "h2" << endl;
	int ind_remove_from_v1, ind_remove_from_v2;
	Vertex* neighbour;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		neighbour = &this->vertices[v1->neighbour_vertices[i]];
		if(contains(common_cell2, neighbour->cells, CELLS_PER_VERTEX) && contains(only_v1, neighbour->cells, CELLS_PER_VERTEX) && neighbour->ind != v2->ind){
			ind_remove_from_v1 = i;
			remove_from_v1 = neighbour->ind;
			neighbour->neighbour_vertices[which(v1->ind, neighbour->neighbour_vertices, CELLS_PER_VERTEX)] = v2->ind; //change neighbour's neighbours also
			break;
		}

	}
	//cout << "h2" << endl;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		neighbour = &this->vertices[v2->neighbour_vertices[i]];
		//cout << "nei: " << neighbour->ind << " cc1 " << common_cell1 << " cc2 " << common_cell2 << " sp2 " << only_v2 << endl;
		if(contains(common_cell1, neighbour->cells, CELLS_PER_VERTEX) && contains(only_v2, neighbour->cells, CELLS_PER_VERTEX) && neighbour->ind != v1->ind){
			ind_remove_from_v2 = i;
			remove_from_v2 = neighbour->ind;
			neighbour->neighbour_vertices[which(v2->ind ,neighbour->neighbour_vertices, CELLS_PER_VERTEX)] = v1->ind; //change neighbour's neighbours also
			break;
		}
	}
	//cout << "h3 " << remove_from_v1 << " " << ind_remove_from_v1 << " " << remove_from_v2 << " " << ind_remove_from_v2 << endl;
	v1->neighbour_vertices[ind_remove_from_v1] = remove_from_v2;
	v2->neighbour_vertices[ind_remove_from_v2] = remove_from_v1;
	//cout << "h4" << endl;

}// End update neighbours


void Tissue::t1_rotate_edge90degrees(Vertex* v1, Vertex* v2, int edge){
	double center_x = (v1->x + v2->x)*0.5;
	double center_y = (v1->y + v2->y)*0.5;
	double angle = atan2(v1->y - v2->y, v1->x - v2->x);
	//Rotate angle 90ºC and calculate new positions (new edge will have a length of 1.5*T1_TRANSITION_CRITICAL_DISTANCE)
	v1->x = center_x + cos(angle - 0.5*M_PI)*LENGTH_ROTATED_EDGE;
	v2->x = center_x - cos(angle - 0.5*M_PI)*LENGTH_ROTATED_EDGE;

	v1->y= center_y + sin(angle - 0.5*M_PI)*LENGTH_ROTATED_EDGE;
	v2->y = center_y - sin(angle - 0.5*M_PI)*LENGTH_ROTATED_EDGE; 
	edges[edge].length = distance(v1->ind, v2->ind);

}//end rotate_edge


void Tissue::t1_rotate_edge(Vertex* v1, Vertex* v2, int edge, Cell* c1, Cell* c2, Cell*sp1, Cell* sp2){


	int v1_in_sp1 = which(v1->ind, sp1->vertices, sp1->num_vertices);
	int v2_in_sp2 = which(v2->ind, sp2->vertices, sp2->num_vertices);
	double y1 = abs( vertices[sp1->vertices[(v1_in_sp1 + 1)%sp1->num_vertices]].y - vertices[sp1->vertices[(v1_in_sp1 - 1 + sp1->num_vertices)%sp1->num_vertices]].y);
	double x1 = abs( vertices[sp1->vertices[(v1_in_sp1 + 1)%sp1->num_vertices]].x - vertices[sp1->vertices[(v1_in_sp1 - 1 + sp1->num_vertices)%sp1->num_vertices]].x);
	double angle1 = atan2(y1, x1);

	double y2 = abs( vertices[sp2->vertices[(v2_in_sp2 + 1)%sp2->num_vertices]].y - vertices[sp2->vertices[(v2_in_sp2 - 1 + sp2->num_vertices)%sp2->num_vertices]].y);
	double x2 = abs( vertices[sp2->vertices[(v2_in_sp2 + 1)%sp2->num_vertices]].x - vertices[sp2->vertices[(v2_in_sp2 - 1 + sp2->num_vertices)%sp2->num_vertices]].x);
	double angle2 = atan2(y2, x2);

	double center_x = (v1->x + v2->x)*0.5;
	double center_y = (v1->y + v2->y)*0.5;
	double angle = 0.5*(angle1  + angle2);
	//Rotate angle 90ºC and calculate new positions (new edge will have a length of 1.5*T1_TRANSITION_CRITICAL_DISTANCE)
	v1->x = center_x + cos(angle)*LENGTH_ROTATED_EDGE;
	v2->x = center_x - cos(angle)*LENGTH_ROTATED_EDGE;

	v1->y= center_y + sin(angle)*LENGTH_ROTATED_EDGE;
	v2->y = center_y - sin(angle)*LENGTH_ROTATED_EDGE;

	edges[edge].length = distance(v1->ind, v2->ind);
}//end rotate_edge

void Tissue::t1_getCellRelationships(Vertex* v1, Vertex* v2, int& common_cell1, int& common_cell2, int& only_v1, int& only_v2){
	for(int i= 0; i < CELLS_PER_VERTEX; i++){
		if(contains(v1->cells[i], v2->cells, CELLS_PER_VERTEX)){
			if(common_cell1 == -1){
				common_cell1 = v1->cells[i];
			}else{
				common_cell2 = v1->cells[i];
			}
		}else{
			only_v1 = v1->cells[i];
			//only_v1_indv1 = i;
		}
	}
	//Get cell that is exclussive of v2
	for(int i= 0; i < CELLS_PER_VERTEX; i++){
		if(v2->cells[i] != common_cell1 && v2->cells[i] != common_cell2){
			only_v2 = v2->cells[i];
			//only_v2_indv2 = i;
		}
	}
}// end t1_getCellRelationships

//dist t1_inwards
double Tissue::t1_inwards_get_dist_sum(Vertex* v1, Vertex* v2, Cell* c1, Cell* s1, Cell* s2){

	int v1_in_cell1 = which(v1->ind, c1->vertices, c1->num_vertices);
	int v2_in_cell1 = which(v2->ind, c1->vertices, c1->num_vertices);
	int v1_in_s1 = which(v1->ind, s1->vertices, s1->num_vertices);
	int v2_in_s2 = which(v2->ind, s2->vertices, s2->num_vertices);

	//These integers will hold the two neighboring vertices (n1 and n2) of vertices v1 and v2, in c1 and c2 respectively, after transition, assuming c1 becomes specific to v1
	int v1_n1 = - 1, v1_n2 = - 1, v2_n1 = - 1, v2_n2 = - 1; 
	float dist;
	if(c1->vertices[(v1_in_cell1 + 1)%c1->num_vertices] == v2->ind){
		v1_n1 = c1->vertices[(v1_in_cell1 + 2)%c1->num_vertices];
		v1_n2 = c1->vertices[(v1_in_cell1 - 1 + c1->num_vertices)%c1->num_vertices];
	}else{
		v1_n1 = c1->vertices[(v1_in_cell1 + 1)%c1->num_vertices];
		v1_n2 = c1->vertices[(v1_in_cell1 - 2 + c1->num_vertices)%c1->num_vertices];
	}

	//get neighbours of v2 without using common_cell2
	if(s2->vertices[(v2_in_s2 + 1)%s2->num_vertices] == v1_n1 || s2->vertices[(v2_in_s2 + 1)%s2->num_vertices] == v1_n2){
		v2_n2 = s2->vertices[(v2_in_s2 - 1 + s2->num_vertices)%s2->num_vertices];
	}else{
		v2_n2 = s2->vertices[(v2_in_s2 + 1)%s2->num_vertices];
	}

	if(s1->vertices[(v1_in_s1 + 1)%s1->num_vertices] == v1_n1 || s1->vertices[(v1_in_s1 + 1)%s1->num_vertices] == v1_n2){
		v2_n1 = s1->vertices[(v1_in_s1 - 1 + s1->num_vertices)%s1->num_vertices];
	}else{
		v2_n1 = s1->vertices[(v1_in_s1 + 1)%s1->num_vertices];
	}

	dist = distance(v1->ind, v1_n1) + distance(v1->ind, v1_n2) + distance(v2->ind, v2_n1) + distance(v2->ind, v2_n2);

	StraightLine lv1_c1n1, lv1_c1n2, lv2_c2n1, lv2_c2n2, lv1v2, cell_edge;
	lv1_c1n1 = getLineFromEdge(v1, &this->vertices[v1_n1]);
	lv1_c1n2 = getLineFromEdge(v1, &this->vertices[v1_n2]);
	lv2_c2n1 = getLineFromEdge(v2, &this->vertices[v2_n1]); 
	lv2_c2n2 = getLineFromEdge(v2, &this->vertices[v2_n2]); 
	lv1v2 = getLineFromEdge(v1, v2);
	vector<StraightLine> lines = {lv1_c1n1, lv1_c1n2, lv2_c2n1, lv2_c2n2, lv1v2};
	vector<Cell*> neighbour_cells = {c1, s1, s2};
	int indices_avoid[6] = {v1->ind, v2->ind, v1_n1, v1_n2, v2_n1, v2_n2};
	
	Edge* ee;
	for(StraightLine l : lines){ 
		for(Cell* cellptr : neighbour_cells){
			for(int i = 0; i < cellptr->num_vertices; i++){
				ee = &this->edges[cellptr->edges[i]];
				//can't be one of the edges that are going to move, and edges cannot touch
				//(contains(s1->ind, ee->cells, 2) && contains(c2->ind, ee->cells, 2)) || (contains(s2->ind, ee->cells, 2) && contains(c1->ind, ee->cells, 2)) 
				if(!(contains(ee->vertices[0], indices_avoid, 4) && contains(ee->vertices[1], indices_avoid, 4)) && 
								!contains(l.v1, ee->vertices, 2) && !contains(l.v2, ee->vertices, 2) ){
					cell_edge=getLineFromEdge(ee);				
					if(lines_cross(cell_edge, l)) return -1;
				}
			}
		}
	}
	if( lines_cross(lv1_c1n1, lv2_c2n1) || lines_cross(lv1_c1n1, lv2_c2n2) || lines_cross(lv1_c1n2, lv2_c2n1) || lines_cross(lv1_c1n2, lv2_c2n2)) return -1;
	return dist;
}///dist t1_inwards

//c1 and c2 are cels common to v1 and v2 in original position. s1 and s2 are cells specific to v1 and v2, respectively, in original positions.
double Tissue::t1_get_dist_sum(Vertex* v1, Vertex* v2, Cell* c1 , Cell* c2, Cell* s1, Cell* s2){

	int v1_in_cell1 = which(v1->ind, c1->vertices, c1->num_vertices);
	int v1_in_cell2 = which(v1->ind, c2->vertices, c2->num_vertices);
	int v2_in_cell1 = which(v2->ind, c1->vertices, c1->num_vertices);
	int v2_in_cell2 = which(v2->ind, c2->vertices, c2->num_vertices);
	//These integers will hold the two neighboring vertices (n1 and n2) of vertices v1 and v2, in c1 and c2 respectively, after transition, assuming c1 becomes specific to v1
	int v1_n1 = - 1, v1_n2 = - 1, v2_n1 = - 1, v2_n2 = - 1; 

	float dist;
	if(c1->vertices[(v1_in_cell1 + 1)%c1->num_vertices] == v2->ind){
		v1_n1 = c1->vertices[(v1_in_cell1 + 2)%c1->num_vertices];
		v1_n2 = c1->vertices[(v1_in_cell1 - 1 + c1->num_vertices)%c1->num_vertices];
	}else{
		v1_n1 = c1->vertices[(v1_in_cell1 + 1)%c1->num_vertices];
		v1_n2 = c1->vertices[(v1_in_cell1 - 2 + c1->num_vertices)%c1->num_vertices];
	}

	if(c2->vertices[(v2_in_cell2 + 1)%c2->num_vertices] == v1->ind){
		v2_n1 =  c2->vertices[(v2_in_cell2 + 2)%c2->num_vertices];
		v2_n2 = c2->vertices[(v2_in_cell2 - 1 + c2->num_vertices)%c2->num_vertices];
	}else{
		v2_n1 = c2->vertices[(v2_in_cell2 + 1)%c2->num_vertices];
		v2_n2 = c2->vertices[(v2_in_cell2 - 2 + c2->num_vertices)%c2->num_vertices];
	}

	dist = distance(v1->ind, v1_n1) + distance(v1->ind, v1_n2) + distance(v2->ind, v2_n1) + distance(v2->ind, v2_n2);

	StraightLine lv1_c1n1, lv1_c1n2, lv2_c2n1, lv2_c2n2, lv1v2, cell_edge;
	lv1_c1n1 = getLineFromEdge(v1, &this->vertices[v1_n1]);
	lv1_c1n2 = getLineFromEdge(v1, &this->vertices[v1_n2]);
	lv2_c2n1 = getLineFromEdge(v2, &this->vertices[v2_n1]); 
	lv2_c2n2 = getLineFromEdge(v2, &this->vertices[v2_n2]); 
	lv1v2 = getLineFromEdge(v1, v2);
	vector<StraightLine> lines = {lv1_c1n1, lv1_c1n2, lv2_c2n1, lv2_c2n2, lv1v2};
	vector<Cell*> neighbour_cells = {c1, c2, s1, s2};
	int indices_avoid[6] = {v1->ind, v2->ind, v1_n1, v1_n2, v2_n1, v2_n2};
	
	Edge* ee;
	for(StraightLine l : lines){ 
		for(Cell* cellptr : neighbour_cells){
			for(int i = 0; i < cellptr->num_vertices; i++){
				ee = &this->edges[cellptr->edges[i]];
				//can't be one of the edges that are going to move, and edges cannot touch
				//(contains(s1->ind, ee->cells, 2) && contains(c2->ind, ee->cells, 2)) || (contains(s2->ind, ee->cells, 2) && contains(c1->ind, ee->cells, 2)) 
				if(!(contains(ee->vertices[0], indices_avoid, 6) && contains(ee->vertices[1], indices_avoid, 6)) && 
								!contains(l.v1, ee->vertices, 2) && !contains(l.v2, ee->vertices, 2) ){
					cell_edge=getLineFromEdge(ee);				
					if(lines_cross(cell_edge, l)) return -1;
				}
			}
		}
	}
	if( lines_cross(lv1_c1n1, lv2_c2n1) || lines_cross(lv1_c1n1, lv2_c2n2) || lines_cross(lv1_c1n2, lv2_c2n1) || lines_cross(lv1_c1n2, lv2_c2n2)) return -1;
	return dist;
}

//Assumes that number of cells contacting both vertices <=3:
void Tissue::make_join_limit_edges(Rearrangement& r){
	int edge = r.element_index;
	int aux; //Auxiliary variable used to store indices
	if(edges[edge].length > T1_TRANSITION_CRITICAL_DISTANCE) return; //Check that condition is still true (other rearrangements could have taken place since detection)

	Vertex* v1 = &vertices[ edges[edge].vertices[0] ];
	Vertex* v2 = &vertices[ edges[edge].vertices[1] ];
	//V1 will be preserved, v2 removed
	v1->x = (v1->x + v2->x)*0.5;
	v1->y = (v1->y + v2->y)*0.5;
	v2->dead = true;
	dead_vertices.push(v2->ind);
	this->num_vertices--;
	edges[edge].dead = true;
	dead_edges.push(edge);
	this->num_edges--;

	//1) Remove old edge from v1 and v2
	for(int i = which(edge, v1->edges, CELLS_PER_VERTEX); i < CELLS_PER_VERTEX; i++){
		v1->edges[i] = i < CELLS_PER_VERTEX - 1? v1->edges[i+1] : EMPTY_CONNECTION;
	}
	for(int i = which(edge, v2->edges, CELLS_PER_VERTEX); i < CELLS_PER_VERTEX; i++){
		v2->edges[i] = i < CELLS_PER_VERTEX - 1? v2->edges[i+1] : EMPTY_CONNECTION;
	}
	//2) Remove old vertex from v1.neighbour vertices
	for(int i = which(v2->ind, v1->neighbour_vertices, CELLS_PER_VERTEX); i < CELLS_PER_VERTEX; i++){
		v1->neighbour_vertices[i] = i < CELLS_PER_VERTEX - 1? v1->neighbour_vertices[i+1] : EMPTY_CONNECTION;
	}
	//3) Update cells 
	Cell* cell;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(v2->cells[i] != EMPTY_CONNECTION && !contains(v2->cells[i], v1->cells, CELLS_PER_VERTEX)){ //If cell i of vertex 2 does not touch vertex 1
			aux = v2->cells[i];
			cell = &this->cells[aux];
			for(int j = 0; j < CELLS_PER_VERTEX; j++){
				if(v1->cells[j] == EMPTY_CONNECTION){
					v1->cells[j] = v2->cells[i]; //Add cell i of v2 in the first empty location of v1 cell vector
					break;
				}//end if v1 location is empty
			}//end for v1
		}//end if cell connected to v2 is not in v1
	}//end for cells in v2

	//4) Update edges and vertices in edges
	Edge* ee;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(v2->edges[i] != EMPTY_CONNECTION && !contains(v2->edges[i], v1->edges, CELLS_PER_VERTEX)){ //If edge i of vertex 2 does not touch vertex 1
			aux = v2->edges[i];
			ee = &this->edges[aux];
			for(int j = 0; j < CELLS_PER_VERTEX; j++){
				if(v1->edges[j] == EMPTY_CONNECTION){
					v1->edges[j] = v2->edges[i]; //Add edge i of v2 in the first empty location of v1 edge vector
					break;
				}//end if v1 location is empty
			}//end for v1
			ee->vertices[ which(v2->ind, ee->vertices, MAX_SIDES_PER_CELL) ] = v1->ind; //add index of v1 to cell, in the place where v2 was
		}//end if edge connected to v2 is not in v1
	}//end for edges in v2

	//5) Update neighbour vertices
	Vertex* vv;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(v2->neighbour_vertices[i] != EMPTY_CONNECTION && v2->neighbour_vertices[i] != v1->ind && !contains(v2->neighbour_vertices[i], v1->neighbour_vertices, CELLS_PER_VERTEX)){
			aux = v2->neighbour_vertices[i];
			vv = &this->vertices[aux];
			for(int j = 0; j < CELLS_PER_VERTEX; j++){
				if(v1->neighbour_vertices[j] == EMPTY_CONNECTION){
					v1->neighbour_vertices[j] = v2->neighbour_vertices[i]; //Add cell i of v2 in the first empty location of v1 cell vector
					break;
				}//end if v1 location is empty
			}//end for v1
			vv->neighbour_vertices[ which(v2->ind, vv->neighbour_vertices, CELLS_PER_VERTEX) ] = v1->ind;
		}
	}


	//6) Update vertices in cells
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(v2->cells[i] != EMPTY_CONNECTION){
			aux = v2->cells[i];
			cell = &this->cells[aux];
			//Now remove old vertex while retaining order
			int j = 0;
			while(cell->vertices[j]  != v1->ind && cell->vertices[j] != v2->ind) j++;
			cell->vertices[j]  = v1->ind;
			for(int k = j + 1; k < MAX_SIDES_PER_CELL; k++){
				cell->vertices[k] = k < MAX_SIDES_PER_CELL - 1? cell->vertices[k + 1] : EMPTY_CONNECTION;
			}
		}
	}
	//7) Remove old edge from cells
	for(int i = 0; i <  VERTEX_PER_EDGE; i++){
		if(edges[edge].cells[i] != EMPTY_CONNECTION){
			aux = edges[edge].cells[i];
			cell = &this->cells[aux];
			for(int j = which( edge, cell->edges, MAX_SIDES_PER_CELL); j < MAX_SIDES_PER_CELL; j++){
				cell->edges[j] = j < (MAX_SIDES_PER_CELL-1)? cell->edges[j + 1] : EMPTY_CONNECTION;
			}
			cell->num_vertices--;
		}

	}


	//8) Update edge lengths
	for(int i = 0; i< CELLS_PER_VERTEX; i++){
		if(v1->edges[i] != EMPTY_CONNECTION) edges[v1->edges[i]].length = distance(edges[v1->edges[i]].vertices[0], edges[v1->edges[i]].vertices[1]);	
	}

	//9) Update cell areas and perimeter
	for(int i = 0; i< CELLS_PER_VERTEX; i++){
		if(v1->cells[i] != EMPTY_CONNECTION){
			aux = v1->cells[i];
			cell = &this->cells[aux];
			cell->area = calculateCellArea(*cell);
			cell->perimeter = calculateCellPerimeter(*cell);
		}		
	}

	//10) Energies will be updated when moving vertices, not now
	
	cout << "REMOVED: " << v2->ind << " and replaced by: " << v1->ind << endl;
	

}

void Tissue::make_t2(Rearrangement& r){

	int cell = r.element_index;
	//cout << "A cell: " << cell << " dead: " << cells[cell].dead << endl;
	if(cells[cell].area > T2_TRANSITION_CRITICAL_AREA || cells[cell].num_vertices > 3 || cells[cell].dead) return; 

	int v1 = cells[cell].vertices[0], v2 = cells[cell].vertices[1], v3 = cells[cell].vertices[2];
	vertices[v2].dead = true;
	vertices[v3].dead = true;
	vertices[v1].x = (vertices[v1].x + vertices[v2].x + vertices[v3].x)/3;
	vertices[v1].y = (vertices[v1].y + vertices[v2].y + vertices[v3].y)/3;
	
	//cout << "B" << endl;
	//Remove edges between vertices of dead cell (still not remove; will be replaced later in v1 edges)
	int v2e = -1, v3e = -1, v23e = -1;
	cells[cell].dead = true;
	for(int i = 0; i < cells[cell].num_vertices; i++){
		if(contains(v1, edges[cells[cell].edges[i]].vertices, 2)){
			if(contains(v2, edges[cells[cell].edges[i]].vertices, 2)){
				v2e = cells[cell].edges[i];
			}else if(contains(v3, edges[cells[cell].edges[i]].vertices, 2)){
				v3e = cells[cell].edges[i];
			}
		}else{ //if(contains(v2, edges[cells[cell].edges[i]].vertices, 2) && contains(v3, edges[cells[cell].edges[i]].vertices, 2)){
			v23e = cells[cell].edges[i];
		}
		edges[cells[cell].edges[i]].dead=true;
	}

	//cout << "C" << endl;
	//Replace neighbours
	int v2nei = -1, v3nei = -1;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		v2nei = vertices[v2].neighbour_vertices[i];
		if(v2nei == EMPTY_CONNECTION){
			vertices[v1].neighbour_vertices[which(v2, vertices[v1].neighbour_vertices, CELLS_PER_VERTEX)] = EMPTY_CONNECTION;
			break;
		}else if(v2nei != v1  && v2nei != v3){
			vertices[v1].neighbour_vertices[which(v2, vertices[v1].neighbour_vertices, CELLS_PER_VERTEX)] = vertices[v2].neighbour_vertices[i]; 
			vertices[v2nei].neighbour_vertices[which(v2, vertices[v2nei].neighbour_vertices, CELLS_PER_VERTEX)] = v1;
			break;
		}
	}
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		v3nei = vertices[v3].neighbour_vertices[i];
		if(v3nei == EMPTY_CONNECTION){
			vertices[v1].neighbour_vertices[which(v3, vertices[v1].neighbour_vertices, CELLS_PER_VERTEX)] = EMPTY_CONNECTION;
			break;
		}else if(v3nei != v1  && v3nei != v2){
			vertices[v1].neighbour_vertices[which(v3, vertices[v1].neighbour_vertices, CELLS_PER_VERTEX)] = vertices[v3].neighbour_vertices[i]; 
			vertices[v3nei].neighbour_vertices[which(v3, vertices[v3nei].neighbour_vertices, CELLS_PER_VERTEX)] = v1;
			break;
		}
	}



	//cout << "E" d< endl;
	//change edges connecting neighbours to v2 and v3
	int echange2, echange3;
	if(v2nei != EMPTY_CONNECTION){
		for(int i = 0; i < CELLS_PER_VERTEX; i++){
			if( contains(v2, edges[vertices[v2nei].edges[i]].vertices, 2) ){
				echange2 = vertices[v2nei].edges[i];
				edges[echange2].vertices[which(v2, edges[echange2].vertices, 2)] = v1;
				edges[echange2].length = distance(v1, v2nei);
				vertices[v1].edges[which(v2e, vertices[v1].edges, CELLS_PER_VERTEX)] = echange2;
				cout << "echange2: " << echange2 << " v1: " << v1 << " v2 " << v2 << " v2nei: " << v2nei << " vertices of echange2: " << edges[echange2].vertices[0] << " " << edges[echange2].vertices[1] << endl;
				break;
			}
		}
	}else{
		vertices[v1].edges[which(v2e, vertices[v1].edges, CELLS_PER_VERTEX)] = EMPTY_CONNECTION;
	}
	if(v3nei != EMPTY_CONNECTION){
		for(int i = 0; i < CELLS_PER_VERTEX; i++){
			if( contains(v3, edges[vertices[v3nei].edges[i]].vertices, 2) ){
				echange3 = vertices[v3nei].edges[i];
				edges[echange3].vertices[which(v3, edges[echange3].vertices, 2)] = v1;
				edges[echange3].length = distance(v1, v3nei);
				vertices[v1].edges[which(v3e, vertices[v1].edges, CELLS_PER_VERTEX)] = echange3;
				break;
			}
		}
	}else{
		vertices[v1].edges[which(v3e, vertices[v1].edges, CELLS_PER_VERTEX)] = EMPTY_CONNECTION;
	}
	//cout << "F" << endl;

	//Recalculate distance of v1 neighbour outside of dying cell
	int v1e;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		v1e = vertices[v1].edges[i];
		if(v1e != EMPTY_CONNECTION){
			edges[v1e].length = distance(edges[v1e].vertices[0], edges[v1e].vertices[1]);
		}

	}

	int neicell2, neicell3, neicell23;
	neicell2 = edges[v2e].cells[0] == cell ? edges[v2e].cells[1] : edges[v2e].cells[0];
	neicell3 = edges[v3e].cells[0] == cell ? edges[v3e].cells[1] : edges[v3e].cells[0];
	neicell23 = edges[v23e].cells[0] == cell ? edges[v23e].cells[1] : edges[v23e].cells[0];

	//cout << "G v2: " << v2 << "v3:  " << v3 << " v2e: "  << v2e << " v3e: " << v3e << " v23e: " << v23e << " cell: " << cell << " movs. accepted: " << counter_moves_accepted <<endl;
	//Change vertices in cells
	if(neicell2 != EMPTY_CONNECTION){
		removeConnectionCell(v2e, cells[neicell2].edges, cells[neicell2].num_vertices);
		removeConnectionCell(v2, cells[neicell2].vertices, cells[neicell2].num_vertices);
		cells[neicell2].num_vertices--;
	}
	if(neicell3 != EMPTY_CONNECTION){
		//cout << "cell " << cell << " neicell3 " << neicell3 << " v3 " << v3 << endl;
		removeConnectionCell(v3e, cells[neicell3].edges, cells[neicell3].num_vertices);
		removeConnectionCell(v3, cells[neicell3].vertices, cells[neicell3].num_vertices);
		cells[neicell3].num_vertices--;
	}
	if(neicell23 != EMPTY_CONNECTION && neicell23 != neicell3 && neicell23 != neicell2){
		removeConnectionCell(v23e, cells[neicell23].edges, cells[neicell23].num_vertices);
		removeConnectionCell(v3, cells[neicell23].vertices, cells[neicell23].num_vertices);
		cells[neicell23].num_vertices--;
		cells[neicell23].vertices[which(v2, cells[neicell23].vertices, cells[neicell23].num_vertices)] = v1;
	} else if(neicell23 == neicell3 && neicell23 != EMPTY_CONNECTION){
		removeConnectionCell(v23e, cells[neicell23].edges, cells[neicell23].num_vertices);
		cells[neicell23].vertices[which(v2, cells[neicell23].vertices, cells[neicell23].num_vertices)] = v1;
	}else if(neicell23 == neicell2 && neicell23 != EMPTY_CONNECTION){
		removeConnectionCell(v23e, cells[neicell23].edges, cells[neicell23].num_vertices);
		cells[neicell23].vertices[which(v3, cells[neicell23].vertices, cells[neicell23].num_vertices)] = v1;
	}


	//cout << "H" << endl;
	//Change cell in v1
	vertices[v1].cells[which(cell, vertices[v1].cells, CELLS_PER_VERTEX)] = neicell23;

	//Recalculate areas and perimeters in neighboring cells
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(vertices[v1].cells[i] != EMPTY_CONNECTION){
			cells[vertices[v1].cells[i]].area = calculateCellArea(cells[vertices[v1].cells[i]]);
			cells[vertices[v1].cells[i]].perimeter = calculateCellPerimeter(cells[vertices[v1].cells[i]]);
		}
	}
	//cout << "I" << endl;
	//Add dead elements to dead element queues:
	dead_cells.push(cell);
	dead_vertices.push(v1);
	dead_vertices.push(v2);
	dead_edges.push(v2e);
	dead_edges.push(v3e);
	dead_edges.push(v23e);

	this->num_cells--;
	this->num_vertices-=2;
	this->num_edges-=3;

	this->counter_t2++;

	cout << "T2 transition: cell=" << cell << ", survivor vertex =" << v1 << "; v2: " << v2 << "; v3: " << v3 << "; v2nei: " << v2nei << "; v3nei: " << v3nei << "; mov. accepted: " << this->counter_moves_accepted << "; T1: " << counter_t1 << "; T1 abortions: " << counter_t1_abortions << "; T2: " << counter_t2 << endl;

}//End make transition T2


//Assumes that array elements is at least 1 element longer than length
void Tissue::removeConnectionCell(int elm, int* elements, int length){
	int i = 0;
		while(elements[i] != elm) i++;
		do{
			elements[i] = elements[i+1];
			i++;
		}while(i < length - 1); 
	elements[length - 1] = EMPTY_CONNECTION;
}//remove element from array of cell connections

bool Tissue::check_if_edges_cross(int vertex){
	Vertex *v = &this->vertices[vertex];
	Edge *e1, *e2;
	Cell *c1;
	StraightLine l_e1, l_e2;

	for(int i = 0; i < CELLS_PER_VERTEX; i++){ //For each edge i touching the vertex of interest
		if(v->edges[i] != EMPTY_CONNECTION){
			e1 = &this->edges[v->edges[i]];
			l_e1 = getLineFromEdge(e1);
			for(int j = 0; j < CELLS_PER_VERTEX; j++){ //For each cell j touching the vertex of interest
				if(v->cells[j] != EMPTY_CONNECTION){
					c1 = &this->cells[v->cells[j]];
					for(int k = 0; k < c1->num_vertices; k++){ //For each edge k in cell j
						e2 = &this->edges[c1->edges[k]];
						if(!contains(e1->vertices[0], e2->vertices, 2) && !contains(e1->vertices[1], e2->vertices, 2)){//edges i and j can't touch 
							l_e2 = getLineFromEdge(e2);
							if(lines_cross(l_e1, l_e2)) return true; //Edges cross			
						}
					}
				}
			}

		}
	}

	return false; //Edges do not cross
}

StraightLine Tissue::getLineFromEdge(const Edge* e){
	Vertex *v1 = &this->vertices[e->vertices[0]];
	Vertex *v2 = &this->vertices[e->vertices[1]];
	StraightLine sl;
	sl.x1 = v1->x;
	sl.x2 = v2->x;
	sl.y1 = v1->y;
	sl.y2 = v2->y;
	sl.slope = (sl.y1 - sl.y2)/(sl.x1 - sl.x2);
	sl.intercept = sl.y2 - sl.x2*sl.slope;
	sl.v1 = v1->ind;
	sl.v2 = v2->ind;
	return sl;
}

StraightLine Tissue::getLineFromEdge(const Vertex* v1, const Vertex* v2){
	StraightLine sl;
	sl.x1 = v1->x;
	sl.x2 = v2->x;
	sl.y1 = v1->y;
	sl.y2 = v2->y;
	sl.slope = (sl.y1 - sl.y2)/(sl.x1 - sl.x2);
	sl.intercept = sl.y2 - sl.x2*sl.slope;
	sl.v1 = v1->ind;
	sl.v2 = v2->ind;
	return sl;
}

StraightLine Tissue::getLineFromEdge(double x1, double x2, double y1, double y2){
	StraightLine sl;
	sl.x1 = x1;
	sl.x2 = x2;
	sl.y1 = y1;
	sl.y2 = y2;
	sl.slope = (sl.y1 - sl.y2)/(sl.x1 - sl.x2);
	sl.intercept = sl.y2 - sl.x2*sl.slope;
	return sl;
}

bool Tissue::lines_cross(StraightLine& a, StraightLine& b){
	//if(a.v1 == b.v1 || a.v1 == b.v2 || a.v2 == b.v1 || a.v2 == b.v2) return false;
	double x_intersect = (b.intercept - a.intercept)/(a.slope - b.slope);
	return ((x_intersect < a.x1 && x_intersect > a.x2) || (x_intersect < a.x2 && x_intersect > a.x1) ) && ((x_intersect < b.x1 && x_intersect > b.x2) || (x_intersect < b.x2 && x_intersect > b.x1));
	
}

std::string Tissue::getStats(){
	std::string s = "";
	s += "Move trials: " + std::to_string(counter_move_trials);
	s += ", Moves accepted: " + std::to_string(counter_moves_accepted);
	s += ", prop. accepted: " + std::to_string(counter_moves_accepted/float(counter_move_trials));
	s += ", T1 accepted: " + std::to_string(counter_t1);
	s += ", T1 rejected: "  + std::to_string(counter_t1_abortions);
	return s;
}

int Tissue::getCounterT1(){
	return this->counter_t1;
}

void Tissue::writeCellsFile(std::string fname){
	fname = fname + CELLS_FILE_EXTENSION;
	ofstream fo;
	fo.open(fname);
	fo << num_cells << "\t" << MAX_SIDES_PER_CELL << "\n";
	for(Cell c: cells){
		if(! c.dead){
			for (int v = 0; v < MAX_SIDES_PER_CELL; v++){
				fo << c.vertices[v] << "\t";
			}
			fo << "\n";
		}
	}
	fo.close();
}


void Tissue::writePointsFile(std::string fname){
	fname = fname + VERTEX_FILE_EXTENSION;
	ofstream fo;
	fo.open(fname);
	fo << num_vertices << "\n";
	for(Vertex v: vertices){
		if(!v.dead) fo << v.x << "\t" << v.y << "\t" << v.ind <<"\n";
	}
	fo.close();
}

void Tissue::writeAllData(std::string fname){ //Writes tables with all data (including parameters for each element)
	ofstream of;
	of.open(fname + ".out");
	of << *this;
	of.close();
}//End  writeAll


///////////////////////////////////////////////////////////////////////////////////////////////
//What follows is only useful to print easily, not part of the model
//
//Overloading of << operator for VERTICES. Useful to print
std::ostream& operator<<(std::ostream& out, const Vertex& v){
	out << v.ind << "\t" << v.x << "\t" << v.y << "\t" << v.energy << "\t";
	for(int i = 0; i < sizeof(v.cells)/sizeof(v.cells[0]); i++){ // print cells touching vertex separated by ','
		out << v.cells[i] << ",";
	} 
	out << "\t";
	for(int i = 0; i < sizeof(v.edges)/sizeof(v.edges[0]); i++){ // print edges touching vertex separated by ','
		out << v.edges[i] << ",";
	} 
	out << "\t";
	for(int i = 0; i < sizeof(v.neighbour_vertices)/sizeof(v.neighbour_vertices[0]); i++){ // print neighbour vertices separated by ','
		out << v.neighbour_vertices[i] << ",";
	} 

	return out;
}


//Overloading of << operator for CELLS. Useful to print
std::ostream& operator<<(std::ostream& out, const Cell& c){
	out << c.ind << "\t" << int(c.type) << "\t" << c.area << "\t" << c.preferred_area << "\t" << c.perimeter << "\t" << c.perimeter_contractility <<"\t"<< c.num_vertices << "\t";
	for(int i = 0; i < c.num_vertices; i++){  // print vertices of cell separated by ','
		out << c.vertices[i] << ",";
	} 
	out << "\t";
	for(int i = 0; i < c.num_vertices; i++){ // print edges of cell separated by ','
		out << c.edges[i] << ",";
	} 
	return out;
}

//Overloading of << operator for EDGES. Useful to print
std::ostream& operator<<(std::ostream& out, const Edge& e){
	out << e.ind << "\t" << int(e.type)  << "\t" << e.length << "\t" <<e.tension << "\t";
	for(int i = 0; i < sizeof(e.vertices)/sizeof(e.vertices[0]); i++){  //print vertices touching edge separated by ','
		out << e.vertices[i] << ",";
	} 
	out << "\t";
	for(int i = 0; i < sizeof(e.cells)/sizeof(e.cells[0]); i++){  //print cells touching edge separated by ','
		out << e.cells[i] << ",";
	} 
	return out;
}

//Overloading of << operator for TISSUE. Useful to print
std::ostream& operator<<(std::ostream& out, const Tissue& t){
	out << t.simname << "\t-\t" << "max. accepted moves: " << t.max_accepted_movements << ",\twrite every N moves: " << t.write_every_N_moves <<"\n\n";
	std::string s = "";
	s += "Move trials: " + std::to_string(t.counter_move_trials);
	s += ", Moves accepted: " + std::to_string(t.counter_moves_accepted);
	s += ", prop. accepted: " + std::to_string(t.counter_moves_accepted/float(t.counter_move_trials));
	s += ", T1 accepted: " + std::to_string(t.counter_t1);
	s += ", T1 rejected: "  + std::to_string(t.counter_t1_abortions);
	s+= ", T2: " + std::to_string(t.counter_t2);
	out << s << "\n";
	//out << t.getStats() << "\n";
	out << "VERTICES:\n";
	out << VERTEX_HEADER;
	for(Vertex v: t.vertices){
		if(!v.dead) out << v << "\n";
	}
	out << "\nCELLS:\n";
	out << CELL_HEADER;
	for(Cell c:t.cells){
		if(!c.dead) out << c << "\n";
	}
	out << "\nEDGES:\n";
	out << EDGE_HEADER;
	for(Edge e:t.edges){
		if(!e.dead) out << e << "\n";
	}
	out << "****\n****\n";
	return out;
}




//Helper functions of general use.
//Returns true if integer n is present in array a, otherwise returns false.
bool contains(int n, const int* a, int len){
	for(int i = 0; i<len; i++){
		if(n == a[i]) return true;
	}
	return false;
}
//Returns -1 if integer n is not in array a, otherwise returns index of n in a.
int which(int n, const int* a, int len){
	for(int i = 0; i<len; i++){
		if(n == a[i]) return i;
	}
	return -1;
}

//Returns index of FIRST element in array a that is not present in array b
int element_not_in(const int* a, const int* b, int len1, int len2){
	bool a_in_b;
	for(int i = 0; i<len1; i++){
		a_in_b = false;
		for(int j = 0; j<len2; j++){
			if(a[i] == b[j]){
				a_in_b = true;
				break;
			}
		}
		if(! a_in_b) return i;
	}
	return -1;
}



