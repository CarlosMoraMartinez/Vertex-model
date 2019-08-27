//VertexSystem.h

#include <iostream>
#include <vector>
#include <random>
#include <cstdlib>
#include <cmath>
#include <string>
#include <queue>

// Default model constants. From Ray et al. (2015) Dev. Cell.
const float MIN_RANGE_VERTEX_MOVEMENT = 0; //not used so far
const float MAX_RANGE_VERTEX_MOVEMENT = 0.02;
const float TEMPERATURE_POSITIVE_ENERGY = 0.005;
const float TEMPERATURE_NEGATIVE_ENERGY = 0.1; //Acceptance probability is basically 1
//const float PROB_ACCEPT_NON_FAVORABLE = 0.05;
const float LINE_TENSION_BLADE = -0.02;
const float LINE_TENSION_HINGE = -0.01;
const float LINE_TENSION_VEIN_HINGE = 7.0;
const float LINE_TENSION_VEIN_BLADE = 1.0;
const float LINE_TENSION_TISSUE_BOUNDARY = 4.0;//4.0;
const float SPRING_CONSTANT = 2;
const float PERIMETER_CONTRACT_BLADE = 0.03;
const float PERIMETER_CONTRACT_HINGE = 0.6;
const float T1_TRANSITION_CRITICAL_DISTANCE = 0.1; // 0.1;
const float T2_TRANSITION_CRITICAL_AREA = 0.01;//0.9;//0.01;
const float MAX_CELL_AREA = 30;//30;		
const float MAX_EDGE_LENGTH = 10; //Not used yet; not in paper
const float PREFERRED_AREA_INITIAL = 35.0; //Use >= 30 for constant division
const float PREFERRED_AREA_HINGE = 3.0; //Use >= 30 for constant division
const float PREFERRED_AREA_FINAL = 5.0;
const float DIVISION_ANGLE_RANDOM_NOISE = 0.3; //1 = variation of 360ºC (completely random), 0 =no variation (division alwys orthogonal to the max. lengh)
const float LENGTH_ROTATED_EDGE = 0.5*T1_TRANSITION_CRITICAL_DISTANCE*1.2; //after a t1 transition, rotated edge length is multiplied by twice this constant


// Implementation constants
const unsigned short CELLS_PER_VERTEX = 3;
const unsigned short MAX_SIDES_PER_CELL = 20; //Maybe pass as argument instead of constant?
const unsigned short VERTEX_PER_EDGE = 2;
const int EMPTY_CONNECTION = -999;  //Value to initialize arrays (cells in vertex, etc)
const int RANDOM_SEED = 1234;

const std::string VERTEX_FILE_EXTENSION = ".points";
const std::string CELLS_FILE_EXTENSION = ".cells";
const std::string SPRING_FILE_EXTENSION = ".spr";

const bool T1_ACTIVE = true;
const bool T1_BORDER_INWARDS_ACTIVE = true;
const bool T1_BORDER_OUTWARDS_ACTIVE = true;
const bool DIVISION_ACTIVE = true;
const bool T2_ACTIVE = true;
const bool JOIN_EDGES_ACTIVE = true;
const bool CONTROL_CELLS_2SIDES = true; 
const int MOVE_TRIALS = 100;  //Times it tries to move a vertex before it quids because always makes edges to cross

const bool REPORT_T1 = false;
const bool REPORT_DIV = false;

const std::string VERTEX_HEADER = "ind\tx\ty\tenergy\tmovable\tspring\tcells\tedges\tneighbour_vertices\n";
const std::string CELL_HEADER = "ind\ttype\tarea\tpreferred_area\tperimeter\tperim_contract\tnum_vertices\tvertices\tedges\n";
const std::string EDGE_HEADER = "ind\ttype\tlength\ttension\tvertices\tcells\n";

//Enum class to define types of cells 
enum class CellType{blade = 0, hinge = 1, vein = 2, vein_hinge = 3};

//Enum class to define types of edges
enum class EdgeType{blade = 0, hinge = 1, vein_hinge = 2, vein_blade = 3, tissue_boundary = 4, spring = 5, vein = 4};

//Enum class to define types of transitions
enum class RearrangementType{t1 = 0, t2 = 1, divide_cell = 2, divide_edge = 3, join_limit_edges = 4, t1_at_border_outwards = 5, t1_at_border_inwards = 6, rotate_inwards = 7};
// Structure for Vertex
struct Vertex{
	int ind;
	double x;
	double y;
	double energy;
	int cells[CELLS_PER_VERTEX];
	int edges[CELLS_PER_VERTEX];
	int neighbour_vertices[CELLS_PER_VERTEX];
	bool movable;
	int spring;
	bool dead;
};
// Structure for Cell
struct Cell{
	int ind;
	CellType type;
	double area;
	double preferred_area;
	double perimeter;
	double perimeter_contractility;
	int edges[MAX_SIDES_PER_CELL];
	int vertices[MAX_SIDES_PER_CELL];
	int num_vertices; 
	bool dead;
};
// Structure for Edge
struct Edge{
	int ind;
	EdgeType type;
	double tension;
	double length;
	int vertices[VERTEX_PER_EDGE];
	int cells[VERTEX_PER_EDGE];
	bool can_transition;
	bool dead;
};

struct Rearrangement{
	int element_index;
	RearrangementType type;
};

struct StraightLine{
	double intercept;
	double slope;
	double x1;
	double x2;
	double y1;
	double y2;
	int v1;
	int v2;
};

//Vector types for each structure
typedef std::vector<Vertex> vertex_v;
typedef std::vector<Cell> cell_v;
typedef std::vector<Edge> edge_v;
typedef std::queue<Rearrangement> rearrangement_q;

class Tissue{
	public:
		Tissue();
		Tissue(std::string starting_tissue_file, int max_accepted_movements,  int write_every_N_moves=1000);
		Tissue(std::string starting_tissue_file, std::string params_file, int max_accepted_movements,  int write_every_N_moves=1000); //not implemented

		void simulate();
		double calculateCellArea(const Cell& c);
		double calculateCellPerimeter(const Cell& c);
		double distance(int v1, int v2);
		double calculateEnergy(Vertex& v);
		void moveVertex(Vertex& v, float x, float y);
		bool tryMoveVertex(std::default_random_engine& generator, std::uniform_real_distribution<double>& unif);
		void detectChangesAfterMove(int vertex_moved);
		void performRearrangements();
		int newVertex();
		int newVertex(double x, double y);
		int newCell();
		int newEdge();

		void make_t1(Rearrangement& r);
		void make_t1_at_border_inwards(Rearrangement& r); 
		void make_t1_at_border_outwards(Rearrangement& r);
		void make_t2(Rearrangement& r);
		void make_divide_cell(Rearrangement& r);
		void make_edge_division(Rearrangement& r); //this will be removed
		void make_join_limit_edges(Rearrangement& r); //this will be removed
		
		int addVertex(Vertex v); //Only used for testing
		int addCell(Cell c); //Only used for testing
		int addEdge(Edge e); //Only used for testing

		void writeCellsFile(std::string fname);
		void writePointsFile(std::string fname);
		void writeSpringsFile(std::string fname);
		void writeAllData(std::string fname);
		void produceOutputs(std::string add_to_name);
		std::string getStats();
		int getCounterT1();
		
		friend std::ostream& operator<<(std::ostream& out, const Tissue& t); //Don't worry about these 4 lines, only important for output
		friend std::ostream& operator<<(std::ostream& out, const Vertex& v);
		friend std::ostream& operator<<(std::ostream& out, const Cell& c);
		friend std::ostream& operator<<(std::ostream& out, const Edge& e);
		
	private:
		std::string simname;
		int num_cells, num_vertices, num_edges, num_springs;

		//parameters
		float min_range_vertex_movement; 
		float max_range_vertex_movement;
		float temperature_positive_energy;
		float temperature_negative_energy; 
		float line_tension_blade;
		float line_tension_hinge;
		float line_tension_vein_hinge;
		float line_tension_vein_blade;
		float line_tension_tissue_boundary;
		float spring_constant;
		float perimeter_contract_blade;
		float perimeter_contract_hinge;
		float t1_transition_critical_distance; 
		float t2_transition_critical_area;
		float max_cell_area;	
		float max_edge_length; 
		float preferred_area_initial;
		float preferred_area_final;
		float preferred_area_hinge;
		float division_angle_random_noise; 
		float length_rotated_edge; 

		//Data structures
		vertex_v vertices;
		cell_v cells;
		edge_v edges;
		edge_v springs;
		std::queue<int> dead_vertices;
		std::queue<int> dead_cells;
		std::queue<int> dead_edges;	
		rearrangement_q rearrangements_needed;

		//counters
		int counter_move_trials, counter_moves_accepted, counter_edges_removed, counter_t1, counter_t1_abortions, counter_divisions, counter_t2, counter_t1_outwards, counter_t1_inwards;
		int max_accepted_movements, write_every_N_moves;

		//Methods used by constructors
		void initialize_vertices(std::ifstream& inp);
		void initialize_cells(std::ifstream& inp);  
		void initialize_edges();
		void initialize_springs(std::ifstream&);
		void addEdge(int v1, int v2, int cell);		//Creates single edge
		void set_default_params();			//Sets model parameters for each vertex, cell and edge, from constants defined in this file
		void setEdgeType(int e);
		void addCellToVertex(int vertex, int cell); 	//Adds one cell to the array of cells in a vertex structure
		void addNeighbourVertex(int vertex1, int vertex2);
		void addEdgeToVertex(int vertex, int edge); 	//Adds one edge to the array of edges in a vertex structure
		void addEdgeToCell(int cell, int edge); 	//Adds one edge to the array of edges in a cell structure

		//Methods related to detect edge crossings
		bool check_if_edges_cross(int vertex);
		StraightLine getLineFromEdge(const Edge* e);
		StraightLine getLineFromEdge(const Vertex* v1, const Vertex* v2);
		StraightLine getLineFromEdge(double x1, double x2, double y1, double y2);
		bool lines_cross(StraightLine &a, StraightLine &b);

		//These methods are different steps in t1 transition
		double t1_get_dist_sum(Vertex* v1, Vertex* v2, Cell* c1, Cell* c2, Cell* s1, Cell* s2);
		void t1_rotate_edge90degrees(Vertex* v1, Vertex* v2, int edge); //not used anymore
		void t1_rotate_edge(Vertex* v1, Vertex* v2, int edge, Cell* c1, Cell* c2, Cell*sp1, Cell* sp2);
		void t1_getCellRelationships(Vertex* v1, Vertex* v2, int &common_cell1, int &common_cell2, int &only_v1, int &only_v2);
		void t1_update_neighbours(Vertex* v1, Vertex* v2, int edge, int common_cell1, int common_cell2, int only_v1, int only_v2, int &remove_from_v1, int &remove_from_v2);
		void t1_update_edges(Vertex* v1, Vertex* v2, int edge, int remove_from_v1, int remove_from_v2, int &e_remove_from_v1, int &e_remove_from_v2);
		void t1_add_vertices_to_cells(Vertex* v1, Vertex* v2, Cell* sp2, int remove_from_v2);
		void t1_update_sizes(Vertex* v1, Vertex* v2, int edge);


		//Methods used by T1 transitions in borders
		double t1_inwards_get_dist_sum(Vertex* v1, Vertex* v2, Cell* c1, Cell* s1, Cell* s2);
		void t1_inwards_update_neighbours(Vertex* v1, Vertex* v2, int edge, int common_cell1, int only_v1, int only_v2, int& remove_from_v1, int& remove_from_v2);

		//Methods used by cell division
		void splitEdgeWithVertex(int e, int cell, int  v);

		//methods used by T2 and join_edges
		void removeConnectionCell(int elm, int* elements, int length);
		void make_remove_size2cell(int cell);
};

//functions of general use
bool contains(int n, const int* a, int len);
int which(int n, const int* a, int len);
int count(int n, const int* a, int len);
int element_not_in(const int* a, const int* b, int len1, int len2);

//What follows is only useful to print easily, not part of the model
//
//Overloading of << operator for VERTICES. Useful to print
std::ostream& operator<<(std::ostream& out, const Vertex& v);

//Overloading of << operator for CELLS. Useful to print
std::ostream& operator<<(std::ostream& out, const Cell& c);

//Overloading of << operator for EDGES. Useful to print
std::ostream& operator<<(std::ostream& out, const Edge& e);

//Overloading of << operator for TISSUE. Useful to print
std::ostream& operator<<(std::ostream& out, const Tissue& t);


