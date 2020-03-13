//VertexSystem.h

#include <iostream>
#include <vector>
#include <map>
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
const bool TEMPERATURE_MEANS_PROPORTION_OF_ACCEPTANCE = true;
//const float PROB_ACCEPT_NON_FAVORABLE = 0.05;
const float LINE_TENSION_BLADE = -0.02;
const float LINE_TENSION_HINGE = -0.01;
const float LINE_TENSION_VEIN_HINGE = 7.0;
const float LINE_TENSION_VEIN_BLADE = 1.0;
const float LINE_TENSION_TISSUE_BOUNDARY = 4.0;//4.0;
const float SPRING_CONSTANT = 2;
const float PERIMETER_CONTRACT_BLADE = 0.03;
const float PERIMETER_CONTRACT_HINGE = 0.6;
const float T1_TRANSITION_CRITICAL_DISTANCE = 0.01; // 0.1;
const float T2_TRANSITION_CRITICAL_AREA = 0.01;//0.9;//0.01;
const float MAX_CELL_AREA = 20;//30;		
const float MAX_EDGE_LENGTH = 10; //Not used yet; not in paper
const float PREFERRED_AREA_INITIAL = 10.0; //Use >= 30 for constant division
const float PREFERRED_AREA_HINGE = 3.0; //Use >= 30 for constant division
const float PREFERRED_AREA_FINAL = 5.0;
const float DIVISION_ANGLE_RANDOM_NOISE = 0.3; //1 = variation of 360ºC (completely random), 0 =no variation (division alwys orthogonal to the max. lengh)
const float DIVISION_ANGLE_LONGEST_AXIS = 0.7;
const float DIVISION_ANGLE_EXTERNAL = 0.0; // 
const float DIVISION_ANGLE_EXTERNAL_DEGREES = 0;
const float CELL_CYCLE_LIMIT = 10;
const float TIME_DECREASE_EXPONENT = 0.5; //Only active if TIME_CONTROLS_SIZE
const float XCOORD_DECREASE_EXPONENT = 0.5; //Only active if XCOORD_CONTROLS_SIZE

const float ENERGY_TERM1 = 1;
const float ENERGY_TERM2 = 1;
const float ENERGY_TERM3 = 1;

const float LENGTH_ROTATED_EDGE = 0.5*T1_TRANSITION_CRITICAL_DISTANCE*1.2; //after a t1 transition, rotated edge length is multiplied by twice this constant


// Implementation constants
const unsigned short CELLS_PER_VERTEX = 3;
const unsigned short MAX_SIDES_PER_CELL = 30; //Maybe pass as argument instead of constant?
const unsigned short VERTEX_PER_EDGE = 2;
const int EMPTY_CONNECTION = -999;  //Value to initialize arrays (cells in vertex, etc)
const int RANDOM_SEED = 1234;
const double NUMERIC_THRESHOLD = 1e-16;
const double EXP_FACTOR = 1 - exp(-1); //Used when time or x position changes eq. size, to normalize the function 1 - exp(-(x**exp)) so its max value is 1

const std::string VERTEX_FILE_EXTENSION = ".points";
const std::string CELLS_FILE_EXTENSION = ".cells";
const std::string CELLTAB_FILE_EXTENSION = ".celltab";
const std::string EDGE_FILE_EXTENSION = ".edges";
const std::string SPRING_FILE_EXTENSION = ".spr";
const std::string PARAMS_FILE_EXTENSION = ".vp";

const bool T1_ACTIVE = true;
const bool T1_BORDER_INWARDS_ACTIVE = false;
const bool T1_BORDER_OUTWARDS_ACTIVE = false;
const bool DIVISION_ACTIVE = true;
const bool T2_ACTIVE = true;
const bool JOIN_EDGES_ACTIVE = true;
const bool CONTROL_CELLS_2SIDES = true; 
const bool CHECK_EDGES_CROSS_AFTER_MOVE = false;
const bool AUTONOMOUS_CELL_CYCLE = true;
const bool KEEP_AREA_AFTER_DIVISION = false;
const bool CELL_CYCLE_CONTROLS_SIZE = true;
const bool TIME_CONTROLS_SIZE = false;
const bool XCOORD_CONTROLS_SIZE = false;
const bool STATIC_PRESENT = true;
const int MOVE_TRIALS = 100;  //Times it tries to move a vertex before it quits because always makes edges to cross
const bool UPDATE_EDGE_TENSION_EVERY_MOVE = true; //Update edge tension if tension is dependent on edge angle (if true, do it in every iteration, if false do it only in transitions and start)
const bool REPORT_T1 = false;
const bool REPORT_DIV = false;
const bool REPORT_OUT = true;//SET THIS TO FALSE IF EVOLUTION IS GOING TO BE USED

const std::string VERTEX_HEADER = "ind\tx\ty\tenergy\tmovable\tspring\tcells\tedges\tneighbour_vertices\n";
const std::string CELL_HEADER = "ind\ttype\tarea\tpreferred_area\tperimeter\tperim_contract\tcentroid_x\tcentroid_y\tangle_longest\tangle_signal\tangle_random\tdegrees_signal\tmax_area\tcell_cycle_state\tcell_cycle_limit\tcan_divide\tnum_vertices\tvertices\tedges\tnum_divisions\n";
const std::string EDGE_HEADER = "ind\ttype\tlength\ttension\tvertices\tcells\n";

//Enum class to define types of cells 
enum class CellType{blade = 0, hinge = 1, vein_blade = 2, vein_hinge = 3};

//Enum class to define types of edges
enum class EdgeType{blade = 0, hinge = 1, vein_hinge = 2, vein_blade = 3, tissue_boundary = 4, spring = 5, vein = 6};

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
//Add: preferred angle with respect to origin/longest axis and angle noise
//     cell cycle stage
//     etc
struct Cell{
	int ind;
	CellType type;
	double area;
	double preferred_area;
	double perimeter;
	double perimeter_contractility;
	float division_angle_random_noise; 
	float division_angle_longest;
	float division_angle_external; 
	float division_angle_external_degrees;
	bool vary_line_tension;
	float edge_angle_prop_external;
	float edge_angle_prop_uniform;
	float edge_angle_prop_maxangle;
	float edge_angle_prop_random;
	float edge_tension_external;
	float edge_maxangle;
	float edge_spatialmax_tension;
	float edge_spatialmin_tension;
	float max_area;
	float cell_cycle_state;
	float cell_cycle_limit;
	bool can_divide;
	double centroid_x, centroid_y;
	int edges[MAX_SIDES_PER_CELL];
	int vertices[MAX_SIDES_PER_CELL];
	int num_vertices; 
	bool dead;
	int num_divisions;
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
        bool vertical;
	double intercept;
	double slope;
	double x1;
	double x2;
	double y1;
	double y2;
	int v1;
	int v2;
};
struct DivisionRecord{
	int parent;
	int offspring;
};

//Vector types for each structure
typedef std::vector<Vertex> vertex_v;
typedef std::vector<Cell> cell_v;
typedef std::vector<Edge> edge_v;
typedef std::queue<Rearrangement> rearrangement_q;
typedef std::queue<DivisionRecord> divisionrecord_q;

typedef std::map<CellType, double> cell_type_param;
typedef std::map<int, double> spring_type_param;


class Tissue{
        friend class basicGRN;
	public:
		Tissue();
		Tissue(std::string starting_tissue_file, int max_accepted_movements=1000000,  int write_every_N_moves=1000, std::string simulname="");
		Tissue(std::string starting_tissue_file, std::string params_file, int max_accepted_movements=1000000,  int write_every_N_moves=1000, std::string simulname=""); //not implemented

		void initialize_params(std::string params_file="");
		double read_real_par(std::vector<std::string>::iterator& it);
		cell_type_param read_celltype_par(std::vector<std::string>::iterator& it, std::string::size_type sz);
                spring_type_param read_springtype_par(std::vector<std::string>::iterator& it, std::string::size_type sz);
		void set_default_simulation_params();
		void setHingeMinAndMaxPositions();

		void simulate(std::default_random_engine& generator, std::uniform_real_distribution<double>& unif);
		double calculateCellArea(const Cell& c);
		double calculateCellPerimeter(const Cell& c);
		void calculateCellCentroid(Cell& c);
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
		void addAcceptedMovements(int add);
		void setAcceptedMovements(int mv);
		void setStepMode(bool mode, int steps);

		void writeCellsFile(std::string fname);
		void writeCellDataTable(std::string fname);
		void writeEdgeDataTable(std::string fname);
		void writePointsFile(std::string fname);
		void writeSpringsFile(std::string fname);
		void writeAllData(std::string fname);
		void produceOutputs(std::string add_to_name="moved");
		std::string getStats();
		int getCounterT1();
		void emptyDivisions();

		void addSpringsAutomatically();
		void setStaticVertex(int v);
		bool AddSpringToVertex(int v, float minx, float maxx);
		void readNewParameters(std::string filename);

		std::vector<int> getNeighbourCells(int cell);
		
		friend std::ostream& operator<<(std::ostream& out, const Tissue& t); //Don't worry about these 4 lines, only important for output
		friend std::ostream& operator<<(std::ostream& out, const Vertex& v);
		friend std::ostream& operator<<(std::ostream& out, const Cell& c);
		friend std::ostream& operator<<(std::ostream& out, const Edge& e);
		
	private:
		std::string simname;
		int num_cells, num_vertices, num_edges, num_springs;
		double hinge_min_xpos, hinge_max_xpos, hinge_min_ypos, hinge_max_ypos;
		bool step_mode; //If the simulation is controlled by an external source, set this to true in order to avoid excessive printing of outputs
		//parameters
		float min_range_vertex_movement; 
		float max_range_vertex_movement;
		float temperature_positive_energy;
		float temperature_negative_energy; 
		bool temperature_means_proportion_of_acceptance;
		cell_type_param line_tension;
		cell_type_param line_tension_tissue_boundary;
		float energy_term1, energy_term2, energy_term3;
		//float spring_constant;
                spring_type_param spring_type_constants;
                spring_type_param spring_type_min_positions;
		float add_static_to_hinge;
		cell_type_param perimeter_contract;

		float t1_transition_critical_distance; 
		float t2_transition_critical_area;
		//float max_cell_area;	
		float max_edge_length; 
		cell_type_param preferred_area_initial, preferred_area_final, max_cell_area;

		cell_type_param division_angle_longest_axis, division_angle_random_noise, division_angle_external, division_angle_external_degrees;

		cell_type_param cell_cycle_limit, xcoord_size_control_factor;
		bool autonomous_cell_cycle, start_cell_cycle_at_random, cell_cycle_controls_size, time_controls_size, xcoord_controls_size, keep_area_after_division;
		float time_decrease_exponent, xcoord_decrease_exponent;

		//Parameters related to edge tension modification
		cell_type_param vary_line_tension;
		bool vary_edge_tension_with_time;
		float vary_edge_tension_time_exponent;
		cell_type_param edge_angle_prop_external;
		cell_type_param edge_angle_prop_uniform;
		cell_type_param edge_angle_prop_maxangle;
		cell_type_param edge_angle_prop_random;
		cell_type_param edge_tension_external; //initial values
		cell_type_param edge_maxangle;
		cell_type_param edge_spatialmax_tension;
		cell_type_param edge_spatialmin_tension;
		cell_type_param edge_temporal_angle_efect_max;
		cell_type_param edge_temporal_angle_efect_min;

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
		divisionrecord_q past_divisions;
		//counters
		int counter_move_trials, counter_moves_accepted, counter_favorable_accepted, counter_favorable_rejected, counter_unfav_accepted, counter_unfav_rejected, counter_edges_removed, counter_t1, counter_t1_abortions, counter_divisions, counter_t2, counter_t1_outwards, counter_t1_inwards;
		int max_accepted_movements, write_every_N_moves;

		//Methods used by constructors
		void initialize_vertices(std::ifstream& inp);
		void initialize_cells(std::ifstream& inp);  
		void initialize_edges();
		void initialize_springs(std::ifstream&);
		void addEdge(int v1, int v2, int cell);		//Creates single edge
		void set_default_params();			//Sets model parameters for each vertex, cell and edge, from constants defined in this file
		void setEdgeType(int e);
		void setEdgeTension(int e);
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
		bool getDivisionPoints(const int cell, double &x1, double &x2, double &y1, double &y2, int &e1, int &e2);
		bool findEdgesToCut(const int cell, double x1, double x2, double y1, double y2, int &e1, int &e2, StraightLine &l1, StraightLine &l2, StraightLine &l3, StraightLine &l4);
		void splitEdgeWithVertex(int e, int cell, int  v);

		//methods used by T2 and join_edges
		void removeConnectionCell(int elm, int* elements, int length);
		void make_remove_size2cell(int cell);
		void killCellWith2Vertices(int cell);

		//Other
		double expAdvance(double x, float exponent);
		void advanceCellCycle(int vertex_moved);
		void advanceSizeWithTime(int vertex_moved);
		void advanceSizeWithXcoordAndTime(int vertex_moved);
		void advanceSizeWithXcoord(int vertex_moved);
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


