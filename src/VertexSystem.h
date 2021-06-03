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
const float MAX_RANGE_VERTEX_MOVEMENT = 0.02; //Monte Carlo
const float DEFAULT_H = 0.01; //Euler, etc
const float TEMPERATURE_POSITIVE_ENERGY = 0.005;
const float TEMPERATURE_NEGATIVE_ENERGY = 0.1; //Acceptance probability is basically 1
const bool TEMPERATURE_MEANS_PROPORTION_OF_ACCEPTANCE = true;
//const float PROB_ACCEPT_NON_FAVORABLE = 0.05;
const float DEFAULT_K = 0.001;
const float LINE_TENSION_BLADE = -0.02;
const float LINE_TENSION_HINGE = -0.01;
const float LINE_TENSION_VEIN_HINGE = 7.0;
const float LINE_TENSION_VEIN_BLADE = 1.0;
const float LINE_TENSION_TISSUE_BOUNDARY = 4.0;//4.0;
const float LINE_TENSION_STRINGEDGE = 10.0;
const float HINGE_BLADE_INTERFACE_TENSION = 1.0;
const float LINE_TENSION_INTERSTATIC = -0.5;
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
const bool USE_TERM4 = true;

const float LENGTH_ROTATED_EDGE = 0.5*T1_TRANSITION_CRITICAL_DISTANCE*1.2; //after a t1 transition, rotated edge length is multiplied by twice this constant


// Implementation constants
const unsigned short CELLS_PER_VERTEX = 3;
const unsigned short MAX_SIDES_PER_CELL = 100; //Maybe pass as argument instead of constant?
const unsigned short VERTEX_PER_EDGE = 2;
const int EMPTY_CONNECTION = -999;  //Value to initialize arrays (cells in vertex, etc)
const int RANDOM_SEED = 1234;
const double NUMERIC_THRESHOLD = 1e-5;
const double EXP_FACTOR = 1 - exp(-1); //Used when time or x position changes eq. size, to normalize the function 1 - exp(-(x**exp)) so its max value is 1
const int INTEGR_MONTECARLO = 0;
const int INTEGR_EULER = 1;
const int INTEGR_MIXTURE = 2;
const int INTEGR_RUNGEKUTTA4 = 3;
const int RESTORE_VEIN_ITERS = 1; //iterations in algorithm to re-fill veins after many movements 

const std::string VERTEX_FILE_EXTENSION = ".points";
const std::string CELLS_FILE_EXTENSION = ".cells";
const std::string CELLTAB_FILE_EXTENSION = ".celltab";
const std::string EDGE_FILE_EXTENSION = ".edges";
const std::string POINTTAB_FILE_EXTENSION = ".ptab";
const std::string SPRING_FILE_EXTENSION = ".spr";
const std::string STRINGEDGE_FILE_EXTENSION = ".stre";
const std::string SPRINGTAB_FILE_EXTENSION = ".sprtab";
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
const bool CELL_CUTICLE_CAN_TRANSITION = false; //Be careful! true, after transitions tension is set according to cell param, not cuticle param
const bool BORDER_CAN_TRANSITION = false;
const bool BLADE_HINGE_INTERFACE_CAN_TRANSITION = false;
const bool VEIN_TO_BLADE_CAN_TRANSITION = true;
const bool VEIN_TO_HINGE_CAN_TRANSITION = true;
const bool VEIN_CAN_TRANSITION = true;

const bool STATIC_PRESENT = true;
const int MOVE_TRIALS = 100;  //Times it tries to move a vertex before it quits because always makes edges to cross
const bool UPDATE_EDGE_TENSION_EVERY_MOVE = true; //Update edge tension if tension is dependent on edge angle (if true, do it in every iteration, if false do it only in transitions and start)
const bool REPORT_T1 = false;
const bool REPORT_DIV = false;
const int REPORT_OUT = 0;//SET THIS TO 0 IF EVOLUTION IS GOING TO BE USED OR IN CASE OF LONG SIMULATIONS
const bool WRITE_DATA_TABLES = true;
const bool STRING_EQUILIBRIUM_DISTANCE = false;
const bool ONLY_ONE_CONTACT = true; //Only 1 contact edge is allowed between 2 cells; T1 is not performed when it will compromise this condition
const bool RECALCULATE_CENTROIDS_FOR_PRINTING = true;
const int USE_HINGE_BLADE_FRONTIER = 1;
const int USE_MAX_HINGE_POSITION = 2;
const int USE_PROPORTION_OF_WING = 3;
const int USE_BARE_EXPONENTIAL_GRADIENT = 4;
const int GRADIENT_FROM_CENTER_ALL = 5;
const int GRADIENT_FROM_CENTER_HINGE = 6;
const float DEFAULT_PROPORTION_FOR_GRADIENT = 0.7;
const int REFERENCE_FOR_GRADIENT = USE_PROPORTION_OF_WING;
const bool NORMALIZE_EDGE_TENSION = false;

const std::string VERTEX_HEADER = "ind\tx\ty\tenergy\tmovable\ttype\tspring\tmoves_accepted\tmoves_rejected\tcells\tedges\tneighbour_vertices\n";
const std::string CELL_HEADER = "ind\ttype\tarea\tpreferred_area\tK\tperimeter\tperim_contract\t" +
	std::string("centroid_x\tcentroid_y\tbase_eq_area\tbase_eq_perimcontr\tangle_longest\tangle_signal\t")+
	std::string("angle_random\tdegrees_signal\tmax_area\tcell_cycle_state\tcell_cycle_limit\tcan_divide\tnum_vertices\t")+
	std::string("vertices\tedges\tnum_divisions\t")+
	std::string("vary_line_tension\tedge_angle_prop_external\tedge_angle_prop_uniform\tedge_angle_prop_maxangl\t")+
	std::string("edge_angle_prop_random\tedge_tension_external\t")+
	std::string("edge_maxangle\tedge_spatialmax_tension\tedge_spatialmin_tension\tfather\tinitial_xpos\n");
const std::string EDGE_HEADER = "ind\ttype\tlength\ttension\tbase_tension\tvertices\tcells\n";
const std::string SPR_HEADER = "ind\ttype\ttension\tlength\tcompartment\tstatic_vertex\tmovable_vertex\tx_static\ty_static\tx_movable\ty_movable\n";

const std::string T1_HEADER = ">T1_name\tmoves_accepted\ttype\tv1\tv2\tx1_before\tx1_after\tx2_before\tx2_after\ty1_before\ty1_after\ty2_before\ty2_after\tangle";
const std::string T2_HEADER = ">T2_name\tcell\tvert_survivor\tx\ty\ttype";
const std::string DIV_HEADER = ">DIV_name\tind1\tind2\tcentroid_1x\tcentroid_1y\tcentroid_2x\tcentroid_2y\ttype\tcounter_moves_accepted";

//Enum class to define types of cells 
enum class CellType{blade = 0, hinge = 1, vein_blade = 2, vein_hinge = 3, cuticle=4};

//Enum class to define types of edges
enum class EdgeType{blade = 0, hinge = 1, vein_hinge = 2, vein_blade = 3, tissue_boundary = 4, spring = 5, vein = 6, stringedge = 7, cellular_cuticle = 8};

//Enum class to define types of vertices
enum class VertexType{tissue=0, cuticle=1, spring_only=2, cuticle_2layers=3}; //1 accounts for cell cuticle and for string cuticle 

//Enum class to define types of transitions
enum class RearrangementType{t1 = 0, t2 = 1, divide_cell = 2, divide_edge = 3, join_limit_edges = 4, t1_at_border_outwards = 5, t1_at_border_inwards = 6, rotate_inwards = 7, random_t1 = 8};

//Enum class to define cuticle type
enum class CuticleType{only_springs = 0, strings_1layer = 1, strings_2layers = 2, cells = 3};

// Structure for Vertex
const int MAX_NODES_PER_CUTICLE_VERTEX = 8;
struct Vertex{
	int ind;
	double x;
	double y;
	double energy;
	int cells[CELLS_PER_VERTEX];
	union{	//Unions just save some memory; edges2 and neigobours2 hold data for cuticle vertices when cuticle is made of 2 layers of strings
		int edges[CELLS_PER_VERTEX];
		int edges2[MAX_NODES_PER_CUTICLE_VERTEX];
	};//In this case it would be easier to declare just a vector of size 5, but I did not want to change other parts of the code
	union{
		int neighbour_vertices[CELLS_PER_VERTEX];
		int neighbour_vertices2[MAX_NODES_PER_CUTICLE_VERTEX];
	};
	int moves_accepted;
	int moves_rejected;
	bool movable; //>0  in param_file
	bool movable_x; //== 2 in param_file, == 1 movable both
	bool movable_y; // == 3  in param_file, == 1 movable both
	int spring;
	bool dead;
	double x_drag;
	double y_drag;
	VertexType type;
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
	double base_preferred_area;
	double K;
	double perimeter;
	double base_perimeter;
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
	int father;
	float x_at_start;
};
// Structure for Edge
struct Edge{
	int ind;
	EdgeType type;
	double tension;
	float base_tension;
	double length;
	float optimal_length;
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

//typedef std::map<CellType, double> cell_type_param;
//typedef std::map<int, double> spring_type_param;

const int NUM_CELL_TYPES = 5; //SHOULD BE A MAXIMUM BECAUSE STRUCT TYPES ARE DECLARED USING THIS NUMBER
const int NUM_SPRING_TYPES = 4;
struct cell_type_param{
	double val[NUM_CELL_TYPES]; //If you use a cell type number > NUM_CELL_TYPES, an error will occur, even when inside the program the actual cell number is known
};
struct spring_type_param{
	double val[NUM_SPRING_TYPES];
};

struct bufferRejectMovement{
	double edge_lengths[MAX_NODES_PER_CUTICLE_VERTEX];
	double edge_tensions[MAX_NODES_PER_CUTICLE_VERTEX];
	//double edge_tensions2[MAX_NODES_PER_CUTICLE_VERTEX];
	//double edge_lengths2[MAX_NODES_PER_CUTICLE_VERTEX];
	double cell_areas[CELLS_PER_VERTEX];
	double cell_perimeters[CELLS_PER_VERTEX];
	double spring_length;
	double x;
	double y;
	double energy;

};
struct pointDerivative{
	double x;
	double y;
};
typedef std::vector<pointDerivative> pointDerivative_v;

class Tissue{
        friend class basicGRN;
	public:
		Tissue();
		Tissue(std::string starting_tissue_file, int max_accepted_movements=1000000,  int write_every_N_moves=1000, std::string simulname="");
		Tissue(std::string starting_tissue_file, std::string params_file, int max_accepted_movements=1000000,  int write_every_N_moves=1000, std::string simulname=""); //not implemented

		void initialize_params(std::string params_file="");
		double read_real_par(std::vector<std::string>::iterator& it);
		double read_longint_par(std::vector<std::string>::iterator& it);
		cell_type_param read_celltype_par(std::vector<std::string>::iterator& it, std::string::size_type sz);
        spring_type_param read_springtype_par(std::vector<std::string>::iterator& it, std::string::size_type sz);
		void set_default_simulation_params();
		void setHingeMinAndMaxPositions();
		void setMinAndMaxPositions();

		void simulate(std::default_random_engine& generator, std::uniform_real_distribution<double>& unif);
		void simulateEuler();
		void simulateMonteCarlo();	
		double calculateCellArea(const Cell& c);
		double calculateCellPerimeter(const Cell& c);
		void calculateCellCentroid(Cell& c);
		double distance(int v1, int v2);
		double calculateEnergy(Vertex& v);
		double calculateEnergyCuticle2(Vertex& v);
		double calculateEnergy2(Vertex& v);
		double calculateEnergy_term4(Vertex& v);
		double calculateEnergy_term4_cuticle(Vertex& v);
		void derivativeVertexPos(const Vertex &v, pointDerivative & pd);
		void moveVertex(Vertex& v, float x, float y);
		void moveVertexBack(Vertex& v);
		void moveVertexCuticle(Vertex& v, float x, float y);
		void moveVertexBackCuticle(Vertex& v);
		bool tryMoveVertex();
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
		void emptyDivisions();
		//Methods for printing results
		void writeCellsFile(std::string fname);
		void writeCellDataTable(std::string fname);
		void writeEdgeDataTable(std::string fname);
		void writeSpringDataTable(std::string fname);
		void writePointsDataTable(std::string fname);
		void writePointsFile(std::string fname);
		void writeSpringsFile(std::string fname);
		void writeAllData(std::string fname);
		void produceOutputs(std::string add_to_name="moved");
		std::string getStats();
		void printLine(std::string name, int ind1, int ind2, double centroid_x, double centroid_y, int type);
		void printLine_t1(std::string name, int ind1, int ind2, 
			double v1_x0, double v2_x0, double v1_y0, double v2_y0, 
			double v1_x1, double v2_x1, double v1_y1, double v2_y1, int type);
		void printLine_t2(std::string name, int ind1, int ind2, double centroid_x, double centroid_y, int type);
		void printLine_div(std::string name, int ind1, int ind2, double centroid_1x, double centroid_1y, 
		double centroid_2x, double centroid_2y, int type);
		int getCounterT1();
		void printCelltypeParam(cell_type_param par, std::string name);
	 	//Methods to transition from expansion to hinge contraction (adding springs etc)
		void addSpringsAutomatically();
		void setStaticVertex(int v);
		bool AddSpringToVertex(int v, float minx, float maxx);
		void readNewParameters(std::string filename);
		void makeVeinsThinner(int iters=1);
		void restoreHinge();
		void restoreVeins();
		void restoreShape();
		//Methods to set spring tensions according to gradients
		void setSpringTension();
		void setSpringTension_mode1(); // A-P compartments (uses tension determined by spring type)
		void setSpringTension_mode2(); // P-D gradient
		void setSpringTension_mode3(); //P-D gradient multiplied by a factor in each compartment
		void setSpringTension_mode4(); // A different (independent) gradient in each compartment
 
		std::vector<int> getNeighbourCells(int cell);
		
		friend std::ostream& operator<<(std::ostream& out, const Tissue& t); //Don't worry about these 4 lines, only important for output
		friend std::ostream& operator<<(std::ostream& out, const Vertex& v);
		friend std::ostream& operator<<(std::ostream& out, const Cell& c);
		friend std::ostream& operator<<(std::ostream& out, const Edge& e);
		
	private:
		std::string simname;
		int num_cell_types;
		int num_cells, num_vertices, num_edges, num_springs;
		int integration_mode;
		bool step_mode; //If the simulation is controlled by an external source, set this to true in order to avoid excessive printing of outputs
		//parameters
		bool t1_active, t1_inwards_active, t1_outwards_active, division_active, t2_active, join_edges_active, control_cells_2sides, check_if_edges_cross_opt;
		bufferRejectMovement bufferMovement;
		float min_range_vertex_movement; 
		float max_range_vertex_movement; //In case of Monte Carlo integration
		float h; //In case of Euler/Runge-Kutta integration
		float temperature_positive_energy;
		float temperature_negative_energy; 
		bool temperature_means_proportion_of_acceptance;
		cell_type_param K, K_final;
		bool K_gradient_x, K_gradient_y;
		float K_grad_exponent;
		cell_type_param line_tension;
		cell_type_param line_tension_tissue_boundary;
		float energy_term1, energy_term2, energy_term3;
		cell_type_param energy_term4_posterior, energy_term4_anterior;
		bool use_term4;
		//float spring_constant;
        spring_type_param spring_type_constants; //Used if spring_tension_mode ==0
        spring_type_param spring_type_min_positions; //Used when going from expansion to hinge contraction
		float add_static_to_hinge; //Used when going from expansion to hinge contraction, to set static vertices in anterior hinge
		int spring_tension_mode; //0: with tension for each type, 1: A-P compartments, 2: PD gradient, 3: AP compartments and PD gradient, 4:Different gradients in A or P
		float posterior_comparment_region; //used if spring_tension_mode is 1 or 3
		float spring_posterior_compartment_factor; //used if spring_tension_mode is 1 or 3
		float spring_gradient_min_tension; //used if spring_tension_mode is 2 or 3
		float spring_gradient_max_tension; //used if spring_tension_mode is 2 or 3
		float spring_gradient_exponent; //used if spring_tension_mode is 2 or 3
		float spring_gradient_min_tension_P; //used if spring_tension_mode is 4
		float spring_gradient_max_tension_P; //used if spring_tension_mode is 4
		float spring_gradient_exponent_P; //used if spring_tension_mode is 4
		CuticleType cuticle_type;
		float line_tension_interstatic;
		float AP_compartment_limit;	
		int mode_to_order_springs_PD; //used if spring_tension_mode is 2 or 3
		float line_tension_stringedge, tension_stringedge_posterior_prop, string_edge_tension_min, string_edge_tension_exponent;
		float string_distal_transition_tension, string_distal_transition_prop, hinge_string_tension;
		bool string_anterior_gradient, string_posterior_gradient, set_hinge_string_tension, include_hinge_in_spring_gradient;
		
		cell_type_param perimeter_contract, perimeter_contract_final;
		bool time_controls_perim, xcoord_controls_perim, ycoord_controls_perim; //Exponent will be the same as for area
		bool string_equilibrium_distance;
		float t1_transition_critical_distance; 
		float t2_transition_critical_area;
		float active_t1_prob,min_angle_for_active_t1, max_angle_for_active_t1, minsin2rant1, maxsin2rant1; //To make T1s at random
		//float max_cell_area;	
		float max_edge_length; 
		cell_type_param preferred_area_initial, preferred_area_initial_gradient, preferred_area_final, max_cell_area;
		cell_type_param division_angle_longest_axis, division_angle_random_noise, division_angle_external, division_angle_external_degrees;
		cell_type_param cell_cycle_limit, xcoord_size_control_factor;
		bool autonomous_cell_cycle, start_cell_cycle_at_random, cell_cycle_controls_size;
		bool time_controls_size, xcoord_controls_size, ycoord_controls_size, use_blade_area_for_coord_gradient, gradient_with_same_final_area;
		bool keep_area_after_division;
		int reference_for_gradient;
		float time_decrease_exponent, xcoord_decrease_exponent, wing_proportion_in_gradient;
		int random_seed;
		//float difference_flow_rate;

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
		float hinge_blade_interface_tension;
		float length_rotated_edge; 
		double hinge_min_xpos, hinge_max_xpos, hinge_min_ypos, hinge_max_ypos, max_xpos, max_ypos, min_xpos, min_ypos;
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
		long unsigned int counter_move_trials, counter_moves_accepted, counter_favorable_accepted, counter_favorable_rejected, counter_unfav_accepted, counter_unfav_rejected, counter_edges_removed, counter_t1, counter_t1_abortions, counter_divisions, counter_t2, counter_t1_outwards, counter_t1_inwards;
		long unsigned int max_accepted_movements, write_every_N_moves, upper_bound_movements;
		int written_files;
		std::default_random_engine generator;
		std::uniform_real_distribution<double> unif;
		//drag / momentum parameters:
		float momentum_term_cuticle, momentum_term_tissue, momentum_ponderate_current;
		//Methods used by constructors
		void initialize_vertices(std::ifstream& inp);
		void initialize_cells(std::ifstream& inp);  
		void initialize_edges();
		void initialize_springs(std::ifstream&);
		void initialize_stringedges(std::ifstream&);
		void addEdge(int v1, int v2, int cell);		//Creates single edge
		void set_default_params();			//Sets model parameters for each vertex, cell and edge, from constants defined in this file
		void setEdgeType(int e);
		void setEdgeTension(int e);
		void setStringTension(int e); //For outer strings (cuticle)
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
		void calculateBasePrefAreaAndPerim(Cell& cell);
		void advanceCellCycle(int vertex_moved);
		void advanceSizeWithTime(int vertex_moved);
		void advanceSizeAndPerimWithTime_new(int vertex_moved);
		void advancePerimWithTime(int vertex_moved);
		void advanceSizeWithCoordsAndTime(int vertex_moved);
		void advancePerimWithCoordsAndTime(int vertex_moved);
		void advanceSizeWithCoordsAndTime_v2(int vertex_moved);
		void advanceSizeWithCoords(int vertex_moved);
		void advancePerimWithCoords(int vertex_moved);
		void advanceKWithCoords(int vertex_moved);
		void setKWithCoords(int cell);
		std::vector<float> getSpringGradientFactor_mode1(std::vector<int> &vert_indices);
		std::vector<float> getSpringGradientFactor_mode2(std::vector<int> &vert_indices);
		std::vector<int> getHingeBladeFrontier();
		double getXgradFromFrontier(float xpos, float ypos, std::vector<int> & frontier);
		void addRandomTransition();
		bool checkCellConsistency(int cell);
		//double calculateTerm4Energy(Vertex &v, double old_x, double old_y);
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


