

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <random>
#include <limits>
#include "VertexSystem.h"
//#include <chrono>

using namespace std;

//Default constructor
Tissue::Tissue() : num_cells(0), num_vertices(0), num_edges(0), counter_move_trials(0), written_files(0), counter_moves_accepted(0), counter_favorable_accepted(0), counter_favorable_rejected(0), counter_unfav_accepted(0), counter_unfav_rejected(0), counter_t1(0), counter_t1_abortions(0), counter_edges_removed(0), counter_divisions(0), counter_t2(0), counter_t1_outwards(0), counter_t1_inwards(0)
{
	cout << "Default constructor\n";
	simname = "";
	max_accepted_movements = 0;
	upper_bound_movements = 0;
	write_every_N_moves = 0;
	step_mode = false;

	setMinAndMaxPositions();
	setHingeMinAndMaxPositions();
	set_default_simulation_params();
}

//Constructor that reads vertex positions and cells from two different files, and initializes all variables using constants defined in VertexSystem.h
Tissue::Tissue(std::string starting_tissue_file, int max_accepted_movements, int write_every_N_moves, string simulname) : Tissue()
{
	cout << "Initializing Tissue object without params_file...\n";
	this->simname = simulname == "" ? starting_tissue_file : simulname;
	this->max_accepted_movements = max_accepted_movements;
	this->upper_bound_movements = max_accepted_movements;
	this->write_every_N_moves = write_every_N_moves;
	step_mode = false;

	//Read file of vertices (indicates coordinates for each vertex)
	string vertexfile = starting_tissue_file + VERTEX_FILE_EXTENSION;
	std::ifstream fin_vertex;
	if (REPORT_OUT > 0)
		cout << "reading .points file...\n";
	fin_vertex.open(vertexfile);
	initialize_vertices(fin_vertex);
	fin_vertex.close();
	if (REPORT_OUT > 0)
		cout << ".points file read...\n";

	//Read file of cells (indicates vertices for each cell)
	try
	{
		string cellfile = starting_tissue_file + CELLS_FILE_EXTENSION;
		if (REPORT_OUT > 0)
			cout << "reading .cells file...\n";
		ifstream fin_cells(cellfile);
		initialize_cells(fin_cells);
		fin_cells.close();
		if (REPORT_OUT > 0)
			cout << ".cells file read...\n";
	}
	catch (const char *msg)
	{
		exit(1);
	}
	//Initialize edges from vertices and cells

	initialize_edges();

	try
	{
		string springsfile = starting_tissue_file + SPRING_FILE_EXTENSION;
		if (REPORT_OUT > 0)
			cout << "reading .spr file...\n";
		ifstream fin_springs(springsfile);
		if (fin_springs.good())
		{
			initialize_springs(fin_springs);
			if (REPORT_OUT > 0)
				cout << ".spr file read...\n";
		}
		else
		{
			num_springs = 0;
		}
		fin_springs.close();
	}
	catch (const char *msg)
	{
		num_springs = 0;
		if (REPORT_OUT > 0)
			cout << msg << endl;
		if (REPORT_OUT > 0)
			cout << "No spring file\n";
	}
	try
	{
		string strefile = starting_tissue_file + STRINGEDGE_FILE_EXTENSION;
		ifstream fin_stre(strefile);
		if (fin_stre.good())
		{
			if (REPORT_OUT > 0)
				cout << "reading .stre file...\n";
			initialize_stringedges(fin_stre);
			if (REPORT_OUT > 0)
				cout << ".stre file read...\n";
		}
		else
		{
			cout << "No stre file\n";
		}
		fin_stre.close();
	}
	catch (const char *msg)
	{
		if (REPORT_OUT > 0){
			cout << msg << endl;
			cout << "No stre file\n";
		}
	}
	//No file with parameters, therefore use constants defined in VertexSystem.h.
	//Also initializes cell area and vertex energy
	setMinAndMaxPositions();
	setHingeMinAndMaxPositions();
	set_default_params();
	
	if (REPORT_OUT > 0)
		cout << "parameters set to cells\n";
		cout << T1_HEADER << endl << DIV_HEADER << endl << T2_HEADER << endl;
}

Tissue::Tissue(std::string starting_tissue_file, std::string params_file, int max_accepted_movements, int write_every_N_moves, string simulname) : num_cells(0), num_vertices(0), num_edges(0), counter_move_trials(0), written_files(0), counter_moves_accepted(0), counter_favorable_accepted(0), counter_favorable_rejected(0), counter_unfav_accepted(0), counter_unfav_rejected(0), counter_t1(0), counter_t1_abortions(0), counter_edges_removed(0), counter_divisions(0), counter_t2(0), counter_t1_outwards(0), counter_t1_inwards(0)
{
	cout << "Initializing Tissue object with params_file...\n";
	std::flush(cout);
	this->simname = simulname == "" ? starting_tissue_file : simulname;
	this->max_accepted_movements = 0; //max_accepted_movements; //Will be overwritten by Parameter File unless value in file is 0!! (only conserved for compatibility with elder scripts)
	this->upper_bound_movements = max_accepted_movements; //Will be overwritten
	this->write_every_N_moves = 0; //write_every_N_moves; //Will be overwritten by Parameter File unless value in file is 0!!
	step_mode = false;

	try
	{
		if (REPORT_OUT > 0){
			cout << "reading param file...\n";
			std::flush(cout);
		}
		initialize_params(params_file);
		if (REPORT_OUT > 0){
			cout << "param file read...\n";
			std::flush(cout);
		}
	}
	catch (const char *msg)
	{
		cout << "Error in param file" << endl;
		cout << msg << endl;
		exit(1);
	}
	//Read file of vertices (indicates coordinates for each vertex)
	try
	{
		string vertexfile = starting_tissue_file + VERTEX_FILE_EXTENSION;
		std::ifstream fin_vertex(vertexfile);
		if (REPORT_OUT > 0){
			cout << "reading .points file...\n";
			std::flush(cout);
		}
		//fin_vertex.open(vertexfile);
		if(! fin_vertex.good()) throw ".points file not present";
		initialize_vertices(fin_vertex);
		fin_vertex.close();
		if (REPORT_OUT > 0){
			cout << ".points file read...\n";
			std::flush(cout);
		}
	}
	catch (const char *msg)
	{
		cout << "Error in vertices file" << endl;
		cout << msg << endl;
		exit(1);
	}
	//Read file of cells (indicates vertices for each cell)
	try
	{
		string cellfile = starting_tissue_file + CELLS_FILE_EXTENSION;
		if (REPORT_OUT > 0){
			cout << "reading .cells file...\n";
			std::flush(cout);
		}
		ifstream fin_cells(cellfile);
		if(! fin_cells.good()) throw ".cells file not present";
		initialize_cells(fin_cells);
		fin_cells.close();
		if (REPORT_OUT > 0){
			cout << ".cells file read...\n";
			std::flush(cout);
		}
	}
	catch (const char *msg)
	{
		cout << "Error in cell file" << endl;
		cout << msg << endl;
		exit(1);
	}
	//Initialize edges from vertices and cells
	if (REPORT_OUT > 0)
		cout << "Initializing edges...\n";
	initialize_edges();
	try
	{
		string springsfile = starting_tissue_file + SPRING_FILE_EXTENSION;
		ifstream fin_springs(springsfile);
		if (fin_springs.good())
		{
			if (REPORT_OUT > 0){
				cout << "reading .spr file...\n";
				std::flush(cout);
			}
			initialize_springs(fin_springs);
			if (REPORT_OUT > 0){
				cout << ".spr file read...\n";
				std::flush(cout);
			}
		}
		else
		{
			num_springs = 0;
			cout << "No spring file\n";
			std::flush(cout);
		}
		fin_springs.close();
	}
	catch (const char *msg)
	{
		num_springs = 0;
		if (REPORT_OUT > 0){
			cout << msg << endl;
			cout << "No spring file\n";
			std::flush(cout);
		}
	}
	//Get string edges
	try
	{
		string strefile = starting_tissue_file + STRINGEDGE_FILE_EXTENSION;
		ifstream fin_stre(strefile);
		if (fin_stre.good())
		{
			if (REPORT_OUT > 0)
				cout << "reading .stre file...\n";
			initialize_stringedges(fin_stre);
			if (REPORT_OUT > 0)
				cout << ".stre file read...\n";
		}
		else
		{
			cout << "No stre file\n";
		}
		fin_stre.close();
	}
	catch (const char *msg)
	{
		if (REPORT_OUT > 0){
			cout << msg << endl;
			cout << "No stre file\n";
		}
	}
	//No file with parameters, therefore use constants defined in VertexSystem.h.
	//Also initializes cell area and vertex energy

	setMinAndMaxPositions();
	setHingeMinAndMaxPositions();
	set_default_params();
	
	if (REPORT_OUT > 0)
		cout << "parameters set to cells\n";
	if (REPORT_OUT > 1)
		cout << T1_HEADER << endl << DIV_HEADER << endl << T2_HEADER << endl;

}

void Tissue::set_default_simulation_params()
{
	t1_active = T1_ACTIVE;
	t1_inwards_active = T1_BORDER_INWARDS_ACTIVE;
	t1_outwards_active = T1_BORDER_OUTWARDS_ACTIVE;
	division_active = DIVISION_ACTIVE;
	t2_active = T2_ACTIVE;
	join_edges_active = JOIN_EDGES_ACTIVE;
	check_if_edges_cross_opt = CHECK_EDGES_CROSS_AFTER_MOVE;
	control_cells_2sides = CONTROL_CELLS_2SIDES;
	integration_mode = INTEGR_MONTECARLO;
	min_range_vertex_movement = MIN_RANGE_VERTEX_MOVEMENT;
	max_range_vertex_movement = MAX_RANGE_VERTEX_MOVEMENT;
	h = DEFAULT_H;
	temperature_positive_energy = TEMPERATURE_POSITIVE_ENERGY;
	temperature_negative_energy = TEMPERATURE_NEGATIVE_ENERGY;
	temperature_means_proportion_of_acceptance = TEMPERATURE_MEANS_PROPORTION_OF_ACCEPTANCE;

	autonomous_cell_cycle = AUTONOMOUS_CELL_CYCLE;
	cell_cycle_controls_size = CELL_CYCLE_CONTROLS_SIZE;
	start_cell_cycle_at_random = true;
	keep_area_after_division = KEEP_AREA_AFTER_DIVISION;
	time_controls_size = TIME_CONTROLS_SIZE;
	time_decrease_exponent = TIME_DECREASE_EXPONENT;
	xcoord_controls_size = XCOORD_CONTROLS_SIZE;
	use_blade_area_for_coord_gradient = false;
	xcoord_decrease_exponent = XCOORD_DECREASE_EXPONENT;
	K_gradient_x = false;
	K_gradient_y = false;
	K_grad_exponent = TIME_DECREASE_EXPONENT;

	energy_term1 = ENERGY_TERM1;
	energy_term2 = ENERGY_TERM2;
	energy_term3 = ENERGY_TERM3;
	//energy_term4 = ENERGY_TERM4;
	//difference_flow_rate = 0;

	for(int i = 0; i < NUM_CELL_TYPES; i++) line_tension.val[i] = i % 2 == 0 ? LINE_TENSION_BLADE : LINE_TENSION_HINGE;
	for(int i = 0; i < NUM_CELL_TYPES; i++) line_tension_tissue_boundary.val[i] = LINE_TENSION_TISSUE_BOUNDARY;
	for(int i = 0; i < NUM_CELL_TYPES; i++) perimeter_contract.val[i] = i % 2 == 0 ? PERIMETER_CONTRACT_BLADE : PERIMETER_CONTRACT_HINGE;

	perimeter_contract_final = perimeter_contract;
	xcoord_controls_perim = false;
	ycoord_controls_perim = false;
	time_controls_perim = false;
	t1_transition_critical_distance = T1_TRANSITION_CRITICAL_DISTANCE;
	t2_transition_critical_area = T2_TRANSITION_CRITICAL_AREA;
	max_edge_length = MAX_EDGE_LENGTH;

	for(int i = 0; i < NUM_CELL_TYPES; i++) K.val[i] =  DEFAULT_K;
	for(int i = 0; i < NUM_CELL_TYPES; i++) K_final.val[i] =  DEFAULT_K;
	for(int i = 0; i < NUM_CELL_TYPES; i++) preferred_area_initial.val[i] =  PREFERRED_AREA_INITIAL;
	for(int i = 0; i < NUM_CELL_TYPES; i++) preferred_area_initial_gradient.val[i] =  PREFERRED_AREA_INITIAL;
	for(int i = 0; i < NUM_CELL_TYPES; i++) preferred_area_final.val[i] = PREFERRED_AREA_FINAL;
	for(int i = 0; i < NUM_CELL_TYPES; i++) division_angle_random_noise.val[i] = DIVISION_ANGLE_RANDOM_NOISE;
	for(int i = 0; i < NUM_CELL_TYPES; i++) division_angle_longest_axis.val[i] = DIVISION_ANGLE_LONGEST_AXIS;
	for(int i = 0; i < NUM_CELL_TYPES; i++) division_angle_external.val[i] = DIVISION_ANGLE_EXTERNAL;
	for(int i = 0; i < NUM_CELL_TYPES; i++) division_angle_external_degrees.val[i] = DIVISION_ANGLE_EXTERNAL_DEGREES;
	for(int i = 0; i < NUM_SPRING_TYPES; i++) spring_type_constants.val[i] = SPRING_CONSTANT;
	for(int i = 0; i < NUM_CELL_TYPES; i++) cell_cycle_limit.val[i] = CELL_CYCLE_LIMIT;
	for(int i = 0; i < NUM_CELL_TYPES; i++) vary_line_tension.val[i] = -1.0;
	for(int i = 0; i < NUM_CELL_TYPES; i++) line_tension.val[i] = LINE_TENSION_BLADE;
	for(int i = 0; i < NUM_CELL_TYPES; i++) line_tension.val[i] = LINE_TENSION_BLADE;

	vary_edge_tension_with_time = false;
	vary_edge_tension_time_exponent = 0;
	edge_angle_prop_external = vary_line_tension;
	edge_angle_prop_uniform = vary_line_tension;
	edge_angle_prop_maxangle = vary_line_tension;
	edge_angle_prop_random = vary_line_tension;
	edge_tension_external = vary_line_tension;
	edge_maxangle = vary_line_tension;
	edge_spatialmax_tension = vary_line_tension;
	edge_spatialmin_tension = vary_line_tension;
	edge_temporal_angle_efect_max = vary_line_tension;
	edge_temporal_angle_efect_min = vary_line_tension;

	spring_tension_mode = 0; //0: with tension for each type, 1: A-P compartments, 2: PD gradient, 3: AP compartments and PD gradient, 4:Different gradients in A or P
	posterior_comparment_region = 0; //used if spring_tension_mode is 1 or 3
	spring_posterior_compartment_factor = 0; //used if spring_tension_mode is 1 or 3
	spring_gradient_min_tension = 0; //used if spring_tension_mode is 2 or 3
	spring_gradient_max_tension = 0; //used if spring_tension_mode is 2 or 3
	spring_gradient_exponent = 0; //used if spring_tension_mode is 2 or 3
	spring_gradient_min_tension_P = 0; //used if spring_tension_mode is 4
	spring_gradient_max_tension_P = 0; //used if spring_tension_mode is 4
	spring_gradient_exponent_P = 0; //used if spring_tension_mode is 4
	mode_to_order_springs_PD = 0; //used if spring_tension_mode is 2 or 3
	include_hinge_in_spring_gradient = false;

	line_tension_interstatic  = LINE_TENSION_INTERSTATIC;
	hinge_blade_interface_tension = HINGE_BLADE_INTERFACE_TENSION;
	for(int i = 0; i < NUM_SPRING_TYPES; i++) spring_type_min_positions.val[i] = 0.5;

	use_term4 = USE_TERM4;
	energy_term4_posterior = vary_line_tension;
	energy_term4_anterior = vary_line_tension;

	add_static_to_hinge = -1;
	line_tension_interstatic = LINE_TENSION_INTERSTATIC;

	line_tension_stringedge = LINE_TENSION_TISSUE_BOUNDARY;
	tension_stringedge_posterior_prop = 1;
	string_edge_tension_min = LINE_TENSION_TISSUE_BOUNDARY;
	string_edge_tension_exponent = 1;
	string_anterior_gradient = false;
	string_posterior_gradient = false;
	string_distal_transition_tension = LINE_TENSION_TISSUE_BOUNDARY;
	set_hinge_string_tension = false;
	hinge_string_tension = LINE_TENSION_TISSUE_BOUNDARY;
	string_distal_transition_prop = 0.0;
	reference_for_gradient = REFERENCE_FOR_GRADIENT;
	wing_proportion_in_gradient = DEFAULT_PROPORTION_FOR_GRADIENT;
	string_equilibrium_distance = STRING_EQUILIBRIUM_DISTANCE;
	random_seed = RANDOM_SEED;
}

void Tissue::readNewParameters(std::string filename)
{
	initialize_params(filename);
	setMinAndMaxPositions(); //uses vertices
	setHingeMinAndMaxPositions(); //calculates and uses centroids, only hinge
	set_default_params(); //sets also cell centroids, in all cells
}

double Tissue::read_real_par(std::vector<std::string>::iterator &it)
{
	cout << "Reading real parameter... " << endl;
	cout << *it << endl;
	while (it->at(0) != '>'){
		it++;
		cout << "line: ";
		cout << *it << endl;
	}
		
	it++;
	double res = stod(*it);
	cout << res << endl;
	return res;
}
double Tissue::read_longint_par(std::vector<std::string>::iterator &it)
{
	while (it->at(0) != '>')
		it++;
	it++;
	long unsigned int res = stoul(*it);
	//cout << res << endl;
	return res;
}

cell_type_param Tissue::read_celltype_par(std::vector<std::string>::iterator &it, std::string::size_type sz)
{
	std::string s;
	int celltype;
	double auxd;
	static cell_type_param res;

	while (it->at(0) != '>')
		it++;
	it++;

	while (it->at(0) != '<')
	{
		s = *it;
		//cout << "s: " << s <<", ";
		celltype = stoi(s, &sz);
		s = s.substr(sz);
		s.erase(0, 1);
		auxd = stod(s, &sz);
		//cout << "celltype: " << celltype << ", val: " << auxd << endl;
		res.val[celltype] =auxd;
		it++;
	}
	//cout << "Inside: " << res.val[0] << " "<< res.val[1] << " "<< res.val[2] << " "<< res.val[3] << " \n";
	return res;
}

spring_type_param Tissue::read_springtype_par(std::vector<std::string>::iterator &it, std::string::size_type sz)
{
	std::string s;
	int sprtype;
	double auxd;
	spring_type_param res;

	while (it->at(0) != '>')
		it++;
	it++;

	while (it->at(0) != '<')
	{
		s = *it;
		sprtype = stoi(s, &sz);
		s = s.substr(sz);
		s.erase(0, 1);
		auxd = stod(s, &sz);
		res.val[sprtype] = auxd;
		it++;
	}
	return res;
}

void Tissue::initialize_params(std::string params_file)
{

	params_file += PARAMS_FILE_EXTENSION;
	std::ifstream fin;
	fin.open(params_file);
	std::string line;
	std::vector<std::string> inp;
	std::string::size_type sz;

	while (getline(fin, line))
		if (!(line.empty() || line.find_first_not_of(' ') == std::string::npos))
			if (line.at(0) != '#')
				inp.push_back(line);
	std::vector<std::string>::iterator it = inp.begin();
	//READ PARAMETERS IN ORDER FROM HERE:
	cout << "Param file reading... breakpoint 1" << endl;
	std::flush(cout);
	t1_active = read_real_par(it) > 0;
	t1_inwards_active = read_real_par(it) > 0;
	t1_outwards_active = read_real_par(it) > 0;
	division_active = read_real_par(it) > 0;
	t2_active = read_real_par(it) > 0;
	join_edges_active = read_real_par(it) > 0;
	control_cells_2sides = static_cast<bool>(read_real_par(it));
	check_if_edges_cross_opt = read_real_par(it) > 0;
	cout << "Param file reading... breakpoint 1b" << endl;
	std::flush(cout);
	long unsigned int temp_num_moves = read_longint_par(it);
	upper_bound_movements = max_accepted_movements;
	if(temp_num_moves > 0)
		max_accepted_movements += temp_num_moves;
	long unsigned int temp_write_freq = read_longint_par(it);
	if(temp_write_freq > 0)
		write_every_N_moves = temp_write_freq;
	long unsigned int temp_bound = read_longint_par(it);
	if(temp_bound > 0)
		upper_bound_movements += temp_bound;
	integration_mode = static_cast<int>(read_real_par(it));
	min_range_vertex_movement = read_real_par(it);
	max_range_vertex_movement = read_real_par(it);
	h = read_real_par(it);
	temperature_positive_energy = read_real_par(it);
	temperature_negative_energy = read_real_par(it);
	temperature_means_proportion_of_acceptance = read_real_par(it) > 0;
cout << "Param file reading... breakpoint 1c" << endl;
	std::flush(cout);
	energy_term1 = read_real_par(it);
	energy_term2 = read_real_par(it);
	energy_term3 = read_real_par(it);
	double sum_terms = (energy_term1 + energy_term2 + energy_term3)/3;
	energy_term1 /= sum_terms;
	energy_term2 /= sum_terms;
	energy_term3 /= sum_terms;
	//energy_term4 = read_real_par(it);
	cout << "Param file reading... breakpoint 2" << endl;
	std::flush(cout);
	//spring_constant = read_real_par(it);
	t1_transition_critical_distance = read_real_par(it);
	length_rotated_edge = read_real_par(it)*t1_transition_critical_distance;
	t2_transition_critical_area = read_real_par(it);
	//max_cell_area = read_real_par(it);
	max_edge_length = read_real_par(it);
	autonomous_cell_cycle = read_real_par(it) > 0;
	start_cell_cycle_at_random = read_real_par(it) > 0;
	cell_cycle_controls_size = read_real_par(it) > 0;
	keep_area_after_division = read_real_par(it) > 0;
	time_controls_size = read_real_par(it) > 0;
	time_controls_size = cell_cycle_controls_size ? false : time_controls_size;
	time_decrease_exponent = read_real_par(it);
	int aux = read_real_par(it);
	xcoord_controls_size = (aux == 1 || aux == 3);
	ycoord_controls_size = (aux == 2 || aux == 3);
	gradient_with_same_final_area = aux == 4;
	use_blade_area_for_coord_gradient = read_real_par(it) > 0;
	xcoord_decrease_exponent = read_real_par(it);
	aux = read_real_par(it);
	K_gradient_x = (aux == 1 || aux == 3);
	K_gradient_y = (aux == 2 || aux == 3);
	K_grad_exponent = read_real_par(it);
	cout << "Param file reading... breakpoint 3" << endl;
	std::flush(cout);
	K = read_celltype_par(it, sz);
	K_final = read_celltype_par(it, sz);
	line_tension = read_celltype_par(it, sz);
	line_tension_tissue_boundary = read_celltype_par(it, sz);

	time_controls_perim = read_real_par(it) > 0;
	aux = read_real_par(it);
	xcoord_controls_perim = (aux == 1 || aux == 3);
	ycoord_controls_perim = (aux == 2 || aux == 3);
	cout << "B" << endl;
	perimeter_contract = read_celltype_par(it, sz);
	perimeter_contract_final = read_celltype_par(it, sz);
	preferred_area_initial = read_celltype_par(it, sz);
	preferred_area_initial_gradient = read_celltype_par(it, sz);
	preferred_area_final = read_celltype_par(it, sz);
	division_angle_random_noise = read_celltype_par(it, sz);
	division_angle_longest_axis = read_celltype_par(it, sz);
	division_angle_external = read_celltype_par(it, sz);
	division_angle_external_degrees = read_celltype_par(it, sz);
	spring_type_constants = read_springtype_par(it, sz);
	max_cell_area = read_celltype_par(it, sz);
	cell_cycle_limit = read_celltype_par(it, sz);
	cout << "Param file reading... breakpoint 4" << endl;
	std::flush(cout);
	vary_line_tension = read_celltype_par(it, sz);
	vary_edge_tension_with_time = read_real_par(it) > 0;
	vary_edge_tension_time_exponent = read_real_par(it);
	edge_angle_prop_external = read_celltype_par(it, sz);
	edge_angle_prop_uniform = read_celltype_par(it, sz);
	edge_angle_prop_maxangle = read_celltype_par(it, sz);
	edge_angle_prop_random = read_celltype_par(it, sz);
	edge_tension_external = read_celltype_par(it, sz); //Initial values, then gene expression will change it in cells
	edge_maxangle = read_celltype_par(it, sz);
	edge_spatialmax_tension = read_celltype_par(it, sz);
	edge_spatialmin_tension = read_celltype_par(it, sz);
	edge_temporal_angle_efect_max = read_celltype_par(it, sz);
	edge_temporal_angle_efect_min = read_celltype_par(it, sz);

	spring_type_min_positions = read_springtype_par(it, sz);
	add_static_to_hinge = read_real_par(it);
	//difference_flow_rate = read_real_par(it);
	std::cout << "Param file reading... breakpoint 5" << endl;
	std::flush(cout);
	spring_tension_mode = read_real_par(it); //0: with tension for each type, 1: A-P compartments, 2: PD gradient, 3: AP compartments and PD gradient, 4:Different gradients in A or P
	posterior_comparment_region = read_real_par(it); //used if spring_tension_mode is 1 or 3
	spring_posterior_compartment_factor = read_real_par(it); //used if spring_tension_mode is 1 or 3
	spring_gradient_min_tension = read_real_par(it); //used if spring_tension_mode is 2 or 3
	spring_gradient_max_tension = read_real_par(it); //used if spring_tension_mode is 2 or 3
	spring_gradient_exponent = read_real_par(it); //used if spring_tension_mode is 2 or 3
	spring_gradient_min_tension_P = read_real_par(it); //used if spring_tension_mode is 4
	spring_gradient_max_tension_P = read_real_par(it); //used if spring_tension_mode is 4
	spring_gradient_exponent_P = read_real_par(it); //used if spring_tension_mode is 4
	mode_to_order_springs_PD = read_real_par(it);
	include_hinge_in_spring_gradient = read_real_par(it) > 0;
	cout << "C" << endl;
	//For random ("active") T1 transitions
	active_t1_prob = read_real_par(it);
	min_angle_for_active_t1 = read_real_par(it);
	max_angle_for_active_t1 = read_real_par(it);
	minsin2rant1 = sin(M_PI * min_angle_for_active_t1/180);
	maxsin2rant1 = sin(M_PI * max_angle_for_active_t1/180);
	std::cout << "Param file reading... breakpoint 6" << endl;
	std::flush(cout);
	hinge_blade_interface_tension = read_real_par(it);
	line_tension_interstatic = read_real_par(it);
	cout << "Param file reading... breakpoint 7" << endl;
	std::flush(cout);
	use_term4 = read_real_par(it) > 0;
	energy_term4_anterior = read_celltype_par(it, sz);
	std::cout << "Param file reading... breakpoint 8" << endl;
	std::flush(cout);
	energy_term4_posterior = read_celltype_par(it, sz);
	line_tension_stringedge = read_real_par(it);
	tension_stringedge_posterior_prop = read_real_par(it);
	string_edge_tension_min =  read_real_par(it);
	string_edge_tension_exponent =  read_real_par(it);
	string_anterior_gradient = read_real_par(it) > 0;
	string_posterior_gradient = read_real_par(it) > 0;
	set_hinge_string_tension = read_real_par(it) > 0;
	hinge_string_tension = read_real_par(it);
	string_distal_transition_tension = read_real_par(it);
	std::cout << "Param file reading... breakpoint 9" << endl;
	std::flush(cout);
	string_distal_transition_prop = read_real_par(it);
	reference_for_gradient = read_real_par(it);
	wing_proportion_in_gradient = read_real_par(it);
	string_equilibrium_distance = read_real_par(it) > 0;

	momentum_term_cuticle = read_real_par(it);
	momentum_term_tissue = read_real_par(it);
	momentum_ponderate_current = read_real_par(it);

	random_seed = read_real_par(it);
	std::cout << "Param file reading... final" << endl;
	std::flush(cout);
	cout << "D" << endl;
}

/*
InitiDlizes vertices from file
Input: ifstreDm pointing to a file defining vertex coordinates
*/
void Tissue::initialize_vertices(std::ifstream &inp)
{
	string s;
	getline(inp, s);
	this->num_vertices = stoi(s);

	Vertex v;
	v.dead = false;
	//v.movable = true;
	v.spring = EMPTY_CONNECTION;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		v.cells[i] = EMPTY_CONNECTION;
		v.edges[i] = EMPTY_CONNECTION;
		v.neighbour_vertices[i] = EMPTY_CONNECTION;
	}

	std::string::size_type sz;
	while (getline(inp, s))
	{
		v.x = stod(s, &sz);
		s = s.substr(sz);
		v.y = stod(s, &sz);
		s = s.substr(sz);
		v.ind = stoi(s, &sz);
		s = s.substr(sz);
		v.movable = s.empty() ? true : stoi(s);
		v.movable_x = s.empty() ? true : stoi(s) == 1 || stoi(s) == 2;
		v.movable_y = s.empty() ? true : stoi(s) == 1 || stoi(s) == 3;
		v.moves_accepted = 1;
		v.moves_rejected = 1;
		//cout << "movable: " << s << ", " << static_cast<int>(v.movable) << static_cast<int>(v.movable_x) << static_cast<int>(v.movable_y) << endl;
		this->vertices.push_back(v);
	}
}

//Initializes springs from file
//Input: ifstream pointing to a file defining springs (pairs of vertex indices; first index connected to cell, second index outside cells)

void Tissue::initialize_springs(std::ifstream &inp)
{
	string s;

	getline(inp, s);
	this->num_springs = stoi(s);
	if (this->num_springs == 0)
		return;
	Edge e;
	e.dead = false;
	e.type = EdgeType::spring;
	e.cells[0] = EMPTY_CONNECTION;
	e.cells[1] = EMPTY_CONNECTION;
	e.tension = spring_type_constants.val[0];

	std::string::size_type sz;
	int i = 0;
	int sprtype = 0;
	while (getline(inp, s))
	{
		e.vertices[0] = stoi(s, &sz);
		s = s.substr(sz);
		e.vertices[1] = stoi(s, &sz);
		s = s.substr(sz);
		try
		{
			sprtype = stoi(s, &sz);
			e.tension = this->spring_type_constants.val[sprtype];
		}
		catch (const std::invalid_argument &ia)
		{
			std::cerr << "Invalid argument (No third spring column): " << ia.what() << '\n';
			e.tension = spring_type_constants.val[0];
		}
		e.ind = i;
		e.length = distance(e.vertices[0], e.vertices[1]);

		//Avoid errors in input file (some vertices may have more than one spring)
		if (vertices[e.vertices[0]].spring == EMPTY_CONNECTION && vertices[e.vertices[1]].spring == EMPTY_CONNECTION)
		{
			this->springs.push_back(e);
			vertices[e.vertices[0]].spring = e.ind;
			vertices[e.vertices[1]].spring = e.ind;
			i++;
		}
		else
		{
			this->num_springs--; //The input file has more springs
		}
	}

}

void Tissue::setSpringTension(){
		//Now set spring tension in case it is not directly set by spring_type_constants
	cout << "Setting spring tension... \n";
	switch (spring_tension_mode)
	{
	case 1:
		if(REPORT_OUT > 0) cout << "Spring tension mode : 1\n";
		setSpringTension_mode1();
		break;
	case 2:
		if(REPORT_OUT > 0) cout << "Spring tension mode : 2\n";
		setSpringTension_mode2();
		break;
	case 3:
		if(REPORT_OUT > 0) cout << "Spring tension mode : 3\n";
		setSpringTension_mode3();
		break;
	case 4:
		if(REPORT_OUT > 0) cout << "Spring tension mode : 4\n";
		setSpringTension_mode4();
		break;
	case 0:
		if(REPORT_OUT > 0) cout << "Spring tension mode : 0\n";
	default:
		break;
	}
	cout << "Spring tension set\n";
}

void Tissue::initialize_stringedges(std::ifstream& inp){
//Create new edges
	string s;
	getline(inp, s);
	int num_stringedges = stoi(s);
	if (num_stringedges == 0)
		return;
	Edge e;
	e.dead = false;
	e.type = EdgeType::stringedge;
	e.cells[0] = EMPTY_CONNECTION;
	e.cells[1] = EMPTY_CONNECTION;
	e.tension = line_tension_stringedge;
	e.can_transition = false;
	e.base_tension = e.tension;

	std::string::size_type sz;
	int i = num_edges;
	while (getline(inp, s))
	{
		e.vertices[0] = stoi(s, &sz);
		if(e.vertices[0] < 0)
			continue;
		s = s.substr(sz);
		e.vertices[1] = stoi(s, &sz);
		s = s.substr(sz);
		e.ind = i;
		e.length = distance(e.vertices[0], e.vertices[1]);
		vertices[e.vertices[0]].neighbour_vertices[which(EMPTY_CONNECTION, vertices[e.vertices[0]].neighbour_vertices, CELLS_PER_VERTEX)] = e.vertices[1];
		vertices[e.vertices[1]].neighbour_vertices[which(EMPTY_CONNECTION, vertices[e.vertices[1]].neighbour_vertices, CELLS_PER_VERTEX)] = e.vertices[0];
		vertices[e.vertices[0]].edges[which(EMPTY_CONNECTION, vertices[e.vertices[0]].edges, CELLS_PER_VERTEX)] = i;
		vertices[e.vertices[1]].edges[which(EMPTY_CONNECTION, vertices[e.vertices[1]].edges, CELLS_PER_VERTEX)] = i;

		this->edges.push_back(e);
		i++;
	}
}
/*
Initializes cells from file; adds cells to vertices
Input: ifstream pointing to a file defining vertex identifiers (starting from 0) for each cell
*/
void Tissue::initialize_cells(std::ifstream &inp)
{
	string s, s2;

	std::string::size_type sz;
	getline(inp, s);
	this->num_cells = stoi(s, &sz); //stoi(s.substr(0, split_ind));
	if (stoi(s.substr(sz)) > MAX_SIDES_PER_CELL)
	{
		throw "Cell file has more vertices by cell than allowed by MAX_SIDES_PER_CELL constant";
	}
	int vertex, vert_count;

	Cell c;
	c.dead = false;
	c.can_divide = !this->autonomous_cell_cycle;
	c.centroid_x = 0;
	c.centroid_y = 0;
	c.num_divisions = 0;
	int cell_count = 0;
	while (getline(inp, s))
	{					//Each row in file represents a cell
		vert_count = 0; //Store number of vertices read so far for this cell
		c.ind = cell_count;
		c.num_vertices = 0;
		for (int i = 0; i < MAX_SIDES_PER_CELL; i++)
		{ //Reset all connections of cell c to EMPTY_CONNECTIONS
			c.vertices[i] = EMPTY_CONNECTION;
			c.edges[i] = EMPTY_CONNECTION; //edges are initialized here but will be added later, in method initialize_edges
		}

		vertex = stoi(s, &sz);
		while (vertex != EMPTY_CONNECTION)
		{
			c.vertices[vert_count] = vertex;
			addCellToVertex(vertex, cell_count);
			vert_count++;
			s = s.substr(sz);
			vertex = stoi(s, &sz);
		}
		c.num_vertices = vert_count;
		while (vertex == EMPTY_CONNECTION)
		{ //Last number in row, if present, represents cell type
			s = s.substr(sz);
			if (!s.empty())
			{
				vertex = stoi(s, &sz);
			}
			else
			{
				break;
			}
		}
		c.type = vertex == EMPTY_CONNECTION ? CellType::blade : static_cast<CellType>(vertex);
		c.cell_cycle_state = start_cell_cycle_at_random ? std::rand() % static_cast<int>(cell_cycle_limit.val[static_cast<int>(c.type)]) : 0;
		this->cells.push_back(c);
		cell_count++;
	} //while rows in file
}

/*
Initializes edges from vertex and cell data already read from file
Assumes that vertices in cell are in order (i.e., consecutive vertices are bound by an edge, and last vertex is bound to first)
Calls: Tissue::addEdge
*/

void Tissue::initialize_edges()
{

	int previous_vertex;
	for (int cell = 0; cell < this->num_cells; cell++)
	{ //for each cell
		previous_vertex = this->cells[cell].vertices[this->cells[cell].num_vertices - 1];
		for (int v = 0; v < this->cells[cell].num_vertices; v++)
		{ //For each vertex in a cell
			addEdge(previous_vertex, this->cells[cell].vertices[v], cell);
			previous_vertex = this->cells[cell].vertices[v];
		} //for vertex in cell
	}	  //for each cell
}
/*
Initializes a single edge. If edge between two vertices already exists, adds the second cell.
If the edge does not exist, it creates one and adds the two vertices and the first cell.
Inputs:
index of vertex1
index ofvertex2
index of cell to which v1 and v2 belong
*/
void Tissue::addEdge(int v1, int v2, int cell)
{
	bool edge_found = false;
	//If edge between v1 and v2 already exists, add Cell to it
	for (int i = 0; i < this->num_edges; i++)
	{
		if ((this->edges[i].vertices[0] == v1 && this->edges[i].vertices[1] == v2) ||
			(this->edges[i].vertices[1] == v1 && this->edges[i].vertices[0] == v2))
		{
			this->edges[i].cells[1] = cell;
			addEdgeToCell(cell, i);
			edge_found = true;
			break;
		}
	}
	//If edge does not exist, create edge and add the two vertices and the cell passed as argument to it. The other cell will be added later
	if (!edge_found)
	{
		Edge e;
		e.dead = false;
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

void Tissue::set_default_params()
{
	for (int c = 0; c < cells.size(); c++)
	{
		int celltype = static_cast<int>(cells[c].type);
		cells[c].perimeter_contractility = perimeter_contract.val[celltype];
		//cells[c].preferred_area = preferred_area_initial.val[celltype];
		cells[c].area = calculateCellArea(cells[c]);
		cells[c].perimeter = calculateCellPerimeter(cells[c]);
		calculateCellCentroid(cells[c]); 
		calculateBasePrefAreaAndPerim(cells[c]);
		setKWithCoords(c);
		cells[c].father = c;
		cells[c].x_at_start = cells[c].centroid_x;

		cells[c].division_angle_random_noise = division_angle_random_noise.val[celltype];
		cells[c].division_angle_longest = division_angle_longest_axis.val[celltype];
		cells[c].division_angle_external = division_angle_external.val[celltype];
		cells[c].division_angle_external_degrees = division_angle_external_degrees.val[celltype];

		cells[c].cell_cycle_limit = cell_cycle_limit.val[celltype];
		cells[c].max_area = max_cell_area.val[celltype];

		cells[c].vary_line_tension = vary_line_tension.val[celltype] > 0;
		cells[c].edge_angle_prop_external = edge_angle_prop_external.val[celltype];
		cells[c].edge_angle_prop_uniform = edge_angle_prop_uniform.val[celltype];
		cells[c].edge_angle_prop_maxangle = edge_angle_prop_maxangle.val[celltype];
		cells[c].edge_angle_prop_random = edge_angle_prop_random.val[celltype];
		cells[c].edge_tension_external = edge_tension_external.val[celltype];
		cells[c].edge_maxangle = edge_maxangle.val[celltype];
		cells[c].edge_spatialmax_tension = edge_spatialmax_tension.val[celltype];
		cells[c].edge_spatialmin_tension = edge_spatialmin_tension.val[celltype];
	}

	for (int e = 0; e < edges.size(); e++)
	{
		if(edges[e].type != EdgeType::stringedge){
			setEdgeType(e);
			setEdgeTension(e);
		}else{
			setStringTension(e);
		}
	}

	setSpringTension();

	for (int v = 0; v < vertices.size(); v++)
	{
		vertices[v].x_drag = vertices[v].x;
		vertices[v].y_drag = vertices[v].y;
		vertices[v].energy = calculateEnergy(vertices[v]);
	}
}

void Tissue::setStringTension(int e){
			double position = 0.5*(vertices[edges[e].vertices[0]].x + vertices[edges[e].vertices[1]].x);
			double position_factor = (position - min_xpos)/(max_xpos - min_xpos);
			if((1 - position_factor) < string_distal_transition_prop && string_distal_transition_prop > 0){
				//cout << "A\n";
				//Transition zone between static vertices at the tip and the rest of the cuticle
				edges[e].tension = string_distal_transition_tension;
			}else if(set_hinge_string_tension && position <= hinge_max_xpos){
				//cout << "B\n";
				//Set a specific tension to hinge
				edges[e].tension = hinge_string_tension;
				if(0.5*(vertices[edges[e].vertices[0]].y + vertices[edges[e].vertices[1]].y) < AP_compartment_limit ){
					//Posterior hinge tension is also multiplied by some proportion
					edges[e].tension *= tension_stringedge_posterior_prop;
				}
			}else{
				//Set string tension in the rest of edges
				position_factor = expAdvance(position_factor, string_edge_tension_exponent);
				//position_factor =  0.5*(vertices[edges[e].vertices[0]].x + vertices[edges[e].vertices[1]].x) >= 600?  1: 0;// PROVISIONAL
				//cout << "C " << position_factor << endl;
				double pos_tension = string_edge_tension_min + (line_tension_stringedge - string_edge_tension_min)*position_factor;
				if(0.5*(vertices[edges[e].vertices[0]].y + vertices[edges[e].vertices[1]].y) >= AP_compartment_limit ){
					//If in anterior compartment
					if(string_anterior_gradient){ 
						edges[e].tension = pos_tension; //Redundant at start, but useful when changing parameter files
					}else{ 
						edges[e].tension = line_tension_stringedge; 
					}
				}else{
					//In posterior compartment
					if(momentum_term_cuticle > 0){
						float pos = (position - min_xpos)/(max_xpos - min_xpos); //Same as position_factor before applying exp
						edges[e].tension = line_tension_stringedge * (tension_stringedge_posterior_prop + (1-tension_stringedge_posterior_prop)*position_factor);
					}else if(string_posterior_gradient){
						edges[e].tension = pos_tension*tension_stringedge_posterior_prop;
					}else{
						edges[e].tension = line_tension_stringedge * tension_stringedge_posterior_prop;
					}
				}
			}
			edges[e].base_tension = edges[e].tension;
			edges[e].optimal_length = edges[e].length;
}

void Tissue::setEdgeType(int e)
{
	int c1 = edges[e].cells[0];
	int c2 = edges[e].cells[1];
	int aux;
	if (contains(EMPTY_CONNECTION, edges[e].cells, VERTEX_PER_EDGE))
	{
		aux = edges[e].cells[0] == EMPTY_CONNECTION ? edges[e].cells[1] : edges[e].cells[0];
		edges[e].type = EdgeType::tissue_boundary;
		if(!(vertices[edges[e].vertices[0]].movable_x && vertices[edges[e].vertices[0]].movable_y) ){
			edges[e].tension = line_tension_interstatic;
			edges[e].can_transition = false;
		}else{
			edges[e].tension = line_tension_tissue_boundary.val[static_cast<int>(cells[aux].type)];
			edges[e].can_transition = true;
		}
	}
	else if (cells[c1].type == CellType::blade && cells[c2].type == CellType::blade)
	{
		edges[e].type = EdgeType::blade;
		edges[e].tension = line_tension.val[static_cast<int>(CellType::blade)];
		edges[e].can_transition = true;
	}
	else if (cells[c1].type == CellType::hinge && cells[c2].type == CellType::hinge)
	{
		edges[e].type = EdgeType::hinge;
		edges[e].tension = line_tension.val[static_cast<int>(CellType::hinge)];
		edges[e].can_transition = true;
	}
	else if ((cells[c1].type == CellType::hinge && cells[c2].type == CellType::blade) ||
	(cells[c1].type == CellType::blade && cells[c2].type == CellType::hinge) ||
	(cells[c1].type == CellType::vein_blade && cells[c2].type == CellType::vein_hinge) ||
	(cells[c1].type == CellType::vein_hinge && cells[c2].type == CellType::vein_blade) ||
	(cells[c1].type == CellType::blade && cells[c2].type == CellType::vein_hinge) ||
	(cells[c1].type == CellType::vein_hinge && cells[c2].type == CellType::blade) ||
	(cells[c1].type == CellType::vein_blade && cells[c2].type == CellType::hinge) ||
	(cells[c1].type == CellType::hinge && cells[c2].type == CellType::vein_blade))
	{
		edges[e].type = EdgeType::hinge;
		edges[e].tension = hinge_blade_interface_tension;
		edges[e].can_transition = false;
	}
	else if ((cells[c1].type == CellType::blade && cells[c2].type == CellType::vein_blade) ||
			 (cells[c2].type == CellType::blade && cells[c1].type == CellType::vein_blade) ||
			 (cells[c1].type == CellType::blade && cells[c2].type == CellType::vein_hinge) ||
			 (cells[c2].type == CellType::blade && cells[c1].type == CellType::vein_hinge))
	{
		edges[e].type = EdgeType::vein_blade;
		edges[e].tension = line_tension.val[static_cast<int>(CellType::vein_blade)];
		edges[e].can_transition = true; //false; true
	}
	else if ((cells[c1].type == CellType::hinge && cells[c2].type == CellType::vein_hinge) ||
			 (cells[c2].type == CellType::hinge && cells[c1].type == CellType::vein_hinge) ||
			 (cells[c1].type == CellType::hinge && cells[c2].type == CellType::vein_blade) ||
			 (cells[c2].type == CellType::hinge && cells[c1].type == CellType::vein_blade))
	{
		edges[e].type = EdgeType::vein_hinge;
		edges[e].tension = line_tension.val[static_cast<int>(CellType::vein_hinge)];
		edges[e].can_transition = true; //false; true
	}
	else if (cells[c1].type == CellType::vein_blade && cells[c2].type == CellType::vein_blade)
	{
		edges[e].type = EdgeType::vein;
		edges[e].tension = line_tension.val[static_cast<int>(CellType::blade)]; //vein_blade is only for vein border
		edges[e].can_transition = true;					  //false;true
	}
	else if (cells[c1].type == CellType::vein_hinge && cells[c2].type == CellType::vein_hinge)
	{
		edges[e].type = EdgeType::vein;
		edges[e].tension = line_tension.val[static_cast<int>(CellType::hinge)]; //vein_hinge is only for vein border
		edges[e].can_transition = true;					  //false;true
	}
	else if ((cells[c1].type == CellType::vein_blade && cells[c2].type == CellType::vein_hinge) || (cells[c1].type == CellType::vein_hinge && cells[c2].type == CellType::vein_blade))
	{
		edges[e].type = EdgeType::vein;
		edges[e].tension = 0.5 * (line_tension.val[static_cast<int>(CellType::hinge)] + line_tension.val[static_cast<int>(CellType::blade)]); //vein_hinge and vein_blade are only for vein border
		edges[e].can_transition = false;
	}
	edges[e].base_tension = edges[e].tension; //whereas tension can change, base_tension can't
}

void Tissue::setEdgeTension(int e)
{
	if (edges[e].type == EdgeType::tissue_boundary || edges[e].type == EdgeType::stringedge)
		return;
	if (!(cells[edges[e].cells[0]].vary_line_tension && cells[edges[e].cells[1]].vary_line_tension))
		return;
	float mins, maxs, maxangle, angle, pmaxan, punif; //pex, prand, tensionext, tensionrand
	int cellvar;
	mins = 0.5 * (cells[edges[e].cells[0]].edge_spatialmin_tension + cells[edges[e].cells[1]].edge_spatialmin_tension);
	maxs = 0.5 * (cells[edges[e].cells[0]].edge_spatialmax_tension + cells[edges[e].cells[1]].edge_spatialmax_tension);
	maxangle = 0.5 * (cells[edges[e].cells[0]].edge_maxangle + cells[edges[e].cells[1]].edge_maxangle);

	if (vary_edge_tension_with_time)
	{
		//If proportion determined by angle (pmaxan) varies with time, then proportion determined by angle is a value between
		// mint and maxt, which depends on current time time
		double time_factor = static_cast<double>(counter_moves_accepted) / upper_bound_movements;
		time_factor = time_factor >= 1? 1.0 : expAdvance(time_factor, vary_edge_tension_time_exponent);
		//time_factor = expAdvance(time_factor, vary_edge_tension_time_exponent);
		double mint = (edge_temporal_angle_efect_min.val[static_cast<int>(cells[edges[e].cells[0]].type)] + 
		edge_temporal_angle_efect_min.val[static_cast<int>(cells[edges[e].cells[1]].type)]) * 0.5;
		double maxt = (edge_temporal_angle_efect_max.val[static_cast<int>(cells[edges[e].cells[0]].type)] + 
		edge_temporal_angle_efect_max.val[static_cast<int>(cells[edges[e].cells[1]].type)]) * 0.5;
		pmaxan = mint + (maxt - mint) * time_factor;
		punif = abs(1.0 - pmaxan);
/* 		if(counter_moves_accepted % 1000000 == 0) {
			cout << "time_factor: "<< time_factor << ", min_time_effect: " << mint << ", max_time_effect: " << maxt << ", so maxang effect is: " << pmaxan << endl;
		} */
	}
	else
	{
		//If proportion determined by angle does not vary with time, then pmaxan is coded directly as a parameter in cells touched by edge.
		pmaxan = 0.5 * (cells[edges[e].cells[0]].edge_angle_prop_maxangle + cells[edges[e].cells[1]].edge_angle_prop_maxangle);
		punif = abs(1.0 - pmaxan);
	} //end if time determines influence of angle

	maxangle *= M_PI / 180; //maxangle is coded in degrees
	angle = atan2(vertices[edges[e].vertices[1]].y - vertices[edges[e].vertices[0]].y, vertices[edges[e].vertices[1]].x - vertices[edges[e].vertices[0]].x); //angle of edge
	angle = abs(sin(0.5 * M_PI + abs(angle - maxangle)));																								 //Point of sin wave (from 0 to 1)
	angle = mins + angle * (maxs - mins);																													 //Value of tension at this point of sin wave
	edges[e].tension = (punif * edges[e].base_tension + pmaxan * angle ) / (punif + pmaxan);
}

// A-P compartments (uses tension determined by spring type)
void Tissue::setSpringTension_mode1(){
	//1) Calculate division point (Now it is calculated in the set min/max points function)
/* 	float miny = vertices[0].y, maxy = vertices[0].y;
	for(Vertex v: vertices){
		if(v.y > maxy) maxy = v.y;
		if(v.y < miny) miny = v.y;
	}
	AP_compartment_limit = miny + (maxy - miny)*posterior_comparment_region; */
	int vaux;
	for(Edge &s: springs){
		vaux = vertices[s.vertices[0]].movable ? s.vertices[1] : s.vertices[0];
		if(vertices[vaux].y <= AP_compartment_limit){
			s.tension*=spring_posterior_compartment_factor;
		}
	}
}
// P-D gradient
void Tissue::setSpringTension_mode2(){
	//AP_compartment_limit = 0;//To avoid future problems
	std::vector<int> vert_indices, spr_indices;
	std::vector<float> gradient_factor;
	int vaux, vcell;
	bool include[num_springs];
	for(Edge &s: springs){
		if(vertices[s.vertices[0]].movable){
			vaux = s.vertices[1];
			vcell = s.vertices[0];
		}else{
			 vaux = s.vertices[0];
			 vcell = s.vertices[1];
		}
		include[s.ind] = true;
		if(!include_hinge_in_spring_gradient){
 			for(int i=0; i < CELLS_PER_VERTEX; i++){
				if(vertices[vcell].cells[i] != EMPTY_CONNECTION){
					if(cells[vertices[vcell].cells[i]].type == CellType::hinge || cells[vertices[vcell].cells[i]].type == CellType::vein_hinge){
						include[s.ind] = false;
						break;
					}
				}
			}
		}
		if(include[s.ind]){
			vert_indices.push_back(vaux);
			spr_indices.push_back(s.ind);
		}
	}

	switch(mode_to_order_springs_PD){
		case 1:
			gradient_factor = getSpringGradientFactor_mode1(vert_indices);
			break;
		case 2:
			gradient_factor = getSpringGradientFactor_mode2(vert_indices);
			break;
		default:
			break;
	}
	for(int i = 0; i < spr_indices.size(); i++){
			//cout << i <<": tension before: " << springs[spr_indices[i]].tension;
			springs[spr_indices[i]].tension = spring_gradient_min_tension + (spring_gradient_max_tension - spring_gradient_min_tension)*expAdvance(gradient_factor[i], spring_gradient_exponent); //gradient_factor[i]; //*
			//cout << ", tension after: " << springs[spr_indices[i]].tension << ", grad fact: " << gradient_factor[i] << ", exp: " << expAdvance(gradient_factor[i], spring_gradient_exponent)<< endl;

	}
}
//P-D gradient multiplied by a factor in each compartment
void Tissue::setSpringTension_mode3(){
	setSpringTension_mode2();
	setSpringTension_mode1();
}
// A different (independent) gradient in each compartment
void Tissue::setSpringTension_mode4(){
/* 	//1) Calculate division point
	float miny = vertices[0].y, maxy = vertices[0].y;
	for(Vertex v: vertices){
		if(v.y > maxy) maxy = v.y;
		if(v.y < miny) miny = v.y;
	}
	AP_compartment_limit = miny + (maxy - miny)*posterior_comparment_region; */
	int vaux;
	std::vector<int> vert_indices_anterior, vert_indices_posterior, anterior_springs, posterior_springs;
	std::vector<float> gradient_factor_anterior, gradient_factor_posterior;
	for(Edge &s: springs){
		vaux = vertices[s.vertices[0]].movable ? s.vertices[1] : s.vertices[0];
		if(vertices[vaux].y <= AP_compartment_limit){
			vert_indices_posterior.push_back(vaux);
			posterior_springs.push_back(s.ind);
		}else{
			vert_indices_anterior.push_back(vaux);
			anterior_springs.push_back(s.ind);
		}
	}
	switch (mode_to_order_springs_PD)
	{
	case 1:
		gradient_factor_anterior = getSpringGradientFactor_mode1(vert_indices_anterior);
		gradient_factor_posterior = getSpringGradientFactor_mode1(vert_indices_posterior);
		break;
	case 2:
		gradient_factor_anterior = getSpringGradientFactor_mode2(vert_indices_anterior);
		gradient_factor_posterior = getSpringGradientFactor_mode2(vert_indices_posterior);
		break;
	default:
		break;
	}
	float aux;
	int spr;
	for(int i = 0; i < anterior_springs.size(); i++){
		spr = anterior_springs[i];
		aux = expAdvance(gradient_factor_anterior[i], spring_gradient_exponent);
		aux=gradient_factor_anterior[i];
		springs[spr].tension = spring_gradient_min_tension + (spring_gradient_max_tension - spring_gradient_min_tension)*aux;
	}
	for(int i = 0; i < posterior_springs.size(); i++){
		spr = posterior_springs[i];
		gradient_factor_anterior[i];
		aux = expAdvance(gradient_factor_posterior[i], spring_gradient_exponent_P);
		springs[spr].tension = spring_gradient_min_tension_P + (spring_gradient_max_tension_P - spring_gradient_min_tension_P)*aux;
	}
}
//Factor = proportional to position in x
std::vector<float> Tissue::getSpringGradientFactor_mode1(std::vector<int> &vert_indices){
	std::vector<float> factors;
	//1) Calculate min and max points to compare
	float minx = vertices[vert_indices[0]].x, maxx = vertices[vert_indices[0]].x, threshold;
	for(int i: vert_indices){
		if(vertices[i].x > maxx) maxx = vertices[i].x;
		if(vertices[i].x < minx) minx = vertices[i].x;
	}
	float aux;
	for(int i = 0; i < vert_indices.size(); i++ ){
		//aux = vertices[vert_indices[i]].x > 600? 1:0; // PROVISIONAL
		//factors.push_back(aux);
		factors.push_back((vertices[vert_indices[i]].x - minx)/(maxx - minx));
	}
	return factors;
}
//Factor = proportional to vertex order in x
std::vector<float> Tissue::getSpringGradientFactor_mode2(std::vector<int> &vert_indices){
	std::vector<int> ordered(vert_indices.size());
	std::vector<float> factor(vert_indices.size());
	//Initialize new vector
	for(int i = 0; i < vert_indices.size();i++){
		ordered[i] = vert_indices[i];
	}
	//Order
	int aux;
	for(int i = 0; i < ordered.size() - 1;i++){
			for(int j = i+1; j < ordered.size();j++){
				if(vertices[ordered[j]].x < vertices[ordered[i]].x){
					aux = ordered[j];
					ordered[j] = ordered[i];
					ordered[i] = aux;
				}
		}
	}
	//Get factor according to Order
	for(int i = 0; i < ordered.size();i++){
			for(int j = 0; j < ordered.size();j++){
				if(ordered[j] ==  vert_indices[i]){
					factor[i] = static_cast<float>(j)/ordered.size();
				}
		}
	}
	return factor;
}


/*
Adds a cell index to the vector of cells in a vertex structure
Input:
- vertex: the index of vertex to which cell must be added
- cell: the index of the cell to add
*/
void Tissue::addCellToVertex(int vertex, int cell)
{
	int i = 0;
	while (this->vertices[vertex].cells[i] != EMPTY_CONNECTION)
		i++;
	this->vertices[vertex].cells[i] = cell;
}

/*
Adds an edge index to the vector of edges in a vertex structure
Input:
- vertex: the index of vertex to which cell must be added
- edge: the index of the edge to add
*/
void Tissue::addEdgeToVertex(int vertex, int edge)
{
	int i = 0;
	while (this->vertices[vertex].edges[i] != EMPTY_CONNECTION)
		i++;
	this->vertices[vertex].edges[i] = edge;
}

/*
Adds an edge index to the vector of edges in a cell structure
Input:
- cell: the index of edge to which cell must be added
- edge: the index of the edge to add
*/
void Tissue::addEdgeToCell(int cell, int edge)
{
	int i = 0;
	while (this->cells[cell].edges[i] != EMPTY_CONNECTION)
		i++;
	this->cells[cell].edges[i] = edge;
}

void Tissue::addNeighbourVertex(int vertex1, int vertex2)
{
	int i = 0;
	while (this->vertices[vertex1].neighbour_vertices[i] != EMPTY_CONNECTION)
		i++;
	this->vertices[vertex1].neighbour_vertices[i] = vertex2;
	i = 0;
	while (this->vertices[vertex2].neighbour_vertices[i] != EMPTY_CONNECTION)
		i++;
	this->vertices[vertex2].neighbour_vertices[i] = vertex1;
}

//dead_vertices, dead_cells, dead_edges

int Tissue::newVertex()
{
	int v;
	if (false)
	{ //if (!dead_vertices.empty()){
		v = dead_vertices.front();
	}
	else
	{
		Vertex vert;
		vert.ind = static_cast<int>(this->vertices.size());
		v = vert.ind;
		this->vertices.push_back(vert);
	}
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		vertices[v].cells[i] = EMPTY_CONNECTION;
		vertices[v].edges[i] = EMPTY_CONNECTION;
		vertices[v].neighbour_vertices[i] = EMPTY_CONNECTION;
	}
	vertices[v].movable = true;
	vertices[v].movable_x = true;
	vertices[v].movable_y = true;
	vertices[v].spring = EMPTY_CONNECTION;
	vertices[v].dead = false;
	vertices[v].moves_accepted = 1;
	vertices[v].moves_rejected = 1;
	this->num_vertices++;
	//if (!dead_vertices.empty()) dead_vertices.pop();
	return v;
}

int Tissue::newVertex(double x, double y)
{
	int ind = newVertex();
	vertices[ind].x = x;
	vertices[ind].y = y;
	return ind;
}
int Tissue::newCell()
{
	int c;
	if (false)
	{ //if(!dead_cells.empty()){
		c = dead_cells.front();
	}
	else
	{
		Cell cl;
		cl.ind = static_cast<int>(this->cells.size());
		c = cl.ind;
		this->cells.push_back(cl);
	}
	for (int i = 0; i < MAX_SIDES_PER_CELL; i++)
	{
		cells[c].vertices[i] = EMPTY_CONNECTION;
		cells[c].edges[i] = EMPTY_CONNECTION;
	}
	cells[c].dead = false;
	cells[c].father = -1;
	cells[c].x_at_start = 0;
	cells[c].num_vertices = 0;
	cells[c].cell_cycle_state = 0;
	cells[c].can_divide = !autonomous_cell_cycle;
	cells[c].centroid_x = 0;
	cells[c].centroid_y = 0;
	cells[c].num_divisions = 0;
	this->num_cells++;
	//if (!dead_cells.empty()) dead_cells.pop();
	return c;
}
int Tissue::newEdge()
{
	int e;
	if (false)
	{ //if(!dead_edges.empty()){ //if(false){//
		e = dead_edges.front();
	}
	else
	{
		Edge ed;
		ed.ind = static_cast<int>(this->edges.size());
		e = ed.ind;
		edges.push_back(ed);
	}
	edges[e].vertices[0] = EMPTY_CONNECTION;
	edges[e].vertices[1] = EMPTY_CONNECTION;
	edges[e].cells[0] = EMPTY_CONNECTION;
	edges[e].cells[1] = EMPTY_CONNECTION;
	edges[e].can_transition = true;
	edges[e].dead = false;
	this->num_edges++;
	//if (!dead_edges.empty()) dead_edges.pop();
	return e;
}
///////////////////////////////////////////////////////////////////////////////////////////////
//public functions only useful for testing
int Tissue::addVertex(Vertex v)
{
	vertices.push_back(v);
	num_vertices++;
	return num_vertices - 1;
}
int Tissue::addCell(Cell c)
{
	cells.push_back(c);
	num_cells++;
	return num_cells - 1;
}
int Tissue::addEdge(Edge e)
{
	edges.push_back(e);
	num_edges++;
	return num_edges;
}

void Tissue::addAcceptedMovements(int add)
{
	this->max_accepted_movements += add;
}

void Tissue::setAcceptedMovements(int mv)
{
	this->max_accepted_movements = mv;
}

void Tissue::setStepMode(bool mode, int steps)
{
	this->write_every_N_moves = steps;
	this->step_mode = mode;
}

void Tissue::setHingeMinAndMaxPositions()
{
	calculateCellCentroid(cells[0]);
	double min, max, miny, maxy; //initial value must be that of the first hinge cell, not that of cell 0
	bool set = false;
	for (Cell c : cells)
	{
		if (c.type == CellType::hinge || c.type == CellType::vein_hinge)
		{
			c.perimeter = calculateCellPerimeter(c);
			c.area = calculateCellArea(c);
			calculateCellCentroid(c);
			//cout <<"in hingeMinMax, centroid calculated: " << c.centroid_x << ", " << c.centroid_y << endl;
			if (set)
			{
				if (c.centroid_x < min)
				{
					min = c.centroid_x;
					//cout << "cent_x menor: " <<c.centroid_x << endl;
				}
				else if (c.centroid_x > max)
				{
					max = c.centroid_x;
					//cout << "cent_x mayor: " <<c.centroid_x << endl;
				}
				if (c.centroid_y < miny)
				{
					miny = c.centroid_y;
				}
				else if (c.centroid_y > maxy)
				{
					maxy = c.centroid_y;
				}
			}
			else
			{
				c.perimeter = calculateCellPerimeter(c);
				c.area = calculateCellArea(c);
				calculateCellCentroid(c);
				min = c.centroid_x;
				max = c.centroid_x;
				miny = c.centroid_y;
				maxy = c.centroid_y;
				set = true;
			}
		}
	}
	hinge_min_xpos = min;
	hinge_max_xpos = max;
	hinge_min_ypos = miny;
	hinge_max_ypos = maxy;
	cout << "Min hinge position: " << hinge_min_xpos << "; Max hinge position: " << hinge_max_xpos << endl;
}//end set hinge max and min positions

/* void Tissue::setMinAndMaxPositions()
{
	calculateCellCentroid(cells[0]);
	double min, max, miny, maxy;
	bool set = false;
	for (Cell c : cells)
	{
		calculateCellCentroid(c);
		if (set)
		{
			if (c.centroid_x < min)
			{
				min = c.centroid_x;
			}
			else if (c.centroid_x > max)
			{
				max = c.centroid_x;
			}
			if (c.centroid_y < miny)
			{
				miny = c.centroid_y;
			}
			else if (c.centroid_y > maxy)
			{
				maxy = c.centroid_y;
			}
		}
		else
		{
			min = c.centroid_x;
			max = c.centroid_x;
			miny = c.centroid_y;
			maxy = c.centroid_y;
			set = true;
		}
	}
	min_xpos = min;
	max_xpos = max;
	min_ypos = miny;
	max_ypos = maxy;
	AP_compartment_limit = miny + (maxy - miny)*posterior_comparment_region;
	cout << "Min wing x position: " << min_xpos << "; Max wing xposition: " << max_xpos << "; A-P limit: " << AP_compartment_limit<< endl;
} */
void Tissue::setMinAndMaxPositions()
{
	calculateCellCentroid(cells[0]);
	double min, max, miny, maxy;
	bool set = false;
	for (Vertex v : vertices)
	{
		if(v.cells[0] == EMPTY_CONNECTION)
			continue;
		if (set)
		{
			if (v.x < min)
			{
				min = v.x;
			}
			else if (v.x > max)
			{
				max = v.x;
			}
			if (v.y < miny)
			{
				miny = v.y;
			}
			else if (v.y > maxy)
			{
				maxy = v.y;
			}
		}
		else
		{
			min = v.x;
			max = v.x;
			miny = v.y;
			maxy = v.y;
			set = true;
		}
	}
	min_xpos = min;
	max_xpos = max;
	min_ypos = miny;
	max_ypos = maxy;
	AP_compartment_limit = miny + (maxy - miny)*posterior_comparment_region;
	cout << "Min wing x position: " << min_xpos << "; Max wing xposition: " << max_xpos << "; A-P limit: " << AP_compartment_limit<< endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//Area of a polygon using Shoelace formula
double Tissue::calculateCellArea(const Cell &c)
{
	/* auto start = chrono::high_resolution_clock::now();  */
	double area = 0.0;
	int previous = c.num_vertices - 1;
	for (int i = 0; i < c.num_vertices; i++)
	{
		area += (vertices[c.vertices[previous]].x * vertices[c.vertices[i]].y - vertices[c.vertices[previous]].y * vertices[c.vertices[i]].x);
		previous = i;
	}
/* 	auto stop = chrono::high_resolution_clock::now();
	auto duration1 = chrono::duration_cast<chrono::nanoseconds>(stop - start);
	cout << duration1.count() << "\t"; */
	return abs(0.5 * area);
}

// Calculates centroid, (some of the terms are the same as in area, but did not join both functions since area will be called much more often
//From n to 0 because vertices are ordered clock-wise and some terms become negative
void Tissue::calculateCellCentroid(Cell &c)
{
	double xc = 0, yc = 0, a;
	int previous = 0;
	for (int i = c.num_vertices - 1; i >= 0; i--)
	{
		//cout << "centroid: num_vertices "<< c.num_vertices << ", v_i: "<< vertices[c.vertices[i]] <<endl;
		a = vertices[c.vertices[previous]].x * vertices[c.vertices[i]].y - vertices[c.vertices[previous]].y * vertices[c.vertices[i]].x;
		xc += a * (vertices[c.vertices[previous]].x + vertices[c.vertices[i]].x);
		yc += a * (vertices[c.vertices[previous]].y + vertices[c.vertices[i]].y);
		previous = i;
	}
	c.centroid_x = xc / (6 * c.area);
	c.centroid_y = yc / (6 * c.area);
	//cout << "cell centroid inside: " << c.centroid_x << ", " << c.centroid_y << ",  "<< c.area << endl;
	return;
}

inline double Tissue::calculateCellPerimeter(const Cell &c)
{
	double perim = 0.0;
	for (int i = 0; i < c.num_vertices; i++)
		perim += edges[c.edges[i]].length;
	return perim;
}

inline double Tissue::distance(int v1, int v2)
{
	return sqrt(pow(this->vertices[v1].x - this->vertices[v2].x, 2) + pow(this->vertices[v1].y - this->vertices[v2].y, 2));
}
inline double Tissue::calculateEnergy2(Vertex &v)
{
	double term1 = 0, term2 = 0, term3 = 0;
	double pref_area;
	int aux_ind;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.cells[i] != EMPTY_CONNECTION){
			aux_ind = v.cells[i];
			term1 += pow(cells[aux_ind].area / cells[aux_ind].preferred_area - 1, 2);
			term3 += 0.5 * cells[aux_ind].perimeter_contractility * cells[aux_ind].perimeter * cells[aux_ind].perimeter/ cells[aux_ind].preferred_area;
		}
	}
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.edges[i] != EMPTY_CONNECTION)
		{
			aux_ind = v.edges[i];
			pref_area = edges[aux_ind].type == EdgeType::tissue_boundary ?
								edges[aux_ind].cells[0] == EMPTY_CONNECTION ?
											cells[edges[aux_ind].cells[1]].preferred_area :
											cells[edges[aux_ind].cells[0]].preferred_area :
								(cells[edges[aux_ind].cells[0]].preferred_area + cells[edges[aux_ind].cells[1]].preferred_area) * 0.5;
			term2 += edges[aux_ind].tension * edges[aux_ind].length / sqrt(pref_area);
		}
	}
	if (v.spring != EMPTY_CONNECTION){
		term2 += springs[v.spring].tension * springs[v.spring].length / sqrt(pref_area);
	}
	return 0.5 * term1 * energy_term1 + term2 * energy_term2 + term3 * energy_term3;
}
//
inline double Tissue::calculateEnergy(Vertex &v) {
  double term1 = 0, term2 = 0, term3 = 0, momentum_term = 0;
  double aux;
  int aux_ind;
  bool is_cuticle = false;
  for (int i = 0; i < CELLS_PER_VERTEX; i++) {
    if (v.cells[i] != EMPTY_CONNECTION) {
      aux_ind = v.cells[i];
      aux = cells[aux_ind].area - cells[aux_ind].preferred_area;
      // term1 += aux*aux*K.val[static_cast<int>(cells[aux_ind].type)];
      term1 += aux * aux * cells[aux_ind].K;
      term3 += 0.5 * cells[aux_ind].perimeter_contractility *
               cells[aux_ind].perimeter * cells[aux_ind].perimeter;
    }
  }
  for (int i = 0; i < CELLS_PER_VERTEX; i++) {
    if (v.edges[i] != EMPTY_CONNECTION) {
      aux_ind = v.edges[i];
      if (edges[aux_ind].type == EdgeType::stringedge) {
		is_cuticle = true;
        if (string_equilibrium_distance) {
          // cout << "Eq. Distance" << endl;
          term2 +=
              edges[aux_ind].tension *
              pow(edges[aux_ind].length - edges[aux_ind].optimal_length, 2);

        } else {
          term2 += edges[aux_ind].tension * edges[aux_ind].length;
        }
      }
    }
  }
  if (v.spring != EMPTY_CONNECTION) {
    term2 += springs[v.spring].tension * springs[v.spring].length;
  }

  if((is_cuticle && momentum_term_cuticle > 0.0)){
	  momentum_term = momentum_term_cuticle*sqrt(pow(v.x - v.x_drag, 2) + pow(v.y - v.y_drag, 2));
  }else if( momentum_term_tissue > 0.0){
	  //Also applies to cuticle if no momentum in cuticle but there is momentum in tissue
	  momentum_term = momentum_term_tissue*sqrt(pow(v.x - v.x_drag, 2) + pow(v.y - v.y_drag, 2));
  }
  return term1 * energy_term1 + term2 * energy_term2 + term3 * energy_term3 + momentum_term;
} // calcEnergy2
inline double Tissue::calculateEnergy_term4(Vertex& v){
	int n = 0, aux_nei, aux_ind, aux_cell;
	float products[CELLS_PER_VERTEX];
	double term4 = 0;
	cell_type_param *term4coef = v.y > AP_compartment_limit ? &energy_term4_anterior: &energy_term4_posterior;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.edges[i] != EMPTY_CONNECTION)
		{
			aux_ind = v.edges[i];
			if(edges[aux_ind].type == EdgeType::tissue_boundary){
				aux_nei = edges[aux_ind].vertices[0] == v.ind ? edges[aux_ind].vertices[1] : edges[aux_ind].vertices[0];
				aux_cell = edges[aux_ind].cells[0] == EMPTY_CONNECTION ? edges[aux_ind].cells[1] : edges[aux_ind].cells[0];
				products[n] = ((vertices[aux_nei].x - v.x)/edges[aux_ind].length)*((vertices[aux_nei].x -bufferMovement.x)/bufferMovement.edge_lengths[i]);
				products[n] += ((vertices[aux_nei].y - v.y)/edges[aux_ind].length)*((vertices[aux_nei].y -bufferMovement.y)/bufferMovement.edge_lengths[i]);
				products[n] = term4coef->val[static_cast<int>(cells[aux_cell].type)]*(2-abs(products[n]));
				products[n] = abs(products[n]);
				n++;
			}
		}
	}
	if(n > 0){
		term4 = products[0];
		for(int i = 1; i < n; n++){
				term4 += products[i];
				//if(products[i]>term4)
				//	term4 = products[i];
		}
		term4 /= n;
		//term4 = term4coef->val[static_cast<int>(cells[aux_cell].type)]*(2 - term4);
	}
	return term4;
}
inline double Tissue::calculateEnergy_term4_cuticle(Vertex& v){
	int n = 0, aux_nei, aux_ind;
	float products[CELLS_PER_VERTEX];
	double term4 = 0;
	float term4coef = energy_term4_anterior.val[static_cast<int>(CellType::blade)]; //a default value
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.edges[i] != EMPTY_CONNECTION)
		{
			aux_ind = v.edges[i];
			if(edges[aux_ind].type == EdgeType::stringedge){
				aux_nei = edges[aux_ind].vertices[0] == v.ind ? edges[aux_ind].vertices[1] : edges[aux_ind].vertices[0];
				products[n] = ((vertices[aux_nei].x - v.x)/edges[aux_ind].length)*((vertices[aux_nei].x -bufferMovement.x)/bufferMovement.edge_lengths[i]);
				products[n] += ((vertices[aux_nei].y - v.y)/edges[aux_ind].length)*((vertices[aux_nei].y -bufferMovement.y)/bufferMovement.edge_lengths[i]);
				products[n] = term4coef*(2-abs(products[n]));
				products[n] = abs(products[n]);
				n++;
			}
		}
	}
	if(n > 0){
		term4 = products[0];
		for(int i = 1; i < n; n++){
				term4 += products[i];
				//if(products[i]>term4)
				//	term4 = products[i];
		}
		term4 /= n;
		//term4 = term4coef->val[static_cast<int>(cells[aux_cell].type)]*(2 - term4);
	}
	return term4;
}
void Tissue::moveVertexBack(Vertex &v)
{
	v.x = bufferMovement.x;
	v.y = bufferMovement.y;
	v.energy = bufferMovement.energy;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.edges[i] != EMPTY_CONNECTION)
		{
			edges[v.edges[i]].length = bufferMovement.edge_lengths[i];
			edges[v.edges[i]].tension = bufferMovement.edge_tensions[i];
		}
	}
	if (v.spring != EMPTY_CONNECTION){
		springs[v.spring].length = bufferMovement.spring_length;
	}
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{ // re-calculate cell areas
		if (v.cells[i] != EMPTY_CONNECTION)
		{
			cells[v.cells[i]].area = bufferMovement.cell_areas[i];
			cells[v.cells[i]].perimeter = bufferMovement.cell_perimeters[i];
		}
	}
}//moveVertexBack
void Tissue::moveVertex(Vertex &v, float x, float y)
{
	v.x = x;
	v.y = y;
	//float tension_sum = 0;  //REMOVE LATER
	//int num_edges = 0;  //REMOVE LATER
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.edges[i] != EMPTY_CONNECTION)
		{
			bufferMovement.edge_lengths[i] = edges[v.edges[i]].length;
			edges[v.edges[i]].length = distance(edges[v.edges[i]].vertices[0], edges[v.edges[i]].vertices[1]);
			if (UPDATE_EDGE_TENSION_EVERY_MOVE)
			{
				//setEdgeType(v.edges[i]); not needed since there is base_tension
				bufferMovement.edge_tensions[i] = edges[v.edges[i]].tension;
				setEdgeTension(v.edges[i]);
				//tension_sum+=edges[v.edges[i]].tension;  //REMOVE LATER
				//num_edges++;  //REMOVE LATER
				//cout <<edges[v.edges[i]].tension<<endl;
			}
		}
	}
	 /* if(NORMALIZE_EDGE_TENSION && UPDATE_EDGE_TENSION_EVERY_MOVE && num_edges==3){  //REMOVE LATER
		for (int i = 0; i < CELLS_PER_VERTEX; i++)  //REMOVE LATER
		{
			if (v.edges[i] != EMPTY_CONNECTION)  //REMOVE LATER
			{
				edges[v.edges[i]].tension /= tension_sum;  //REMOVE LATER
			}
		}
	}  */

	if (v.spring != EMPTY_CONNECTION)
	{
		bufferMovement.spring_length = springs[v.spring].length;
		springs[v.spring].length = distance(springs[v.spring].vertices[0], springs[v.spring].vertices[1]);
	}
/*  	double area, perimeter;
	int previous;
	int c; */
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{ // re-calculate cell areas
		if (v.cells[i] != EMPTY_CONNECTION)
		{
			bufferMovement.cell_areas[i] = cells[v.cells[i]].area;
			bufferMovement.cell_perimeters[i] = cells[v.cells[i]].perimeter;
			cells[v.cells[i]].area = calculateCellArea(this->cells[v.cells[i]]);
			cells[v.cells[i]].perimeter = calculateCellPerimeter(this->cells[v.cells[i]]);
/* 			c = v.cells[i]; //calc area and perimeter outside of function
 			area = 0.0;
			perimeter= 0.0;
			previous = cells[c].num_vertices - 1;
			for (int i = 0; i < cells[c].num_vertices; i++)
			{
				area += (vertices[cells[c].vertices[previous]].x * vertices[cells[c].vertices[i]].y - vertices[cells[c].vertices[previous]].y * vertices[cells[c].vertices[i]].x);
				previous = i;
				perimeter += edges[cells[c].edges[i]].length;
			}
			cells[c].area = area;
			cells[c].perimeter = perimeter;  */
		}
	}
}
//counter_favorable_accepted, counter_favorable_rejected, counter_unfav_accepted, counter_unfav_rejected
bool Tissue::tryMoveVertex()
{
	int vertex_to_move;
	do
	{
		vertex_to_move = std::rand() % static_cast<int>(vertices.size());
	} while (vertices[vertex_to_move].dead || !vertices[vertex_to_move].movable); //(!vertices[vertex_to_move].movable && STATIC_PRESENT || vertices[vertex_to_move].cells[0] == EMPTY_CONNECTION)); // A or static dead vertex cannot be selected
	bufferMovement.x = vertices[vertex_to_move].x;
	bufferMovement.y = vertices[vertex_to_move].y;
	bufferMovement.energy = calculateEnergy(vertices[vertex_to_move]); //Calculate again because it is not updated every time a cell area changes etc
	double angle, radius, new_x, new_y;
	double move_prob;

	angle = 2 * M_PI * unif(generator);
	radius = unif(generator) * max_range_vertex_movement;
	new_x = vertices[vertex_to_move].movable_x ? bufferMovement.x + cos(angle) * radius : bufferMovement.x;
	new_y = vertices[vertex_to_move].movable_y ? bufferMovement.y + sin(angle) * radius : bufferMovement.y;
	moveVertex(vertices[vertex_to_move], new_x, new_y);
	vertices[vertex_to_move].energy = calculateEnergy(vertices[vertex_to_move]);
	//REMOVE LATER
	/*bool in_cuticle = false;
	 for(int i = 0; i < CELLS_PER_VERTEX; i++){
		if(vertices[vertex_to_move].edges[i] != EMPTY_CONNECTION){
			if(edges[vertices[vertex_to_move].edges[i]].type == EdgeType::stringedge){
				in_cuticle = true;
				break;
			}

		}
	} */

	if(use_term4) /// UNTIL HERE
		vertices[vertex_to_move].energy += calculateEnergy_term4(vertices[vertex_to_move]);
	/* else if(in_cuticle)
		vertices[vertex_to_move].energy += calculateEnergy_term4_cuticle(vertices[vertex_to_move]); */
	//Prob of accepting unfavourable movement ins constant
	//move_prob = vertices[vertex_to_move].energy <= bufferMovement.energy ? temperature_negative_energy : temperature_positive_energy;
	if (temperature_means_proportion_of_acceptance)
		{ //Prob of accepting unfavourable movement ins constant
			move_prob = vertices[vertex_to_move].energy <= bufferMovement.energy ? temperature_negative_energy : temperature_positive_energy;
		}
		else
		{ //Prob of accepting unfavourable movement depends on how much unfavourable the movement is
			move_prob = vertices[vertex_to_move].energy <= bufferMovement.energy ? 
				temperature_negative_energy : 
				exp(-(vertices[vertex_to_move].energy - bufferMovement.energy) / temperature_positive_energy);
			//cout << "e_old:e_new:diff\t"<< bufferMovement.energy << "\t" << vertices[vertex_to_move].energy << "\t" << vertices[vertex_to_move].energy - bufferMovement.energy << "\n";
		}
	/* if(counter_move_trials == 2){
		cout << ">idextra\tid\tnew_x\tnew_y\tx\ty\txdif\tydif\te\tnew_e\tedif\tmove_prob\ttime\t" << 
		"old_elen0\told_elen1\told_elen2\told_eten0\told_eten1\told_eten2\t"<<
		"new_elen0\tnew_elen1\tnew_elen2\tnew_eten0\tnew_eten1\tnew_eten2\t"<<
		"old_cellarea0\told_cellarea1\told_cellarea2\told_cellper0\told_cellper1\told_cellper2\t"<<
		"new_cellarea0\tnew_cellarea1\tnew_cellarea2\tnew_cellper0\tnew_cellper1\tnew_cellper2\n";
	} */
 	 
	/* if(vertex_to_move == 440 || vertex_to_move == 441){
		cout << ">idd\t" <<
		vertex_to_move <<
	"\t" << vertices[vertex_to_move].x <<
	"\t" << vertices[vertex_to_move].y <<
	"\t" <<  bufferMovement.x <<
	"\t" <<  bufferMovement.y << 
	"\t" <<  vertices[vertex_to_move].x - bufferMovement.x <<
	"\t" <<  vertices[vertex_to_move].y - bufferMovement.y <<
	"\t" << bufferMovement.energy << 
	"\t" << vertices[vertex_to_move].energy << 
	"\t" << vertices[vertex_to_move].energy - bufferMovement.energy << 
	"\t" << move_prob <<
	"\t" << counter_move_trials <<
	"\t" << bufferMovement.edge_lengths[0] <<
	"\t" << bufferMovement.edge_lengths[1] <<
	"\t" << bufferMovement.edge_lengths[2] <<
	"\t" << bufferMovement.edge_tensions[0] <<
	"\t" << bufferMovement.edge_tensions[1] <<
	"\t" << bufferMovement.edge_tensions[2] <<
	"\t" << edges[vertices[vertex_to_move].edges[0]].length << 
	"\t" << edges[vertices[vertex_to_move].edges[1]].length << 
	"\t" << edges[vertices[vertex_to_move].edges[2]].length << 
	"\t" << edges[vertices[vertex_to_move].edges[0]].tension << 
	"\t" << edges[vertices[vertex_to_move].edges[1]].tension << 
	"\t" << edges[vertices[vertex_to_move].edges[2]].tension << 
	"\t" << bufferMovement.cell_areas[0] <<
	"\t" << bufferMovement.cell_areas[1] <<
	"\t" << bufferMovement.cell_areas[2] <<
	"\t" << bufferMovement.cell_perimeters[0] <<
	"\t" << bufferMovement.cell_perimeters[1] <<
	"\t" << bufferMovement.cell_perimeters[2] <<
	"\t" << cells[vertices[vertex_to_move].cells[0]].area << 
	"\t" << cells[vertices[vertex_to_move].cells[1]].area << 
	"\t" << cells[vertices[vertex_to_move].cells[2]].area << 
	"\t" << cells[vertices[vertex_to_move].cells[0]].perimeter << 
	"\t" << cells[vertices[vertex_to_move].cells[1]].perimeter << 
	"\t" << cells[vertices[vertex_to_move].cells[2]].perimeter << 
	"\n";  
	} */
	double move = unif(generator);
	if (move < move_prob) //|| vertices[vertex_to_move].moves_rejected/vertices[vertex_to_move].moves_accepted > 5
	{
/* 		for(int i = 0; i < CELLS_PER_VERTEX; i++){
			if(vertices[vertex_to_move].edges[i] != EMPTY_CONNECTION){
				setEdgeTension(vertices[vertex_to_move].edges[i]);
			}
		} */
		vertices[vertex_to_move].moves_accepted ++;
		vertices[vertex_to_move].x_drag = (1 - momentum_ponderate_current)*vertices[vertex_to_move].x_drag + momentum_ponderate_current*vertices[vertex_to_move].x;
		vertices[vertex_to_move].y_drag = (1 - momentum_ponderate_current)*vertices[vertex_to_move].y_drag + momentum_ponderate_current*vertices[vertex_to_move].y;
		

		detectChangesAfterMove(vertex_to_move);
		if (autonomous_cell_cycle)
		{ //NOTE: Maybe make this before detecting changes so the effect is immediate? Otherwise cells have to wait until next round
			advanceCellCycle(vertex_to_move);
		}
		/* if(time_controls_size && gradient_with_same_final_area){
			advanceSizeWithCoordsAndTime_v2(vertex_to_move);
		}
		else if (time_controls_size && !(xcoord_controls_size || ycoord_controls_size))
		{
			advanceSizeWithTime(vertex_to_move);
		}
		else if (time_controls_size && (xcoord_controls_size || ycoord_controls_size))
		{
			advanceSizeWithCoordsAndTime(vertex_to_move);
		}
		else if (xcoord_controls_size || ycoord_controls_size)
		{
			advanceSizeWithCoords(vertex_to_move);
		} */

		/* if (time_controls_perim && !(xcoord_controls_perim || ycoord_controls_perim))
		{
			advancePerimWithTime(vertex_to_move);
		}
		else if (time_controls_perim && (xcoord_controls_perim || ycoord_controls_perim))
		{
			advancePerimWithCoordsAndTime(vertex_to_move);
		}
		else if (xcoord_controls_perim || ycoord_controls_perim)
		{
			advancePerimWithCoords(vertex_to_move);
		} */
		/* if(K_gradient_x || K_gradient_y){
			advanceKWithCoords(vertex_to_move);
		} */
		//Counters
		if(time_controls_size || time_controls_perim){
			advanceSizeAndPerimWithTime_new(vertex_to_move);
		}
		if (vertices[vertex_to_move].energy <= bufferMovement.energy)
		{
			counter_favorable_accepted++;
		}
		else
		{
			counter_unfav_accepted++;
		} //Counters
/* 		for(int vv; vv < CELLS_PER_VERTEX; vv++){
			if(vertices[vertex_to_move].cells[vv] == 7998 || vertices[vertex_to_move].cells[vv] == 8397){
				checkCellConsistency(vertices[vertex_to_move].cells[vv]);
			}
		} */
		return true;
	}
	else
	{
		//Counters
		vertices[vertex_to_move].moves_rejected ++;
		if (vertices[vertex_to_move].energy <= bufferMovement.energy)
		{
			counter_favorable_rejected++;
		}
		else
		{
			counter_unfav_rejected++;
		} //Counters
		moveVertexBack(vertices[vertex_to_move]);
		return false;
	}
}

/*double Tissue::calculateTerm4Energy(Vertex &v, double old_x, double old_y)
{
	double ycent = 0.5*(max_ypos - min_ypos);
	double yval = (1 - abs(ycent - old_y)/ycent)*(1 - difference_flow_rate) + difference_flow_rate;
	//cout << "ymin " << min_ypos << ", ymax " << max_ypos <<"ycent " << ycent << ", yval " << yval << ", oldx " << old_x << ", vx " << v.x << ", e: " << (old_x - v.x)*energy_term4*yval << endl;
	return (old_x - v.x)*energy_term4*yval;
}*/
inline double Tissue::expAdvance(double x, float exponent)
{
	return exponent > 0 ?
	(1 - exp(-pow(x, exponent))) / EXP_FACTOR :
	((1 - exp(-pow(x, exponent))) - EXP_FACTOR) / (1 - EXP_FACTOR);
}

void Tissue::advanceSizeWithCoordsAndTime(int vertex_moved)
{
	int caux;
	//double auxprint;
	double time_factor = static_cast<double>(counter_moves_accepted) / upper_bound_movements;
	time_factor = time_factor >= 1? 1.0 : expAdvance(time_factor, time_decrease_exponent);
	//time_factor = expAdvance(time_factor, time_decrease_exponent); //(1 - exp(-pow(time_factor, time_decrease_exponent)))/EXP_FACTOR;
	float pos_factor = 0, pos_factor_y = 0, aux_area, ini, fin;

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		cells[caux].preferred_area = preferred_area_initial.val[static_cast<int>(cells[caux].type)] + (preferred_area_final.val[static_cast<int>(cells[caux].type)] - preferred_area_initial.val[static_cast<int>(cells[caux].type)]) * time_factor; //Caution: problems if step_mode is on
		if (cells[caux].type == CellType::vein_blade || cells[caux].type == CellType::blade)
		{
			if (keep_area_after_division)
			{
				cells[caux].preferred_area /= pow(2, cells[caux].num_divisions);
			}
			continue;
		} //auxprint = cells[caux].preferred_area;
		calculateCellCentroid(cells[caux]);
		if (xcoord_controls_size)
		{
			pos_factor = cells[caux].centroid_x - hinge_min_xpos;
			pos_factor = pos_factor < 0 ? 0 : pos_factor / (hinge_max_xpos - hinge_min_xpos);
		}
		if (ycoord_controls_size)
		{
			pos_factor_y = cells[caux].centroid_y - hinge_min_ypos;
			pos_factor_y = pos_factor_y < 0 ? 0 : pos_factor_y / (hinge_max_ypos - hinge_min_ypos);
			pos_factor = xcoord_controls_size ? 0.5 * (pos_factor_y + pos_factor) : pos_factor_y;
		}
		if(use_blade_area_for_coord_gradient){
			if (cells[caux].type == CellType::vein_hinge)
			{
				ini = preferred_area_initial.val[static_cast<int>(CellType::vein_blade)];
				fin = preferred_area_final.val[static_cast<int>(CellType::vein_blade)];
			}
			else //if (cells[caux].type == CellType::hinge) //unnecessary
			{
				ini = preferred_area_initial.val[static_cast<int>(CellType::blade)];
				fin = preferred_area_final.val[static_cast<int>(CellType::blade)];
			}
			aux_area = ini + (fin - ini) * time_factor;
		}else{
			aux_area = preferred_area_initial_gradient.val[static_cast<int>(cells[caux].type)];
		}
		cells[caux].preferred_area += (aux_area - cells[caux].preferred_area) * expAdvance(pos_factor, xcoord_decrease_exponent);
		if (keep_area_after_division)
		{
			cells[caux].preferred_area /= pow(2, cells[caux].num_divisions);
		}
		//cout << endl;
	}
} //advanceSizeWithTimeAndXcoord

void Tissue::advanceSizeWithCoordsAndTime_v2(int vertex_moved)
{
	int caux;
	//double auxprint;
	double time_factor;// = static_cast<double>(counter_moves_accepted) / upper_bound_movements;
	//time_factor = time_factor >= 1? 1.0 : expAdvance(time_factor, time_decrease_exponent);
	//time_factor = expAdvance(time_factor, time_decrease_exponent); //(1 - exp(-pow(time_factor, time_decrease_exponent)))/EXP_FACTOR;
	float pos_factor = 0, pos_factor_y = 0, aux_area, ini, fin;

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		if (cells[caux].type == CellType::vein_blade || cells[caux].type == CellType::blade)
		{
			time_factor = counter_moves_accepted >= upper_bound_movements ? 1 : static_cast<double>(counter_moves_accepted) / upper_bound_movements;
			cells[caux].preferred_area = preferred_area_initial.val[static_cast<int>(cells[caux].type)] + (preferred_area_final.val[static_cast<int>(cells[caux].type)] - preferred_area_initial.val[static_cast<int>(cells[caux].type)]) * time_factor; //Caution: problems if step_mode is on
			if (keep_area_after_division)
			{
				cells[caux].preferred_area /= pow(2, cells[caux].num_divisions);
			}
			continue;
		} //auxprint = cells[caux].preferred_area;
		calculateCellCentroid(cells[caux]);
		/* if (xcoord_controls_size)
		{
			pos_factor = cells[caux].centroid_x - hinge_min_xpos;
			pos_factor = pos_factor < 0 ? 0 : pos_factor / (hinge_max_xpos - hinge_min_xpos);
		}
		if (ycoord_controls_size)
		{
			pos_factor_y = cells[caux].centroid_y - hinge_min_ypos;
			pos_factor_y = pos_factor_y < 0 ? 0 : pos_factor_y / (hinge_max_ypos - hinge_min_ypos);
			pos_factor = xcoord_controls_size ? 0.5 * (pos_factor_y + pos_factor) : pos_factor_y;
		} */
		pos_factor = cells[caux].centroid_x - hinge_min_xpos;
		pos_factor = pos_factor < 0 ? 1 : pos_factor / (hinge_max_xpos - hinge_min_xpos);
		ini = preferred_area_initial.val[static_cast<int>(cells[caux].type)];
		fin = preferred_area_final.val[static_cast<int>(cells[caux].type)];	
		time_factor = static_cast<double>(counter_moves_accepted) / (upper_bound_movements*pos_factor);
		time_factor = time_factor > 1 ? 1 : time_factor;				
		
		cells[caux].preferred_area =  ini + (fin - ini) * time_factor;
		if (keep_area_after_division)
		{
			cells[caux].preferred_area /= pow(2, cells[caux].num_divisions);
		}
		//cout << endl;
	}
} //advancePerimWithTimeAndXcoord_v2

void Tissue::advancePerimWithCoordsAndTime(int vertex_moved)
{
	int caux;
	//double auxprint;
	double time_factor = static_cast<double>(counter_moves_accepted) / upper_bound_movements;
	time_factor = time_factor >= 1? 1.0 : expAdvance(time_factor, time_decrease_exponent);
	//time_factor = expAdvance(time_factor, time_decrease_exponent); //(1 - exp(-pow(time_factor, time_decrease_exponent)))/EXP_FACTOR;
	float pos_factor = 0, pos_factor_y = 0, aux_perim, ini, fin;

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		cells[caux].perimeter_contractility = perimeter_contract.val[static_cast<int>(cells[caux].type)] + (perimeter_contract_final.val[static_cast<int>(cells[caux].type)] - perimeter_contract.val[static_cast<int>(cells[caux].type)]) * time_factor; //Caution: problems if step_mode is on
		if (cells[caux].type == CellType::vein_blade || cells[caux].type == CellType::blade)
		{
			continue;
		} //auxprint = cells[caux].preferred_area;
		calculateCellCentroid(cells[caux]);
		if (xcoord_controls_perim)
		{
			pos_factor = cells[caux].centroid_x - hinge_min_xpos;
			pos_factor = pos_factor < 0 ? 0 : pos_factor / (hinge_max_xpos - hinge_min_xpos);
		}
		if (ycoord_controls_perim)
		{
			pos_factor_y = cells[caux].centroid_y - hinge_min_ypos;
			pos_factor_y = pos_factor_y < 0 ? 0 : pos_factor_y / (hinge_max_ypos - hinge_min_ypos);
			pos_factor = xcoord_controls_size ? 0.5 * (pos_factor_y + pos_factor) : pos_factor_y;
		}
		if (cells[caux].type == CellType::vein_hinge)
		{
			ini = perimeter_contract.val[static_cast<int>(CellType::vein_blade)];
			fin = perimeter_contract_final.val[static_cast<int>(CellType::vein_blade)];
		}else{
			ini = perimeter_contract.val[static_cast<int>(CellType::blade)];
			fin = perimeter_contract_final.val[static_cast<int>(CellType::blade)];
		}
		aux_perim = ini + (fin - ini) * time_factor;
		cells[caux].perimeter_contractility += (aux_perim - cells[caux].perimeter_contractility) * expAdvance(pos_factor, xcoord_decrease_exponent);
	}
} //advancePerimWithTimeAndXcoord


void Tissue::advanceSizeWithCoords(int vertex_moved)
{
	int caux;
	float pos_factor = 0, pos_factor_y = 0, ini, fin;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		if (cells[caux].type == CellType::vein_blade || cells[caux].type == CellType::blade){
			continue;
		}
		if(use_blade_area_for_coord_gradient){
		 if (cells[caux].type == CellType::hinge)
			{
				ini = preferred_area_initial.val[static_cast<int>(CellType::hinge)];
				fin = preferred_area_initial.val[static_cast<int>(CellType::blade)];
			}
			else //if (cells[caux].type == CellType::vein_hinge) //not necessary
			{
				ini = preferred_area_initial.val[static_cast<int>(CellType::vein_hinge)];
				fin = preferred_area_initial.val[static_cast<int>(CellType::vein_blade)];
				//cout << "HingeVein cell: " << " x abs = " << cells[caux].centroid_x << ", x rel = " << pos_factor << ", pref area = " << cells[caux].preferred_area << endl;
			}
		}else{
			 if (cells[caux].type == CellType::hinge)
			{
				ini = preferred_area_initial.val[static_cast<int>(CellType::hinge)];
			}
			else //if (cells[caux].type == CellType::vein_hinge) //not necessary
			{
				ini = preferred_area_initial.val[static_cast<int>(CellType::vein_hinge)];
				//cout << "HingeVein cell: " << " x abs = " << cells[caux].centroid_x << ", x rel = " << pos_factor << ", pref area = " << cells[caux].preferred_area << endl;
			}
			fin = preferred_area_initial_gradient.val[static_cast<int>(cells[caux].type)];
		}
		calculateCellCentroid(cells[caux]);
		if (xcoord_controls_size)
		{
			pos_factor = cells[caux].centroid_x - hinge_min_xpos;
			pos_factor = pos_factor < 0 ? 0 : pos_factor / (hinge_max_xpos - hinge_min_xpos);
		}
		if (ycoord_controls_size)
		{
			pos_factor_y = cells[caux].centroid_y - hinge_min_ypos;
			pos_factor_y = pos_factor_y < 0 ? 0 : pos_factor_y / (hinge_max_ypos - hinge_min_ypos);
			pos_factor = xcoord_controls_size ? 0.5 * (pos_factor_y + pos_factor) : pos_factor_y;
		}
		cells[caux].preferred_area = ini + (fin - ini) * expAdvance(pos_factor, xcoord_decrease_exponent); //(1 - exp(-pow(pos_factor, xcoord_decrease_exponent)))/EXP_FACTOR;
		if (keep_area_after_division)
		{
			cells[caux].preferred_area /= pow(2, cells[caux].num_divisions);
		}
	}
} //advanceSizeWithXcoord

void Tissue::advancePerimWithCoords(int vertex_moved)
{
	int caux;
	float pos_factor = 0, pos_factor_y = 0, ini, fin;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		if (cells[caux].type == CellType::vein_blade || cells[caux].type == CellType::blade){
			continue;
		}
		if (cells[caux].type == CellType::hinge)
		{
			ini = perimeter_contract.val[static_cast<int>(CellType::hinge)];
			fin = perimeter_contract.val[static_cast<int>(CellType::blade)];
		}
		else //if (cells[caux].type == CellType::vein_hinge) //not necessary
		{
			ini = perimeter_contract.val[static_cast<int>(CellType::vein_hinge)];
			fin = perimeter_contract.val[static_cast<int>(CellType::vein_blade)];
				//cout << "HingeVein cell: " << " x abs = " << cells[caux].centroid_x << ", x rel = " << pos_factor << ", pref area = " << cells[caux].preferred_area << endl;
		}
		
		calculateCellCentroid(cells[caux]);
		if (xcoord_controls_size)
		{
			pos_factor = cells[caux].centroid_x - hinge_min_xpos;
			pos_factor = pos_factor < 0 ? 0 : pos_factor / (hinge_max_xpos - hinge_min_xpos);
		}
		if (ycoord_controls_size)
		{
			pos_factor_y = cells[caux].centroid_y - hinge_min_ypos;
			pos_factor_y = pos_factor_y < 0 ? 0 : pos_factor_y / (hinge_max_ypos - hinge_min_ypos);
			pos_factor = xcoord_controls_size ? 0.5 * (pos_factor_y + pos_factor) : pos_factor_y;
		}
		cells[caux].perimeter_contractility = ini + (fin - ini) * expAdvance(pos_factor, xcoord_decrease_exponent); //(1 - exp(-pow(pos_factor, xcoord_decrease_exponent)))/EXP_FACTOR;
	}
} //advancePerimWithXcoord

void Tissue::setKWithCoords(int cell){
	int celltype = static_cast<int>(cells[cell].type);
	float ini = K.val[celltype];
	float fin = K_final.val[celltype];
	float pos_factor = 0, pos_factor_y = 0;
	if(!(K_gradient_x || K_gradient_y)){
		cells[cell].K = ini;
	}else{
	//calculateCellCentroid(cells[caux]); //Assumes it has been calculated
		if (K_gradient_x)
		{
			if(cells[cell].type == CellType::hinge || cells[cell].type == CellType::vein_hinge){
				pos_factor = cells[cell].centroid_x - hinge_min_xpos;
				pos_factor = pos_factor < 0 ? 0 : pos_factor / (hinge_max_xpos - hinge_min_xpos);
			}else{
				pos_factor = cells[cell].centroid_x - hinge_max_xpos; // assumes that max hinge == min blade
				pos_factor = pos_factor < 0 ? 0 : pos_factor / (max_xpos - hinge_max_xpos);
			}
		}
		if (K_gradient_y){
			pos_factor_y = cells[cell].centroid_y - min_ypos;
			pos_factor_y = pos_factor_y < 0 ? 0 : pos_factor_y / (max_ypos - min_ypos);
			pos_factor = K_gradient_x ? 0.5 * (pos_factor_y + pos_factor) : pos_factor_y;
		}
		cells[cell].K = ini + (fin - ini) * pos_factor;//expAdvance(pos_factor, K_grad_exponent);
	}
}
void Tissue::advanceKWithCoords(int vertex_moved)
{
	int caux;
	float pos_factor = 0, pos_factor_y = 0, ini, fin;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		/* if (cells[caux].type == CellType::vein_blade || cells[caux].type == CellType::blade){
			continue;
		} */
		int celltype = static_cast<int>(cells[caux].type);
		ini = K.val[celltype];
		fin = K_final.val[celltype];
		calculateCellCentroid(cells[caux]);
		if (K_gradient_x)
		{
			pos_factor = cells[caux].centroid_x - hinge_min_xpos;
			pos_factor = pos_factor < 0 ? 0 : pos_factor / (hinge_max_xpos - hinge_min_xpos);
		}
		if (K_gradient_y)
		{
			pos_factor_y = cells[caux].centroid_y - hinge_min_ypos;
			pos_factor_y = pos_factor_y < 0 ? 0 : pos_factor_y / (hinge_max_ypos - hinge_min_ypos);
			pos_factor = K_gradient_x ? 0.5 * (pos_factor_y + pos_factor) : pos_factor_y;
		}
		cells[caux].K = ini + (fin - ini) * pos_factor;//expAdvance(pos_factor, K_grad_exponent); //(1 - exp(-pow(pos_factor, xcoord_decrease_exponent)))/EXP_FACTOR;
	}
} //advanceKWithXcoord
void Tissue::advanceCellCycle(int vertex_moved)
{
	int caux;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		cells[caux].cell_cycle_state += 1; //cell cycle state is float
		if (cells[caux].cell_cycle_state >= cells[caux].cell_cycle_limit)
		{
			cells[caux].can_divide = true;
		}
		if (cell_cycle_controls_size && 
		((cells[caux].type == CellType::blade || cells[caux].type == CellType::vein_blade)))
		{
			cells[caux].preferred_area = preferred_area_initial.val[static_cast<int>(cells[caux].type)] + 
			(preferred_area_final.val[static_cast<int>(cells[caux].type)] - 
			preferred_area_initial.val[static_cast<int>(cells[caux].type)]) * 
			cells[caux].cell_cycle_state / cells[caux].cell_cycle_limit;
			if(keep_area_after_division)
				cells[caux].preferred_area /= pow(2, cells[caux].num_divisions);
		}
	}
} //advanceCellCycle

void Tissue::advanceSizeAndPerimWithTime_new(int vertex_moved)
{
	int caux;
	double time_factor = static_cast<double>(counter_moves_accepted) / upper_bound_movements;
	float initial;
	time_factor = time_factor >= 1? 1.0 : expAdvance(time_factor, time_decrease_exponent);
	//time_factor = expAdvance(time_factor, time_decrease_exponent); //(1 - exp(-pow(time_factor, time_decrease_exponent)))/EXP_FACTOR;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		if(time_controls_size){
			initial = preferred_area_initial.val[static_cast<int>(cells[caux].type)];
			cells[caux].preferred_area = initial + (cells[caux].base_preferred_area - initial) * time_factor; //Caution: problems if step_mode is on
			//cout << "cell="<<caux <<"; type="<<static_cast<int>(cells[caux].type) << "; time factor= " << time_factor << ", init=" << preferred_area_initial[cells[caux].type] << "; final=" << preferred_area_final[cells[caux].type] << "; area=" << cells[caux].preferred_area << endl;
			if (keep_area_after_division)
			{
				cells[caux].preferred_area /= pow(2, cells[caux].num_divisions);
			}
		}
		if(time_controls_perim){
			initial = perimeter_contract.val[static_cast<int>(cells[caux].type)];
			cells[caux].perimeter_contractility = initial + (cells[caux].base_perimeter - initial) * time_factor; //Caution: problems if step_mode is on			
		}

	}
} //advanceSizeAndPerimWithTime_new

void Tissue::advanceSizeWithTime(int vertex_moved)
{
	int caux;
	double time_factor = static_cast<double>(counter_moves_accepted) / upper_bound_movements;
	time_factor = time_factor >= 1? 1.0 : expAdvance(time_factor, time_decrease_exponent);
	//time_factor = expAdvance(time_factor, time_decrease_exponent); //(1 - exp(-pow(time_factor, time_decrease_exponent)))/EXP_FACTOR;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		cells[caux].preferred_area = preferred_area_initial.val[static_cast<int>(cells[caux].type)] + (preferred_area_final.val[static_cast<int>(cells[caux].type)] - preferred_area_initial.val[static_cast<int>(cells[caux].type)]) * time_factor; //Caution: problems if step_mode is on
		//cout << "cell="<<caux <<"; type="<<static_cast<int>(cells[caux].type) << "; time factor= " << time_factor << ", init=" << preferred_area_initial[cells[caux].type] << "; final=" << preferred_area_final[cells[caux].type] << "; area=" << cells[caux].preferred_area << endl;
		if (keep_area_after_division)
		{
			cells[caux].preferred_area /= pow(2, cells[caux].num_divisions);
		}
	}
} //advanceSizeWithTime

void Tissue::advancePerimWithTime(int vertex_moved)
{
	int caux;
	double time_factor = static_cast<double>(counter_moves_accepted) / upper_bound_movements;
	time_factor = time_factor >= 1? 1.0 : expAdvance(time_factor, time_decrease_exponent);
	//time_factor = expAdvance(time_factor, time_decrease_exponent); //(1 - exp(-pow(time_factor, time_decrease_exponent)))/EXP_FACTOR;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		cells[caux].perimeter_contractility = perimeter_contract.val[static_cast<int>(cells[caux].type)] + (perimeter_contract_final.val[static_cast<int>(cells[caux].type)] - perimeter_contract.val[static_cast<int>(cells[caux].type)]) * time_factor; //Caution: problems if step_mode is on
		//cout << "cell="<<caux <<"; type="<<static_cast<int>(cells[caux].type) << "; time factor= " << time_factor << ", init=" << preferred_area_initial[cells[caux].type] << "; final=" << preferred_area_final[cells[caux].type] << "; area=" << cells[caux].preferred_area << endl;
	}
} //advancePerimWithTime

void Tissue::calculateBasePrefAreaAndPerim(Cell& cell){
	float pos_factor_x = (cell.centroid_x - min_xpos)/(max_xpos - min_xpos);
	float pos_factor_y = 0, ini_perim, ini_area, fin_perim, fin_area, factor_area, factor_perim;
	if((cell.type == CellType::blade || cell.type == CellType::vein_blade) && 
			(reference_for_gradient != USE_PROPORTION_OF_WING || pos_factor_x >= wing_proportion_in_gradient) &&
			(reference_for_gradient != USE_BARE_EXPONENTIAL_GRADIENT) &&
			(reference_for_gradient != GRADIENT_FROM_CENTER_ALL) &&
			(reference_for_gradient != GRADIENT_FROM_CENTER_HINGE)){
		cell.base_perimeter = perimeter_contract_final.val[static_cast<int>(cell.type)];
		cell.base_preferred_area = preferred_area_final.val[static_cast<int>(cell.type)];
		cell.perimeter_contractility = perimeter_contract.val[static_cast<int>(cell.type)];		
		cell.preferred_area = preferred_area_initial.val[static_cast<int>(cell.type)];
		return;
	}

	if(reference_for_gradient == USE_BARE_EXPONENTIAL_GRADIENT ||
			(reference_for_gradient == GRADIENT_FROM_CENTER_ALL) ||
			(reference_for_gradient == GRADIENT_FROM_CENTER_HINGE)){
		pos_factor_x = exp(- xcoord_decrease_exponent * pos_factor_x);
		cell.base_perimeter = perimeter_contract_final.val[static_cast<int>(cell.type)];
		cell.base_preferred_area = preferred_area_final.val[static_cast<int>(cell.type)];
		if(xcoord_controls_perim)
			cell.base_perimeter *= pos_factor_x;
		if(xcoord_controls_size)
			cell.base_preferred_area *= (1 - pos_factor_x);

		//If there is gradient from medial to lateral
		if(reference_for_gradient == GRADIENT_FROM_CENTER_ALL ||
		(reference_for_gradient == GRADIENT_FROM_CENTER_HINGE && 
		(cell.type == CellType::hinge || cell.type == CellType::vein_hinge))
		){		
			float dist_from_center = 1 - wing_proportion_in_gradient*abs(cell.centroid_y - AP_compartment_limit)/abs(AP_compartment_limit - (cell.centroid_y < AP_compartment_limit ? min_ypos : max_ypos));
			//float CLperim_value = perimeter_contract.val[static_cast<int>(cell.type)] + dist_from_center*(perimeter_contract_final.val[static_cast<int>(cell.type)] - perimeter_contract.val[static_cast<int>(cell.type)]);
			//cell.base_perimeter = cell.base_perimeter*(1 - wing_proportion_in_gradient) + CLperim_value*wing_proportion_in_gradient;
			//cell.base_perimeter += CLperim_value*wing_proportion_in_gradient;
			cell.base_perimeter = cell.base_perimeter*dist_from_center + perimeter_contract.val[static_cast<int>(cell.type)];
		}
		//Time dependence
		cell.perimeter_contractility = time_controls_perim ? 
				perimeter_contract.val[static_cast<int>(CellType::blade)] : 
				cell.base_perimeter;
		cell.preferred_area = time_controls_size ? 
				preferred_area_initial.val[static_cast<int>(CellType::blade)] : 
				cell.base_preferred_area;
		return;
	}
	
	if(reference_for_gradient == USE_HINGE_BLADE_FRONTIER){
		std::vector<int> edges_frontier = getHingeBladeFrontier();
		pos_factor_x = getXgradFromFrontier(cell.centroid_x, cell.centroid_y, edges_frontier);
	}else if(reference_for_gradient == USE_MAX_HINGE_POSITION){
		pos_factor_x = cell.centroid_x - hinge_min_xpos;
		pos_factor_x = pos_factor_x < 0 ? 0 : pos_factor_x / (hinge_max_xpos - hinge_min_xpos);
	}else{
		pos_factor_x = pos_factor_x > wing_proportion_in_gradient ? 1: pos_factor_x/wing_proportion_in_gradient;
	}
	pos_factor_y = cell.centroid_y - hinge_min_ypos;
	pos_factor_y = pos_factor_y < 0 ? 0 : pos_factor_y / (hinge_max_ypos - hinge_min_ypos);

	if(xcoord_controls_perim & ycoord_controls_perim)
		factor_perim = 0.5*(pos_factor_x + pos_factor_y);
	else if(xcoord_controls_perim)
		factor_perim = pos_factor_x;
	else if(ycoord_controls_perim)
		factor_perim = pos_factor_y;
	else
		factor_perim = 0;

	if(xcoord_controls_size & ycoord_controls_size)
		factor_area = 0.5*(pos_factor_x + pos_factor_y);
	else if(xcoord_controls_size)
		factor_area = pos_factor_x;
	else if(ycoord_controls_size)
		factor_area = pos_factor_y;
	else
		factor_area = 0;		

	if (cell.type == CellType::hinge) //Ini and Fin mean proximo-distal (or antero-posterior)
	{
		ini_perim = perimeter_contract_final.val[static_cast<int>(CellType::hinge)];
		fin_perim = perimeter_contract_final.val[static_cast<int>(CellType::blade)];
		ini_area = preferred_area_final.val[static_cast<int>(CellType::hinge)];
		fin_area = preferred_area_final.val[static_cast<int>(CellType::blade)];
		cell.perimeter_contractility = time_controls_perim ? perimeter_contract.val[static_cast<int>(CellType::blade)] : ini_perim;
		cell.preferred_area = time_controls_size ? preferred_area_initial.val[static_cast<int>(CellType::blade)] : ini_area;
	}
	else //if (cells[caux].type == CellType::vein_hinge) //not necessary
	{
		ini_perim = perimeter_contract_final.val[static_cast<int>(CellType::vein_hinge)];
		fin_perim = perimeter_contract_final.val[static_cast<int>(CellType::vein_blade)];
		ini_area = preferred_area_final.val[static_cast<int>(CellType::vein_hinge)];
		fin_area = preferred_area_final.val[static_cast<int>(CellType::vein_blade)];
		cell.perimeter_contractility = time_controls_perim ? perimeter_contract.val[static_cast<int>(CellType::vein_blade)] : ini_perim;
		cell.preferred_area = time_controls_size ? preferred_area_initial.val[static_cast<int>(CellType::vein_blade)] : ini_area;
	}

	cell.base_perimeter = ini_perim + (fin_perim - ini_perim) * expAdvance(factor_perim, xcoord_decrease_exponent); //(1 - exp(-pow(pos_factor, xcoord_decrease_exponent)))/EXP_FACTOR;
	cell.base_preferred_area = ini_area + (fin_area - ini_area) * expAdvance(factor_area, xcoord_decrease_exponent); //(1 - exp(-pow(pos_factor, xcoord_decrease_exponent)))/EXP_FACTOR;
	cell.perimeter_contractility = time_controls_perim ? ini_perim : cell.base_perimeter;
	cell.preferred_area = time_controls_size ? ini_area : cell.base_preferred_area;
	
}

std::vector<int> Tissue::getHingeBladeFrontier(){
	std::vector<int> edge_frointier;
	CellType c1, c2;
	for(Edge e : edges){
		if(e.cells[0] != EMPTY_CONNECTION && e.cells[1] != EMPTY_CONNECTION){
			c1 = cells[e.cells[0]].type;
			c2 = cells[e.cells[1]].type;
			if((c1 == CellType::blade && c2 == CellType::hinge) || (c2 == CellType::blade && c1 == CellType::hinge) ||
				(c1 == CellType::vein_blade && c2 == CellType::vein_hinge) || (c2 == CellType::vein_blade && c1 == CellType::vein_hinge) ||
				(c1 == CellType::vein_blade && c2 == CellType::hinge) || (c2 == CellType::vein_blade && c1 == CellType::hinge) ||
				(c1 == CellType::blade && c2 == CellType::vein_hinge) || (c2 == CellType::blade && c1 == CellType::vein_hinge)
			){
				edge_frointier.push_back(e.ind);
			}
		}
	}
	return edge_frointier;
}


double Tissue::getXgradFromFrontier(float xpos, float ypos, std::vector<int> & frontier){
	int closest_edge = frontier[0];
	double min_ydist = abs(0.5*(vertices[edges[frontier[0]].vertices[0]].y + vertices[edges[frontier[0]].vertices[1]].y) - ypos);
	double dist;
	for(int e: frontier){
		dist = abs(0.5*(vertices[edges[e].vertices[0]].y + vertices[edges[e].vertices[1]].y) - ypos);
		if(dist < min_ydist){
			min_ydist = dist;
			closest_edge = e;
		}			
	}
	double x_end = 0.5*(vertices[edges[closest_edge].vertices[0]].x + vertices[edges[closest_edge].vertices[1]].x);
	return (xpos - min_xpos)/(x_end - min_xpos);
}




void Tissue::produceOutputs(std::string add_to_name)
{

	std:string fname = simname + "_" + add_to_name + "_" + std::to_string(written_files);
	if(RECALCULATE_CENTROIDS_FOR_PRINTING){
		for(int i = 0; i < cells.size(); i++)
			calculateCellCentroid(cells[i]);
	}
	writeCellsFile(fname);
	writePointsFile(fname);
	if (num_springs > 0)
		writeSpringsFile(fname);
	//writeAllData(fname); //THIS MIGHT USEFUL TO DEBUG, BUT THE FORMAT IS NOT READ BY PLOTTING PROGRAM
	if(WRITE_DATA_TABLES){
		writeEdgeDataTable(fname);
		writeCellDataTable(fname);
		writeSpringDataTable(fname);
		writePointsDataTable(fname);
	}
	if (REPORT_OUT > 0)
		cout << "\nWritting file: " << written_files << " at move " << counter_moves_accepted << endl;
		std::cout << getStats() << endl;
	written_files++;
}

void Tissue::derivativeVertexPos(const Vertex &v, pointDerivative &pd)
{
	double t1x, t2x, t3x, t1y, t2y, t3y;
	t1x = t2x = t3x = t1y = t2y = t3y = 0;
	float pref_area, aux, auxd1, auxd2;
	int aux_this, aux_prev, aux_next, c;
	//TERM 1
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.cells[i] == EMPTY_CONNECTION){
			//cout << "No cell " << i << endl;
			continue;
		}
		c = v.cells[i];
        aux = cells[c].K * (cells[c].area - cells[c].preferred_area);//(cells[c].area/cells[c].preferred_area- 1)/cells[c].preferred_area;
		aux_this = which(v.ind, cells[c].vertices, cells[c].num_vertices); //These three lines are repeated in TERM3 but since TERM3 is probably going to be optimized it doesnt matter
		//IMPORTANT: ASSUMES THAT VERTICES IN CELLS ARE ORDERED CLOCK-WISE
		//Calculate derivative of area according to Shoelace formula from Wikipedia
		aux_next = cells[c].vertices[(aux_this + cells[c].num_vertices - 1)%cells[c].num_vertices ];
		aux_prev = cells[c].vertices[(aux_this + 1)%cells[c].num_vertices ];
		t1x += aux*(vertices[aux_prev].y - vertices[aux_next].y);
		t1y += aux*(vertices[aux_next].x - vertices[aux_prev].x);
		//cout << "Term 1: aux " << aux << ", pref area " << cells[c].preferred_area << ", area " << cells[c].area << ", cell " << c << endl;

	} //Term1
	//TERM 2
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.edges[i] == EMPTY_CONNECTION){
			continue;
		}
		c = v.edges[i];
		aux_next = edges[c].vertices[0] == v.ind ? edges[c].vertices[1] : edges[c].vertices[0];
		//aux =  0.5 * edges[c].tension / (edges[c].length > NUMERIC_THRESHOLD ? edges[c].length : NUMERIC_THRESHOLD);
		aux = 0.5 * edges[c].tension / edges[c].length;
		t2x += (v.x - vertices[aux_next].x) * aux;
		t2y += (v.y - vertices[aux_next].y) * aux;
		//cout << "Term 2: aux " << aux << ", x - next: " << v.x - vertices[aux_next].x << ", y - next: " << v.y - vertices[aux_next].y << ", len: " << edges[c].length << ", tens: " << edges[c].tension  << ", edge " << c  << endl;
	} //Term 2
	if (v.spring != EMPTY_CONNECTION){
		 aux = 0.5 * springs[v.spring].tension / (springs[v.spring].length > NUMERIC_THRESHOLD ? springs[v.spring].length : NUMERIC_THRESHOLD);
		aux_next = springs[v.spring].vertices[0] == v.ind ? springs[v.spring].vertices[1] : springs[v.spring].vertices[0];
		t2x += (v.x - vertices[aux_next].x) * aux; // Here may be a numeric problem when distance between 2 points is already 0
		t2y += (v.y - vertices[aux_next].y) * aux;
		//cout << "Term 2 spring: aux" << aux << ", x - next: " << v.x - vertices[aux_next].x << ", y - next: " << v.y - vertices[aux_next].y << ", len: " << springs[v.spring].length << ", tens: " << springs[v.spring].tension  << endl;
	}//Spring part of term2
	//Term 3 THIS TERM MUST BE OPTIMIZED
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.cells[i] == EMPTY_CONNECTION){
			//cout << "No cell " << i << endl;
			continue;
		}
		c = v.cells[i];
		aux = 0.5 * cells[c].perimeter * cells[c].perimeter_contractility;
		aux_this = which(v.ind, cells[c].vertices, cells[c].num_vertices);
		aux_prev = cells[c].vertices[(aux_this + cells[c].num_vertices - 1)%cells[c].num_vertices ];
		aux_next = cells[c].vertices[(aux_this + 1)%cells[c].num_vertices ];
        auxd1 = distance(v.ind, aux_prev);
        auxd2 = distance(v.ind, aux_next); //Maybe numeric problem
		auxd1 = auxd1 > NUMERIC_THRESHOLD ? auxd1 : NUMERIC_THRESHOLD;
		auxd2 = auxd2 > NUMERIC_THRESHOLD ? auxd2 : NUMERIC_THRESHOLD;
		t3x += aux*( (v.x - vertices[aux_prev].x)/auxd1 + (v.x - vertices[aux_next].x)/auxd2);
		t3y += aux*( (v.y - vertices[aux_prev].y)/auxd1 + (v.y - vertices[aux_next].y)/auxd2);
		//cout << "Term 3 aux " << aux << ", d1: " << auxd1 << ",d2 " << auxd2 << endl;
	} //Term3
	//cout << "v.ind=" << v.ind << ": x1=" << t1x << ", y1=" << t1y << ", x2=" << t2x << ", y2=" << t2y << ", x3=" << t3x << ", y3=" << t3y;

	/*t1x = isnan(t1x) || isinf(t1x) ? 0 : t1x;
	t1y = isnan(t1y) || isinf(t1y) ? 0 : t1y;
	t2x = isnan(t2x) || isinf(t2x) ? 0 : t2x;
	t2y = isnan(t2y) || isinf(t2y) ? 0 : t2y;
	t3x = isnan(t3x) || isinf(t3x) ? 0 : t3x;
	t3y = isnan(t3y) || isinf(t3y) ? 0 : t3y;
*/
	pd.x = 0.5 * t1x * energy_term1 + t2x * energy_term2 + t3x * energy_term3;
	pd.y = 0.5 * t1y * energy_term1 + t2y * energy_term2 + t3y * energy_term3;
	if (max_range_vertex_movement > 0 && (abs(pd.x) > max_range_vertex_movement || abs(pd.y) > max_range_vertex_movement))
	{
		float auxpdx = pd.x;
		pd.x = pd.x / (pd.x + pd.y) * max_range_vertex_movement;
		pd.y = pd.y / (auxpdx + pd.y) * max_range_vertex_movement;
	}
	//Add some noise
	if(temperature_positive_energy > 0){
		double angle = 2 * M_PI * unif(generator);
		double radius = unif(generator) * temperature_positive_energy;
		pd.x += cos(angle) * radius;
		pd.y += sin(angle) * radius;
	}
	//cout << ", dx=" << pd.x << ", dy=" << pd.y << endl << endl;
	cout << "v=" << v.ind << ", dx=" << pd.x << ", dy=" << pd.y << ", t1x=" << t1x << ", t1y=" << t1y << ", t2x=" << t2x << ", t2y=" << t2y << ", t3x=" << t3x << ", t3y=" << t3y << endl;
}//derivativeVertexPos used in Simulate Euler
void Tissue::simulateMonteCarlo()
{
	cout << "SIMULATING with MONTE CARLO method for " << max_accepted_movements << " accepted movements." << endl;
	if (!step_mode)
	{
		produceOutputs(); // Initial print of the state of the
	}
	//MAIN SIMULATION LOOPSimulate
	do
	{
		if (tryMoveVertex())
		{ //Tries a random movement; if accepted:
			counter_moves_accepted++;
			performRearrangements(); // Performs any transition (rearrangement) needed
			if ((counter_moves_accepted % write_every_N_moves == 0 && !step_mode)) 
				produceOutputs(); //Print files that then can be used to plot
		}//End if move accepted
		counter_move_trials++;
	} while (counter_moves_accepted < max_accepted_movements); 
}//simulateMonteCarlo

void Tissue::simulateEuler()
{
	pointDerivative_v derxy(num_vertices, pointDerivative{0.0, 0.0});
	//float total_time = static_cast<float>(max_accepted_movements);
	int h_factor = round(1.0/h);
	int write_every_N_moves2 = write_every_N_moves*h_factor;
	int max_accepted_movements2 = max_accepted_movements*h_factor;
	cout << "SIMULATING with EULER method." << endl;
	cout << "h: " << h << ", h factor: "<< h_factor << ", write every: " << write_every_N_moves2 << ", max: " << max_accepted_movements2 <<endl;
	if (!step_mode)
		produceOutputs();
	for (int t = 1; t <= max_accepted_movements2; t ++)
	{
		//cout << "*****" << endl << endl << t << ", t%h_factor: " << t % h_factor << ", t%write_every_N_moves2: " << t % write_every_N_moves2 << endl;
		//Calculate derivative of x and y positions
		for (int v = 0; v < num_vertices; v++)
		{
			if (vertices[v].dead || !vertices[v].movable)
				continue;
			derivativeVertexPos(vertices[v], derxy[v]);
		}

		//Now update vertices, edges and cells
		for (int v = 0; v < num_vertices; v++)
		{
			if (vertices[v].dead || !vertices[v].movable)
				continue;
			//cout << "     v: " << v << ", der_x: " << h*derxy[v].x << ", der_y: " << h*derxy[v].y << endl;
			if(isnan(derxy[v].x) || isnan(derxy[v].y)){
				cout << "Exiting because of NaN found: " << endl;
				cout << "     v: " << v << ", der_x: " << h*derxy[v].x << ", der_y: " << h*derxy[v].y << endl;
				exit(1);
			}
			vertices[v].x += h*derxy[v].x;
			vertices[v].y += h*derxy[v].y;
		}//end update vertices
		//Update edges
		for (int e = 0; e < num_edges; e++)
		{
			if (edges[e].dead)
				continue;
			edges[e].length = distance(edges[e].vertices[0], edges[e].vertices[1]);
		}//end update edges
		//update springs
		for (int e = 0; e < num_springs; e++)
		{
			if (springs[e].dead)
				continue;
			springs[e].length = distance(springs[e].vertices[0], springs[e].vertices[1]);
		}//end update springs
		//update cells
		for (int c = 0; c < num_cells; c++)
		{
			if (cells[c].dead)
				continue;
			cells[c].perimeter = calculateCellPerimeter(cells[c]);
			cells[c].area = calculateCellArea(cells[c]);
			calculateCellCentroid(cells[c]);
		}//end update cells

		for (int v = 0; v < num_vertices; v++)
			detectChangesAfterMove(v);
		performRearrangements();
		while (vertices.size() > derxy.size())
			derxy.push_back(pointDerivative{0.0, 0.0});
		//Now update parameters that  change with time or space
		//produce outputs
		if(t % h_factor == 0){
			counter_moves_accepted++;
			if (t % write_every_N_moves2 == 0 && !step_mode)
				produceOutputs();
		}
	} //end for time
	//produceOutputs();
} // end simulateEuler

void Tissue::simulate(std::default_random_engine &generator, std::uniform_real_distribution<double> &unif)
{
	srand(random_seed);
	this->generator = generator;
	this->unif = unif;

	switch (this->integration_mode)
	{
	case INTEGR_MONTECARLO:
		simulateMonteCarlo();
		break;
	case INTEGR_EULER:
		simulateEuler(); //Doesn't work
		break;
	case INTEGR_MIXTURE:
		break;
	case INTEGR_RUNGEKUTTA4:
		break;
	default:
		break;
	}
}

void Tissue::detectChangesAfterMove(int vertex_moved)
{

	//Check if T1 transitions are needed
	Edge *ee;
	int caux;
	Vertex *v1;
	Vertex *v2;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (vertices[vertex_moved].edges[i] == EMPTY_CONNECTION)
			continue;
		ee = &edges[vertices[vertex_moved].edges[i]];
		//if(ee->type == EdgeType::stringedge)//Redundant: can_transition = False in these edges (strings/cuticle)
		//	continue;
		v1 = &vertices[ee->vertices[0]];
		v2 = &vertices[ee->vertices[1]];
		if (ee->length <= t1_transition_critical_distance && v1->movable && v2->movable && ee->can_transition)
		{ //if length < critical distance
			if (ee->cells[0] != EMPTY_CONNECTION && ee->cells[1] != EMPTY_CONNECTION)
			{ //if the edge touches two cells, and there are 4 cells involved, then T1
				int cellnum = 0;
				for (int j = 0; j < CELLS_PER_VERTEX; j++){
					if (v1->cells[j] != EMPTY_CONNECTION)
						cellnum++;
					if (v2->cells[j] != EMPTY_CONNECTION && !contains(v2->cells[j], v1->cells, CELLS_PER_VERTEX))
						cellnum++;
				}
				if (cellnum == 4 && cells[ee->cells[0]].num_vertices > 3 && cells[ee->cells[1]].num_vertices > 3)
				{ //empty_connections == 0
					rearrangements_needed.push(Rearrangement{ee->ind, RearrangementType::t1});
					//If less than 4 cells involved, t1 transition outwards //empty_connections == 1
				}
				else if (cellnum == 3 && cells[ee->cells[0]].num_vertices > 3 && cells[ee->cells[1]].num_vertices > 3)
				{
					rearrangements_needed.push(Rearrangement{vertices[vertex_moved].edges[i], RearrangementType::t1_at_border_outwards});
				}
				else if (cellnum == 2 && cells[ee->cells[0]].num_vertices > 3 && cells[ee->cells[1]].num_vertices > 3)
				{
					rearrangements_needed.push(Rearrangement{vertices[vertex_moved].edges[i], RearrangementType::join_limit_edges});
				}
			}
			else
			{ //Otherwise, length of union of edges touching two vertices touching the edge
				int edgenum = 0;
				for (int j = 0; j < CELLS_PER_VERTEX; j++){
					if (v2->edges[j] != EMPTY_CONNECTION)
						edgenum++;
					if (v1->edges[j] != EMPTY_CONNECTION && !contains(v1->edges[j], v2->edges, CELLS_PER_VERTEX))
						edgenum++;
				}
				if (edgenum == 5 && v1->spring == EMPTY_CONNECTION && v2->spring == EMPTY_CONNECTION)
					rearrangements_needed.push(Rearrangement{vertices[vertex_moved].edges[i], RearrangementType::t1_at_border_inwards});
				else if (edgenum < 5)
					rearrangements_needed.push(Rearrangement{vertices[vertex_moved].edges[i], RearrangementType::join_limit_edges});
			}
		}
	}

	//Check if T2 transitions (cell death) or cell divisions are needed
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		if (cells[caux].area <= t2_transition_critical_area && cells[caux].num_vertices == 3)
		{
			rearrangements_needed.push(Rearrangement{caux, RearrangementType::t2});
		}
		else if (cells[caux].area >= max_cell_area.val[static_cast<int>(cells[caux].type)] && cells[caux].can_divide)
		{
			rearrangements_needed.push(Rearrangement{caux, RearrangementType::divide_cell});
		}
	}
	if(active_t1_prob > 0){
		if(unif(generator) < active_t1_prob)
			addRandomTransition();
	}

} //detect rearrangements
void Tissue::addRandomTransition(){
	int edge;
	do{
		edge = std::rand() % num_edges;
	}while(edges[edge].dead && static_cast<int>(edges[edge].type) < 2 && vertices[edges[edge].vertices[0]].movable && vertices[edges[edge].vertices[1]].movable);
	float angle = atan2(vertices[edges[edge].vertices[0]].x - vertices[edges[edge].vertices[1]].x, vertices[edges[edge].vertices[0]].y - vertices[edges[edge].vertices[1]].y);
	float ansin = abs(sin(angle));
	if(ansin < maxsin2rant1 && ansin > minsin2rant1){
		rearrangements_needed.push(Rearrangement{edge, RearrangementType::random_t1});
	}
	//cout << "Random T1 in edge " << edge << endl;
}
void Tissue::performRearrangements()
{

	Rearrangement r;
	while (!rearrangements_needed.empty())
	{
		r = rearrangements_needed.front();
		rearrangements_needed.pop();
		switch (r.type)
		{
		case RearrangementType::random_t1:
		case RearrangementType::t1:
			//cout << "Entering t1 \n";
			if (t1_active)
				make_t1(r);
			break;
		case RearrangementType::t2:
			//cout << "Entering t2 \n";
			if (t2_active)
				make_t2(r);
			break;
		case RearrangementType::divide_cell:
			//cout << "Entering division \n";
			if (division_active)
				make_divide_cell(r);
			break;
		case RearrangementType::divide_edge:
			break;
		case RearrangementType::join_limit_edges:
			//cout << "Entering join edges\n";
			if (join_edges_active)
				make_join_limit_edges(r);
			break;
		case RearrangementType::t1_at_border_outwards:
			//cout << "Entering t1 outwards\n";
			if (t1_outwards_active)
				make_t1_at_border_outwards(r);
			break;
		case RearrangementType::t1_at_border_inwards:
			//cout << "Entering t1 inwards " << counter_move_trials << "\n";
			if (t1_inwards_active)
				make_t1_at_border_inwards(r);
			break;
		default:
			break;
		}
	}
}

void Tissue::make_t1_at_border_inwards(Rearrangement &r)
{
	int edge = r.element_index;
	if (edges[edge].length > T1_TRANSITION_CRITICAL_DISTANCE)
		return; //Check that condition is still true (other rearrangements could have taken place since detection)
	Vertex *v1 = &vertices[edges[edge].vertices[0]];
	Vertex *v2 = &vertices[edges[edge].vertices[1]];
	//cout << "A, v1: " << v1->ind << " v2: " << v2->ind << endl;
	if (!(edges[edge].cells[0] == EMPTY_CONNECTION || edges[edge].cells[1] == EMPTY_CONNECTION))
		return;

	double old_x1 = v1->x, old_x2 = v2->x, old_y1 = v1->y, old_y2 = v2->y, old_length = edges[edge].length;

	//Get the 2 cells that are common to both vertices and cell that is exclusive to v1
	int common_cell1 = -1, common_cell2 = -1, only_v1 = -1, only_v2 = -1;
	t1_getCellRelationships(v1, v2, common_cell1, common_cell2, only_v1, only_v2);
	if (common_cell1 == EMPTY_CONNECTION)
	{
		common_cell1 = common_cell2;
		common_cell2 = EMPTY_CONNECTION;
	}
	//cout << "B " << common_cell1 << " " << only_v1 << " " << only_v2 << "\n";
	Cell *cc1 = &this->cells[common_cell1]; //These pointers to Cell structures will be used later
	//Cell* cc2 = &this->cells[common_cell2];
	Cell *sp1 = &this->cells[only_v1];
	Cell *sp2 = &this->cells[only_v2];
	//cout << "Cn";
	if (cc1->num_vertices < 4)
		return; //One of the cells can't lose any other vertex
	//t1_rotate_edge(v1, v2, edge, cc1, cc2, sp1, sp2);
	t1_rotate_edge90degrees(v1, v2, edge);
	//cout << "D\n";
	//Determine which cell goes with which point. Compare new point locations of v1 and v2 with all the neighbours of v1 and v2.
	double dist1 = t1_inwards_get_dist_sum(v1, v2, cc1, sp1, sp2);
	double dist2 = t1_inwards_get_dist_sum(v2, v1, cc1, sp2, sp1);

	//cout << "E\n";
	if (dist2 > 0 && (dist2 < dist1 || dist1 < 0))
	{ //common_cell1 goes with v1, and v2 goes with EMPTY_CONNECTION

		v1 = &vertices[edges[edge].vertices[1]];
		v2 = &vertices[edges[edge].vertices[0]];
		int aux = only_v1;
		only_v1 = only_v2;
		only_v2 = aux;
		sp1 = &this->cells[only_v1];
		sp2 = &this->cells[only_v2];
	}
	else if (dist2 < 0 && dist1 < 0)
	{ //edges always cross
		v1->x = old_x1;
		v2->x = old_x2;
		v1->y = old_y1;
		v2->y = old_y2;
		edges[edge].length = old_length;
		counter_t1_abortions++;
		return;
	}
	ofstream ff;
	if (REPORT_T1)
	{
		ff.open("T1_REPORT_inwards_" + std::to_string(counter_moves_accepted) + "_" + std::to_string(v1->ind) + "_" + std::to_string(v2->ind) + ".txt");

		ff << "\nBEFORE\n\n";
		ff << "edge: " << edge << "\n\n";
		ff << VERTEX_HEADER << *v1 << "\n"
		   << *v2 << "\n\n";
		ff << CELL_HEADER << *cc1 << "\tcommon_cell1\n"
		   << *sp1 << "\tsp1\n"
		   << *sp2 << "\tsp2\n\n";
		ff << getStats();

		ff << "\n\nDist when cell " << cc1->ind << " goes with vertex " << v1->ind << ": " << dist1 << "\n";
		ff << "Dist when cell " << cc1->ind << " goes with vertex " << v2->ind << ": " << dist2 << "\n";
	}
	//cout << "F\n";
	//Now we know that v1 is closest to common_cell1, and v2 should be closer to common_cell2. Time to update everything
	//1) Cell exclussive of v1 is going to touch v2, and cell exclussive of v2 is going to touch v1
	v1->cells[which(EMPTY_CONNECTION, v1->cells, CELLS_PER_VERTEX)] = only_v2;
	v2->cells[which(common_cell1, v2->cells, CELLS_PER_VERTEX)] = only_v1;
	//cout << "G\n";
	//2)Remove edge from cell that was common to both vertices
	for (int i = which(edge, cc1->edges, MAX_SIDES_PER_CELL); i < cc1->num_vertices; i++)
		cc1->edges[i] = i == MAX_SIDES_PER_CELL - 1 ? EMPTY_CONNECTION : cc1->edges[i + 1];
	//3) Remove one of the vertices from each of the cells that were common to both vertices
	for (int i = which(v2->ind, cc1->vertices, cc1->num_vertices); i < cc1->num_vertices; i++)
		cc1->vertices[i] = i == MAX_SIDES_PER_CELL - 1 ? EMPTY_CONNECTION : cc1->vertices[i + 1];
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
	setEdgeType(edge);
	setEdgeTension(edge);
	//cout << "L\n";

	this->counter_t1_inwards++;

	if (REPORT_OUT > 1){
		//double cent_x = 0.5*(v1->x + v2->x);
		//double cent_y = 0.5*(v1->y + v2->y);
		//printLine("T1_inwards", v1->ind, v2->ind, cent_x, cent_y, static_cast<int>(edges[edge].type));
		printLine_t1("T1_inwards", v1->ind, v2->ind, old_x1, old_x2, old_y1, old_y2, v1->x, v2->x, v1->y, v2->y, static_cast<int>(edges[edge].type));
	
	}
	if (REPORT_T1)
	{
		ff << "\n\nAFTER\n\n";
		ff << VERTEX_HEADER << *v1 << "\n"
		   << *v2 << "\n\n";
		ff << CELL_HEADER << *cc1 << "\tcommon_cell1\n"
		   << *sp1 << "\tsp1\n"
		   << *sp2 << "\tsp2\n\n";
		ff << getStats();
		ff.close();
	}

} // end t1 at border inwards

void Tissue::make_t1_at_border_outwards(Rearrangement &r)
{

	int edge = r.element_index;
	if (edges[edge].length > t1_transition_critical_distance)
		return; //Check that condition is still true (other rearrangements could have taken place since detection)

	Vertex *v1 = &vertices[edges[edge].vertices[0]]; //INSIDE CELL
	Vertex *v2 = &vertices[edges[edge].vertices[1]]; //IN BORDER OF CELL
	//cout << "A" << endl;
	int num_cells = 0;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
		if (v1->cells[i] != EMPTY_CONNECTION)
			num_cells++;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
		if (v2->cells[i] != EMPTY_CONNECTION && !contains(v2->cells[i], v1->cells, CELLS_PER_VERTEX))
			num_cells++;
	if (num_cells != 3)
		return;

	//Check which of the two vertices is facing outwards. The one facing outwards will be v2. Assumes that all vertices inside tissue contact with three cells
	if (!contains(EMPTY_CONNECTION, v2->cells, CELLS_PER_VERTEX))
	{
		if (contains(EMPTY_CONNECTION, v1->cells, CELLS_PER_VERTEX))
		{
			v1 = &vertices[edges[edge].vertices[1]];
			v2 = &vertices[edges[edge].vertices[0]];
		}
	}
	//cout << "B" << endl;
	double old_x1 = v1->x, old_x2 = v2->x, old_y1 = v1->y, old_y2 = v2->y, old_length = edges[edge].length;

	//Get the 2 cells that are common to both vertices and cell that is exclusive to v1
	int common_cell1 = -1, common_cell2 = -1, only_v1 = -1, only_v2 = -1;		   //, only_v1_indv1 = -1, only_v2_indv2 = -1;
	t1_getCellRelationships(v1, v2, common_cell1, common_cell2, only_v1, only_v2); //now only_v2 will be EMPTY_CONNECTION
	Cell *cc1 = &this->cells[common_cell1];										   //These pointers to Cell structures will be used later
	Cell *cc2 = &this->cells[common_cell2];
	Cell *sp1 = &this->cells[only_v1];
	//Cell* sp2 = &this->cells[only_v2];
	//cout << "C" << endl;
	if (cc1->num_vertices < 4 || cc2->num_vertices < 4)
		return; //One of the cells can't lose any other vertex
	//t1_rotate_edge(v1, v2, edge, cc1, cc2, sp1, sp2);
	t1_rotate_edge90degrees(v1, v2, edge);
	//cout << "D" << endl;
	//Determine which cell goes with which point. Compare new point locations of v1 and v2 with all the neighbours of v1 and v2.

	double dist1 = t1_get_dist_sum(v1, v2, cc1, cc2, sp1, sp1);
	if (dist1 >= 0)
	{
		int aux = common_cell1;
		common_cell1 = common_cell2;
		common_cell2 = aux;
		cc1 = &this->cells[common_cell1];
		cc2 = &this->cells[common_cell2];
	}

	//Now we know that v1 is closest to common_cell1, and v2 should be closer to common_cell2. Time to update everything
	//1) Cell exclussive of v1 is going to touch v2, and cell exclussive of v2 is going to touch v1
	v1->cells[which(common_cell2, v1->cells, CELLS_PER_VERTEX)] = EMPTY_CONNECTION;
	v2->cells[which(common_cell1, v2->cells, CELLS_PER_VERTEX)] = only_v1;
	//cout << "G" << endl;
	//2)Remove edge from cells that were common to both vertices
	for (int i = which(edge, cc1->edges, MAX_SIDES_PER_CELL); i < cc1->num_vertices; i++)
		cc1->edges[i] = i == MAX_SIDES_PER_CELL - 1 ? EMPTY_CONNECTION : cc1->edges[i + 1];
	for (int i = which(edge, cc2->edges, MAX_SIDES_PER_CELL); i < cc2->num_vertices; i++)
		cc2->edges[i] = i == MAX_SIDES_PER_CELL - 1 ? EMPTY_CONNECTION : cc2->edges[i + 1];
	//cout << "H" << endl;
	//3) Remove one of the vertices from each of the cells that were common to both vertices
	for (int i = which(v2->ind, cc1->vertices, cc1->num_vertices); i < cc1->num_vertices; i++)
		cc1->vertices[i] = i == MAX_SIDES_PER_CELL - 1 ? EMPTY_CONNECTION : cc1->vertices[i + 1];
	cc1->num_vertices--;
	for (int i = which(v1->ind, cc2->vertices, cc2->num_vertices); i < cc2->num_vertices; i++)
		cc2->vertices[i] = i == MAX_SIDES_PER_CELL - 1 ? EMPTY_CONNECTION : cc2->vertices[i + 1];
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
	edges[edge].tension = line_tension_tissue_boundary.val[static_cast<int>(sp1->type)];
	edges[edge].type = EdgeType::tissue_boundary;
	this->counter_t1_outwards++;

	if (REPORT_OUT > 1){
		//double cent_x = 0.5*(v1->x + v2->x);
		//double cent_y = 0.5*(v1->y + v2->y);
		//printLine("T1_outwards", v1->ind, v2->ind, cent_x, cent_y, static_cast<int>(edges[edge].type));
		printLine_t1("T1_outwards", v1->ind, v2->ind, old_x1, old_x2, old_y1, old_y2, v1->x, v2->x, v1->y, v2->y, static_cast<int>(edges[edge].type));
	
	}

} /// End make_t1_outwards

void Tissue::make_divide_cell(Rearrangement &r)
{
	int cell = r.element_index;
	if (cells[cell].area < cells[cell].max_area || cells[cell].num_vertices < 3 || cells[cell].dead || !cells[cell].can_divide)
		return;

	double x1, x2, y1, y2;
	int e1 = EMPTY_CONNECTION, e2 = EMPTY_CONNECTION;
	if (!getDivisionPoints(cell, x1, x2, y1, y2, e1, e2))
	{
		return;
	}
	if(e1 == EMPTY_CONNECTION || e2 == EMPTY_CONNECTION)
	{
		if(REPORT_OUT > 1){
			cout << "At make_divide_cell: after getDivisionPoints edge is not found"<< endl;
		}
		return;
	}
	int newvind1 = newVertex(x1, y1);
	int newvind2 = newVertex(x2, y2); //create new vertices that are going to be positioned at (x1, y1) and (x2, y2) insideedges e1 and e2
	int newcind = newCell();
	cells[newcind].type = cells[cell].type;

	cells[newcind].division_angle_random_noise = cells[cell].division_angle_random_noise;
	cells[newcind].division_angle_longest = cells[cell].division_angle_longest;
	cells[newcind].division_angle_external = cells[cell].division_angle_external;
	cells[newcind].division_angle_external_degrees = cells[cell].division_angle_external_degrees;

	cells[newcind].vary_line_tension = cells[cell].vary_line_tension;
	cells[newcind].edge_angle_prop_external = cells[cell].edge_angle_prop_external;
	cells[newcind].edge_angle_prop_uniform = cells[cell].edge_angle_prop_uniform;
	cells[newcind].edge_angle_prop_maxangle = cells[cell].edge_angle_prop_maxangle;
	cells[newcind].edge_tension_external = cells[cell].edge_tension_external;
	cells[newcind].edge_maxangle = cells[cell].edge_maxangle;
	cells[newcind].edge_spatialmax_tension = cells[cell].edge_spatialmax_tension;
	cells[newcind].edge_spatialmin_tension = cells[cell].edge_spatialmin_tension;
	cells[newcind].K = cells[cell].K;

	int newe = newEdge();

	splitEdgeWithVertex(e1, cell, newvind1); //Changes old cell but does nothing to new (only splits edges, not cell)
	splitEdgeWithVertex(e2, cell, newvind2);

	//connect new edge
	edges[newe].cells[0] = cell;
	edges[newe].cells[1] = newcind;
	edges[newe].vertices[0] = newvind1;
	edges[newe].vertices[1] = newvind2;
	edges[newe].length = distance(newvind1, newvind2);
	//edges[newe].tension = edges[e1].type == EdgeType::tissue_boundary? edges[e2].tension : edges[e1].tension;//make this more sophisticated to take into account veins etc.
	//edges[newe].type = edges[e1].type == EdgeType::tissue_boundary? edges[e2].type : edges[e1].type; //make this more sophisticated to take into account veins etc.
	setEdgeType(newe);
	setEdgeTension(newe);

	cells[cell].edges[which(EMPTY_CONNECTION, cells[cell].edges, MAX_SIDES_PER_CELL)] = newe;
	cells[newcind].edges[which(EMPTY_CONNECTION, cells[newcind].edges, MAX_SIDES_PER_CELL)] = newe;

	vertices[newvind2].cells[which(EMPTY_CONNECTION, vertices[newvind2].cells, CELLS_PER_VERTEX)] = newcind;
	vertices[newvind1].cells[which(EMPTY_CONNECTION, vertices[newvind1].cells, CELLS_PER_VERTEX)] = newcind;
	vertices[newvind2].edges[which(EMPTY_CONNECTION, vertices[newvind2].edges, CELLS_PER_VERTEX)] = newe;
	vertices[newvind1].edges[which(EMPTY_CONNECTION, vertices[newvind1].edges, CELLS_PER_VERTEX)] = newe;
	addNeighbourVertex(newvind1, newvind2);

	int newv1_in_c1 = which(newvind1, cells[cell].vertices, cells[cell].num_vertices);
	int newv2_in_c1 = which(newvind2, cells[cell].vertices, cells[cell].num_vertices);

	//Move half of the vertices of one cell to the other, and update edges accordingly
	int i = (newv1_in_c1 + 1) % cells[cell].num_vertices;
	int changeCounter = 0;
	while (i != newv2_in_c1)
	{
		for (int j = 0; j < cells[cell].num_vertices; j++) //Move edge with vertex i to new cell
		{
			if (cells[cell].edges[j] == EMPTY_CONNECTION)
				continue; //necessary when the first edge removed is not the last
			newe = cells[cell].edges[j];
			if (edges[newe].vertices[0] == cells[cell].vertices[i] || edges[newe].vertices[1] == cells[cell].vertices[i])
			{

				cells[newcind].edges[which(EMPTY_CONNECTION, cells[newcind].edges, MAX_SIDES_PER_CELL)] = newe;
				cells[cell].edges[j] = EMPTY_CONNECTION;
				edges[newe].cells[which(cell, edges[newe].cells, 2)] = newcind;
			}
		}
		cells[newcind].vertices[cells[newcind].num_vertices] = cells[cell].vertices[i];
		cells[newcind].num_vertices++;
		vertices[cells[cell].vertices[i]].cells[which(cell, vertices[cells[cell].vertices[i]].cells, CELLS_PER_VERTEX)] = newcind;
		cells[cell].vertices[i] = EMPTY_CONNECTION;
		changeCounter++;
		i = (i + 1) % cells[cell].num_vertices; //Do not decrease cell->num_vertices !
	}
	cells[newcind].vertices[cells[newcind].num_vertices] = newvind2;
	cells[newcind].vertices[cells[newcind].num_vertices + 1] = newvind1;
	cells[newcind].num_vertices += 2;

	//Reagrupate edges and vertices in cell
	for (int j = 0; j < MAX_SIDES_PER_CELL - 1; j++)
	{
		if (cells[cell].edges[j] == EMPTY_CONNECTION) //reagrupate edges
		{
			for (int k = j + 1; k < MAX_SIDES_PER_CELL; k++)
			{
				if (cells[cell].edges[k] != EMPTY_CONNECTION)
				{
					cells[cell].edges[j] = cells[cell].edges[k];
					cells[cell].edges[k] = EMPTY_CONNECTION;
					break;
				}
			}
		} //end if reagrupate edges
		if (cells[cell].vertices[j] == EMPTY_CONNECTION) //reagrupate vertices
		{
			for (int k = j + 1; k < MAX_SIDES_PER_CELL; k++)
			{
				if (cells[cell].vertices[k] != EMPTY_CONNECTION)
				{
					cells[cell].vertices[j] = cells[cell].vertices[k];
					cells[cell].vertices[k] = EMPTY_CONNECTION;
					break;
				}
			} //end if reagrupate vertices
		}
	}
	cells[cell].num_vertices -= changeCounter;

	//Now update cell areas and perimeters
	cells[newcind].perimeter_contractility = cells[cell].perimeter_contractility;
	cells[newcind].base_perimeter = cells[cell].base_perimeter; //most proximal perim contract to calculate spatial gradient
	cells[newcind].cell_cycle_limit = cells[cell].cell_cycle_limit;

	cells[newcind].max_area = cells[cell].max_area;

	if (this->autonomous_cell_cycle)
	{
		//cells[newcind].cell_cycle_state = 0; //Already done in newCell()
		cells[cell].cell_cycle_state = 0;
		//cells[newcind].can_divide = 0; //Already done in newCell()
		cells[cell].can_divide = false;
	}

	cells[newcind].area = calculateCellArea(cells[newcind]);
	cells[newcind].perimeter = calculateCellPerimeter(cells[newcind]);
	cells[cell].area = calculateCellArea(cells[cell]);
	cells[cell].perimeter = calculateCellPerimeter(cells[cell]);

	cells[cell].num_divisions++;
	cells[newcind].num_divisions = cells[cell].num_divisions;

	if (keep_area_after_division)
	{
		cells[cell].preferred_area /= 2; //pow(2, cells[cell].num_divisions);
	}
	else if (this->cell_cycle_controls_size)
	{
		cells[newcind].preferred_area = preferred_area_initial.val[static_cast<int>(cells[cell].type)];
	}
	cells[newcind].preferred_area = cells[cell].preferred_area;
	cells[newcind].base_preferred_area = cells[cell].base_preferred_area;
	cells[newcind].father = cells[cell].father;
	cells[newcind].x_at_start = cells[cell].x_at_start;



	//Finally, check whether new vertices should be static
	/*int num_static_neighbours = 0;
	int aux_nei;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		aux_nei = vertices[newvind1].neighbour_vertices[i];
		if(aux_nei != EMPTY_CONNECTION && aux_nei != newvind2) if(!vertices[aux_nei].movable) num_static_neighbours++;
	}
	if(num_static_neighbours >= 2) vertices[newvind1].movable = false;
	num_static_neighbours = 0;
	for(int i = 0; i < CELLS_PER_VERTEX; i++){
		aux_nei = vertices[newvind2].neighbour_vertices[i];
		if(aux_nei != EMPTY_CONNECTION && aux_nei != newvind1) if(!vertices[aux_nei].movable) num_static_neighbours++;
	}
	if(num_static_neighbours >= 2) vertices[newvind2].movable = false;
	*/
	counter_divisions++;

	if (REPORT_DIV)
		writeAllData(simname + "_div_2" + to_string(counter_divisions));

	past_divisions.push(DivisionRecord{cell, newcind});
	if (REPORT_OUT > 1){
		calculateCellCentroid(cells[cell]);
		calculateCellCentroid(cells[newcind]);
		//double cent_x = 0.5 * (vertices[newvind1].x + vertices[newvind2].x);
		//double cent_y = 0.5 * (vertices[newvind1].y + vertices[newvind2].y);
		//printLine("DIVISION", cell, newcind, cent_x, cent_y, static_cast<int>(cells[cell].type));
		printLine_div("DIVISION", cell, newcind, cells[cell].centroid_x, cells[cell].centroid_y,cells[newcind].centroid_x, cells[newcind].centroid_y,static_cast<int>(cells[cell].type));

	}
	//cout << "!!!!!written: "<< written_files << endl;

	/* if(newvind1==23620 || newvind2==23620){
		produceOutputs();
		cout << "FOUND" << endl;
		exit(1);
	} */
} // END make_divide_cell

bool Tissue::getDivisionPoints(const int cell, double &x1, double &x2, double &y1, double &y2, int &e1, int &e2)
{

	float division_angle_longest;
	float division_angle_external;
	float division_angle_external_degrees;

	int mv1, mv2;
	double dist, max_dist = -1;
	//get longest ditance
	for (int i = 0; i < cells[cell].num_vertices - 1; i++)
	{
		for (int j = i + 1; j < cells[cell].num_vertices; j++)
		{
			dist = distance(cells[cell].vertices[i], cells[cell].vertices[j]);
			if (dist > max_dist)
			{ // || max_dist < 0
				max_dist = dist;
				mv1 = cells[cell].vertices[i];
				mv2 = cells[cell].vertices[j];
			}
		}
	}

	//double center_x = (vertices[mv1].x + vertices[mv2].x)*0.5;
	//double center_y = (vertices[mv1].y + vertices[mv2].y)*0.5;
	//Instead, now line will go through centroid
	calculateCellCentroid(cells[cell]);
	double center_x = cells[cell].centroid_x;
	double center_y = cells[cell].centroid_y;

	double angle_hertwig = cells[cell].division_angle_longest == 0 ? 0 : cells[cell].division_angle_longest * (atan2(vertices[mv1].y - vertices[mv2].y, vertices[mv1].x - vertices[mv2].x) + 0.5 * M_PI);
	double random_angle = cells[cell].division_angle_random_noise == 0 ? 0 : cells[cell].division_angle_random_noise * (rand() % 360) * M_PI / 180;
	double externally_controlled_angle = cells[cell].division_angle_external == 0 ? 0 : cells[cell].division_angle_external * cells[cell].division_angle_external_degrees * M_PI / 180;
	double final_angle = angle_hertwig + random_angle + externally_controlled_angle;

	if (abs(final_angle - 0.5 * M_PI) <= NUMERIC_THRESHOLD)
	{
		//WARNING! This is a not very elegant fix of the problem of vertical lines (they have infinite slope), and can happen in other parts of the program
		final_angle += 0.01;
	}

	//Get point positions of new edge (new edge is first set to be longer than the longest distance between vertices, so that it intersects with at least two edges)
	x1 = center_x + cos(final_angle) * max_dist * 1.1;
	y1 = center_y + sin(final_angle) * max_dist * 1.1; // calculate an orthogonal line (big in excess to be sure that it cuts the polygon in 2 pieces) and look which edges it cuts
	x2 = center_x + cos(final_angle + M_PI) * max_dist * 1.1;
	y2 = center_y + sin(final_angle + M_PI) * max_dist * 1.1;

	//cout << "new angle is actually (before): " << 180*atan2(y1 - y2, x1 - x2)/M_PI << endl;
	//Find edges to cut
	StraightLine l1, l2, l3, l4;
	if (!findEdgesToCut(cell, x1, x2, y1, y2, e1, e2, l1, l2, l3, l4))
	{
		return false; //It will send a warning but the program will continue
	}

	x1 = (l3.intercept - l1.intercept) / (l1.slope - l3.slope); //Now calculate the points where new edge actually crosses parent cell edges
	x2 = (l4.intercept - l1.intercept) / (l1.slope - l4.slope);
	y1 = l1.intercept + l1.slope * x1;
	y2 = l1.intercept + l1.slope * x2;

	//cout << "new angle is actually (after): " << 180*atan2(y1 - y2, x1 - x2)/M_PI << endl;
	return true;
}

bool Tissue::findEdgesToCut(const int cell, double x1, double x2, double y1, double y2, int &e1, int &e2, StraightLine &l1, StraightLine &l2, StraightLine &l3, StraightLine &l4)
{
	l1 = getLineFromEdge(x1, x2, y1, y2);
	for (int i = 0; i < cells[cell].num_vertices; i++)
	{
		l2 = getLineFromEdge(&this->edges[cells[cell].edges[i]]);
		if (lines_cross(l1, l2))
		{
			if (e1 < 0)
			{
				e1 = cells[cell].edges[i];
				l3 = l2;
			}
			else
			{
				e2 = cells[cell].edges[i];
				l4 = l2;
				break;
			}
		}
	} //end for find edges that intersect with division edge
	if (REPORT_OUT > 1 && (e1 < 0 || e2 < 0))
	{
		cout << ">> Error: new edge does not cut other cell edges. Cell: " << cell << " ****** " << endl;
		return false;
	}
	return true;
}

void Tissue::splitEdgeWithVertex(int e, int cell, int v)
{

	int i = 0;
	int n1 = edges[e].vertices[0];
	int n2 = edges[e].vertices[1];

	//Insert vertex in cell
	while (cells[cell].vertices[i] != n1)
		i++;
	if (cells[cell].vertices[(i + 1) % cells[cell].num_vertices] == n2)
	{
		i = (i + 1) % cells[cell].num_vertices;
	}
	for (int j = cells[cell].num_vertices; j > i; j--)
	{
		cells[cell].vertices[j] = cells[cell].vertices[(j - 1 + cells[cell].num_vertices) % cells[cell].num_vertices];
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
	edges[e2ind].tension = edges[e].tension; //This holds true because both edges will have the same angle
	edges[e2ind].type = edges[e].type;
	edges[e2ind].base_tension = edges[e].base_tension;
	edges[e].length = distance(v, n1);
	edges[e2ind].length = distance(v, n2);
	edges[e2ind].can_transition = edges[e].can_transition;

	vertices[v].edges[0] = e;
	vertices[v].edges[1] = e2ind;
	vertices[n2].edges[which(e, vertices[n2].edges, CELLS_PER_VERTEX)] = e2ind;
	cells[cell].edges[which(EMPTY_CONNECTION, cells[cell].edges, MAX_SIDES_PER_CELL)] = e2ind; //add edge to dividing cell

	int other_cell = edges[e].cells[0] == cell ? edges[e].cells[1] : edges[e].cells[0]; //neighbour cell
	if (other_cell != EMPTY_CONNECTION)
	{
		cells[other_cell].edges[cells[other_cell].num_vertices] = e2ind;		   //add edge to the neighbour cell
		i = which(n1, cells[other_cell].vertices, cells[other_cell].num_vertices); //neighbour cell also needs the new vertex
		if (cells[other_cell].vertices[(i + 1) % cells[other_cell].num_vertices] == n2)
		{
			i = (i + 1) % cells[other_cell].num_vertices;
		}
		for (int j = cells[other_cell].num_vertices; j > i; j--)
		{
			cells[other_cell].vertices[j] = cells[other_cell].vertices[(j - 1 + cells[other_cell].num_vertices) % cells[other_cell].num_vertices];
		}
		cells[other_cell].vertices[i] = v;
		cells[other_cell].num_vertices++;
	}
} //End splitEdge

void Tissue::make_t1(Rearrangement &r)
{

	int edge = r.element_index;
	if (edges[edge].length > t1_transition_critical_distance && r.type != RearrangementType::random_t1)
		return; //Check that condition is still true (other rearrangements could have taken place since detection)

	Vertex *v1 = &vertices[edges[edge].vertices[0]];
	Vertex *v2 = &vertices[edges[edge].vertices[1]];
/* 	bool scrutiny = false;
	if(counter_moves_accepted == 149070349){scrutiny=true;}
	//if(edge == 23729 || v1->ind == 15595 || v2->ind ==15595 || v1->ind == 15864 || v2->ind ==15864  ){scrutiny=true;}
 	if(scrutiny){
		cout << "scrut A\n"<<endl;
		checkCellConsistency(7998);
		checkCellConsistency(8397);
	}  */
	if (contains(EMPTY_CONNECTION, v1->cells, CELLS_PER_VERTEX) ||
		contains(EMPTY_CONNECTION, v2->cells, CELLS_PER_VERTEX) ||
		contains(EMPTY_CONNECTION, v1->edges, CELLS_PER_VERTEX) ||
		contains(EMPTY_CONNECTION, v2->edges, CELLS_PER_VERTEX))
		return;

	double old_x1 = v1->x, old_x2 = v2->x, old_y1 = v1->y, old_y2 = v2->y, old_length = edges[edge].length;

	//Get the 2 cells that are common to both vertices and cell that is exclusive to v1
	int common_cell1 = -1, common_cell2 = -1, only_v1 = -1, only_v2 = -1; //, only_v1_indv1 = -1, only_v2_indv2 = -1;
	t1_getCellRelationships(v1, v2, common_cell1, common_cell2, only_v1, only_v2);
	Cell *cc1 = &this->cells[common_cell1]; //These pointers to Cell structures will be used later
	Cell *cc2 = &this->cells[common_cell2];
	Cell *sp1 = &this->cells[only_v1];
	Cell *sp2 = &this->cells[only_v2];

	if (cc1->num_vertices < 4 || cc2->num_vertices < 4 ||
		sp1->num_vertices >= MAX_SIDES_PER_CELL || sp2->num_vertices >= MAX_SIDES_PER_CELL)
		return; //One of the cells can't lose any other vertex
	if(ONLY_ONE_CONTACT){
		for(int s1e = 0; s1e < sp1->num_vertices; s1e++){
			for(int s2e = 0; s2e < sp2->num_vertices; s2e++){
				if(sp1->edges[s1e] == sp2->edges[s2e]){
					if (REPORT_OUT > 1)
						cout << ">> Error in t1 transition: Would set more than 1 contact between cells: v1=" << v1->ind <<
						" v2=" << v2->ind << " cell s2:" << sp2->ind << " cell s1:" << sp1->ind <<
						" cell cc2:" << cc2->ind << " cell cc1:" << cc2->ind <<
						" movements:" << counter_moves_accepted << endl;
						return;
				}
			}
		}
	}
	//t1_rotate_edge(v1, v2, edge, cc1, cc2, sp1, sp2);
	t1_rotate_edge90degrees(v1, v2, edge);
	//Determine which cell goes with which point. Compare new point locations of v1 and v2 with all the neighbours of v1 and v2.

	double dist1 = t1_get_dist_sum(v1, v2, cc1, cc2, sp1, sp2);
	if (dist1 >= 0)
	{
		int aux = common_cell1;
		common_cell1 = common_cell2;
		common_cell2 = aux;
		cc1 = &this->cells[common_cell1];
		cc2 = &this->cells[common_cell2];
	}

	//Now we know that v1 is closest to common_cell1, and v2 should be closer to common_cell2. Time to update everything
	//1) Cell exclussive of v1 is going to touch v2, and cell exclussive of v2 is going to touch v1
	v1->cells[which(common_cell2, v1->cells, CELLS_PER_VERTEX)] = only_v2;
	v2->cells[which(common_cell1, v2->cells, CELLS_PER_VERTEX)] = only_v1;

	//2)Remove edge from cells that were common to both vertices
	for (int i = which(edge, cc1->edges, MAX_SIDES_PER_CELL); i < cc1->num_vertices; i++)
		cc1->edges[i] = i == MAX_SIDES_PER_CELL - 1 ? EMPTY_CONNECTION : cc1->edges[i + 1];
	for (int i = which(edge, cc2->edges, MAX_SIDES_PER_CELL); i < cc2->num_vertices; i++)
		cc2->edges[i] = i == MAX_SIDES_PER_CELL - 1 ? EMPTY_CONNECTION : cc2->edges[i + 1];

	//3) Remove one of the vertices from each of the cells that were common to both vertices
	for (int i = which(v2->ind, cc1->vertices, cc1->num_vertices); i < cc1->num_vertices; i++)
		cc1->vertices[i] = i == MAX_SIDES_PER_CELL - 1 ? EMPTY_CONNECTION : cc1->vertices[i + 1];
	cc1->num_vertices--;

	for (int i = which(v1->ind, cc2->vertices, cc2->num_vertices); i < cc2->num_vertices; i++)
		cc2->vertices[i] = i == MAX_SIDES_PER_CELL - 1 ? EMPTY_CONNECTION : cc2->vertices[i + 1];
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

	//8)Update edge lengths, areas and tension.
	t1_update_sizes(v1, v2, edge);

	//9) Update cells in edge
	edges[edge].cells[0] = sp1->ind;
	edges[edge].cells[1] = sp2->ind;

	setEdgeTension(edge); //Tension depends on angle and on values of neighboring cells
	this->counter_t1++;

/*  	if(scrutiny){
		cout << "scrut B\n"<<endl;
		checkCellConsistency(7998);
		checkCellConsistency(8397);
	} */
	if (REPORT_OUT > 1){
		//double cent_x = 0.5*(v1->x + v2->x);
		//double cent_y = 0.5*(v1->y + v2->y);
		//printLine("T1", v1->ind, v2->ind, cent_x, cent_y, static_cast<int>(edges[edge].type));
		printLine_t1("T1", v1->ind, v2->ind, old_x1, old_x2, old_y1, old_y2, v1->x, v2->x, v1->y, v2->y, static_cast<int>(edges[edge].type));
	}
} //End make_t1

void Tissue::t1_update_sizes(Vertex *v1, Vertex *v2, int edge)
{

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v1->edges[i] != edge && v1->edges[i] != EMPTY_CONNECTION)
		{
			this->edges[v1->edges[i]].length = distance(this->edges[v1->edges[i]].vertices[0], this->edges[v1->edges[i]].vertices[1]);
			setEdgeTension(v1->edges[i]);
		}
		if (v2->edges[i] != edge && v2->edges[i] != EMPTY_CONNECTION)
		{
			this->edges[v2->edges[i]].length = distance(this->edges[v2->edges[i]].vertices[0], this->edges[v2->edges[i]].vertices[1]);
			setEdgeTension(v2->edges[i]);
		}
	}

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v1->cells[i] != EMPTY_CONNECTION)
		{
			this->cells[v1->cells[i]].area = calculateCellArea(this->cells[v1->cells[i]]);
			this->cells[v1->cells[i]].perimeter = calculateCellPerimeter(this->cells[v1->cells[i]]);
		}
		if (v2->cells[i] != EMPTY_CONNECTION)
		{
			this->cells[v2->cells[i]].area = calculateCellArea(this->cells[v2->cells[i]]);
			this->cells[v2->cells[i]].perimeter = calculateCellPerimeter(this->cells[v2->cells[i]]);
		}
	}
} // End update_sizes

//7) Add v1 to cell that was specific to v2, and v2 to cell that was specific to v1.
//Insertion must keep order of vertices in each cell
//Search rm_from_v2 vertex in sp2 cell (this vertex now is neighbour v1)
//Since vertex rm_from_v2 is touching common_cell1 (which now is specific to vertex v1), v1 should be added between this vertex and v2.
void Tissue::t1_add_vertices_to_cells(Vertex *v1, Vertex *v2, Cell *sp2, int remove_from_v2)
{

	int move_from;
	for (int i = 0; i < sp2->num_vertices; i++)
	{
		if (sp2->vertices[i] == remove_from_v2)
		{
			if (sp2->vertices[(i + 1) % sp2->num_vertices] == v2->ind)
			{
				move_from = (i + 1) % sp2->num_vertices;
			}
			else if (sp2->vertices[(i - 1 + sp2->num_vertices) % sp2->num_vertices] == v2->ind)
			{
				move_from = i;
			}
			else
			{
				if (REPORT_OUT > 1)
					cout << ">> Error in t1 transition: unable to insert vertex in cell: v1=" << v1->ind << " v2=" << v2->ind << " cell:" << sp2->ind << endl;
			}
			break;
		}
	}

	for (int i = sp2->num_vertices; i > move_from; i--)
	{
		sp2->vertices[i] = sp2->vertices[i - 1];
	}

	sp2->vertices[move_from] = v1->ind;
	sp2->num_vertices++;

} // end t1_add_vertices_to_cells

void Tissue::t1_update_edges(Vertex *v1, Vertex *v2, int edge, int remove_from_v1, int remove_from_v2, int &e_remove_from_v1, int &e_remove_from_v2)
{
	Edge *ee;
	int e_ind_remove_from_v1 = -1, e_ind_remove_from_v2;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		ee = &this->edges[v1->edges[i]];
		if (v1->edges[i] != edge && contains(remove_from_v1, ee->vertices, 2))
		{
			ee->vertices[which(v1->ind, ee->vertices, 2)] = v2->ind;
			e_ind_remove_from_v1 = i;
			e_remove_from_v1 = ee->ind;
			break;
		}
	}
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		ee = &this->edges[v2->edges[i]];
		if (v2->edges[i] != edge && contains(remove_from_v2, ee->vertices, 2))
		{
			ee->vertices[which(v2->ind, ee->vertices, 2)] = v1->ind;
			e_ind_remove_from_v2 = i;
			e_remove_from_v2 = ee->ind;
			break;
		}
	}
	v1->edges[e_ind_remove_from_v1] = e_remove_from_v2;
	v2->edges[e_ind_remove_from_v2] = e_remove_from_v1;
} // end update edges

void Tissue::t1_inwards_update_neighbours(Vertex *v1, Vertex *v2, int edge, int common_cell1, int only_v1, int only_v2, int &remove_from_v1, int &remove_from_v2)
{
	//cout << "h2" << endl;
	int ind_remove_from_v1, ind_remove_from_v2;
	Vertex *neighbour;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		neighbour = &this->vertices[v1->neighbour_vertices[i]];
		if (contains(EMPTY_CONNECTION, neighbour->cells, CELLS_PER_VERTEX) && contains(only_v1, neighbour->cells, CELLS_PER_VERTEX) && !contains(common_cell1, neighbour->cells, CELLS_PER_VERTEX) && neighbour->ind != v2->ind)
		{ //LAST CONDITION IS NEEDED FOR ELEMENTS IN CORNERS, SINCE THE WRONG NEIGHBOUR CAN ALSO BE TOUCHING EMPTY_CONNECTION
			ind_remove_from_v1 = i;
			remove_from_v1 = neighbour->ind;
			neighbour->neighbour_vertices[which(v1->ind, neighbour->neighbour_vertices, CELLS_PER_VERTEX)] = v2->ind; //change neighbour's neighbours also
			break;
		}
	}
	//cout << "h2" << endl;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		neighbour = &this->vertices[v2->neighbour_vertices[i]];
		//cout << "nei: " << neighbour->ind << " cc1 " << common_cell1 << " sp2 " << only_v2 << endl;
		if (contains(common_cell1, neighbour->cells, CELLS_PER_VERTEX) && contains(only_v2, neighbour->cells, CELLS_PER_VERTEX) && neighbour->ind != v1->ind)
		{
			ind_remove_from_v2 = i;
			remove_from_v2 = neighbour->ind;
			neighbour->neighbour_vertices[which(v2->ind, neighbour->neighbour_vertices, CELLS_PER_VERTEX)] = v1->ind; //change neighbour's neighbours also
			break;
		}
	}
	//cout << "h3 " << remove_from_v1 << " " << ind_remove_from_v1 << " " << remove_from_v2 << " " << ind_remove_from_v2 << endl;
	v1->neighbour_vertices[ind_remove_from_v1] = remove_from_v2;
	v2->neighbour_vertices[ind_remove_from_v2] = remove_from_v1;
	//cout << "h4" << endl;

} // End update neighbours t1 inwards

void Tissue::t1_update_neighbours(Vertex *v1, Vertex *v2, int edge, int common_cell1, int common_cell2, int only_v1, int only_v2, int &remove_from_v1, int &remove_from_v2)
{

	//5) Update neighbour vertices
	//V1 has to remove a neighbour touching common_cell2 and only_v1 cells
	//V2 has to remove a neighbour touching common_cell1 and only_v2 cells
	//V1 has to add a neighbour of v2 that is touching common_cell1 and only_v2 cells (the one removed from v2)
	//V2 has to add a neighbour of v1 that is touching common_cell2 and only_v1 cells (the one removed from v1)
	//cout << "h2" << endl;
	int ind_remove_from_v1, ind_remove_from_v2;
	Vertex *neighbour;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		neighbour = &this->vertices[v1->neighbour_vertices[i]];
		if (contains(common_cell2, neighbour->cells, CELLS_PER_VERTEX) && contains(only_v1, neighbour->cells, CELLS_PER_VERTEX) && neighbour->ind != v2->ind)
		{
			ind_remove_from_v1 = i;
			remove_from_v1 = neighbour->ind;
			neighbour->neighbour_vertices[which(v1->ind, neighbour->neighbour_vertices, CELLS_PER_VERTEX)] = v2->ind; //change neighbour's neighbours also
			break;
		}
	}
	//cout << "h2" << endl;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		neighbour = &this->vertices[v2->neighbour_vertices[i]];
		//cout << "nei: " << neighbour->ind << " cc1 " << common_cell1 << " cc2 " << common_cell2 << " sp2 " << only_v2 << endl;
		if (contains(common_cell1, neighbour->cells, CELLS_PER_VERTEX) && contains(only_v2, neighbour->cells, CELLS_PER_VERTEX) && neighbour->ind != v1->ind)
		{
			ind_remove_from_v2 = i;
			remove_from_v2 = neighbour->ind;
			neighbour->neighbour_vertices[which(v2->ind, neighbour->neighbour_vertices, CELLS_PER_VERTEX)] = v1->ind; //change neighbour's neighbours also
			break;
		}
	}
	//cout << "h3 " << remove_from_v1 << " " << ind_remove_from_v1 << " " << remove_from_v2 << " " << ind_remove_from_v2 << endl;
	v1->neighbour_vertices[ind_remove_from_v1] = remove_from_v2;
	v2->neighbour_vertices[ind_remove_from_v2] = remove_from_v1;
	//cout << "h4" << endl;

} // End update neighbours

void Tissue::t1_rotate_edge90degrees(Vertex *v1, Vertex *v2, int edge)
{
	double center_x = (v1->x + v2->x) * 0.5;
	double center_y = (v1->y + v2->y) * 0.5;
	double angle = atan2(v1->y - v2->y, v1->x - v2->x);
	//Rotate angle 90ºC and calculate new positions (new edge will have a length of 1.5*T1_TRANSITION_CRITICAL_DISTANCE)
	v1->x = center_x + cos(angle - 0.5 * M_PI) * length_rotated_edge;
	v2->x = center_x - cos(angle - 0.5 * M_PI) * length_rotated_edge;

	v1->y = center_y + sin(angle - 0.5 * M_PI) * length_rotated_edge;
	v2->y = center_y - sin(angle - 0.5 * M_PI) * length_rotated_edge;
	edges[edge].length = distance(v1->ind, v2->ind);

} //end rotate_edge

void Tissue::t1_rotate_edge(Vertex *v1, Vertex *v2, int edge, Cell *c1, Cell *c2, Cell *sp1, Cell *sp2)
{

	int v1_in_sp1 = which(v1->ind, sp1->vertices, sp1->num_vertices);
	int v2_in_sp2 = which(v2->ind, sp2->vertices, sp2->num_vertices);
	double y1 = abs(vertices[sp1->vertices[(v1_in_sp1 + 1) % sp1->num_vertices]].y - vertices[sp1->vertices[(v1_in_sp1 - 1 + sp1->num_vertices) % sp1->num_vertices]].y);
	double x1 = abs(vertices[sp1->vertices[(v1_in_sp1 + 1) % sp1->num_vertices]].x - vertices[sp1->vertices[(v1_in_sp1 - 1 + sp1->num_vertices) % sp1->num_vertices]].x);
	double angle1 = atan2(y1, x1);

	double y2 = abs(vertices[sp2->vertices[(v2_in_sp2 + 1) % sp2->num_vertices]].y - vertices[sp2->vertices[(v2_in_sp2 - 1 + sp2->num_vertices) % sp2->num_vertices]].y);
	double x2 = abs(vertices[sp2->vertices[(v2_in_sp2 + 1) % sp2->num_vertices]].x - vertices[sp2->vertices[(v2_in_sp2 - 1 + sp2->num_vertices) % sp2->num_vertices]].x);
	double angle2 = atan2(y2, x2);

	double center_x = (v1->x + v2->x) * 0.5;
	double center_y = (v1->y + v2->y) * 0.5;
	double angle = 0.5 * (angle1 + angle2);
	//Rotate angle 90ºC and calculate new positions (new edge will have a length of 1.5*T1_TRANSITION_CRITICAL_DISTANCE)
	v1->x = center_x + cos(angle) * length_rotated_edge;
	v2->x = center_x - cos(angle) * length_rotated_edge;

	v1->y = center_y + sin(angle) * length_rotated_edge;
	v2->y = center_y - sin(angle) * length_rotated_edge;

	edges[edge].length = distance(v1->ind, v2->ind);
} //end rotate_edge

void Tissue::t1_getCellRelationships(Vertex *v1, Vertex *v2, int &common_cell1, int &common_cell2, int &only_v1, int &only_v2)
{
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (contains(v1->cells[i], v2->cells, CELLS_PER_VERTEX))
		{
			if (common_cell1 == -1)
			{
				common_cell1 = v1->cells[i];
			}
			else
			{
				common_cell2 = v1->cells[i];
			}
		}
		else
		{
			only_v1 = v1->cells[i];
			//only_v1_indv1 = i;
		}
	}
	//Get cell that is exclussive of v2
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v2->cells[i] != common_cell1 && v2->cells[i] != common_cell2)
		{
			only_v2 = v2->cells[i];
			//only_v2_indv2 = i;
		}
	}
} // end t1_getCellRelationships

//dist t1_inwards
double Tissue::t1_inwards_get_dist_sum(Vertex *v1, Vertex *v2, Cell *c1, Cell *s1, Cell *s2)
{

	int v1_in_cell1 = which(v1->ind, c1->vertices, c1->num_vertices);
	int v2_in_cell1 = which(v2->ind, c1->vertices, c1->num_vertices);
	int v1_in_s1 = which(v1->ind, s1->vertices, s1->num_vertices);
	int v2_in_s2 = which(v2->ind, s2->vertices, s2->num_vertices);

	//These integers will hold the two neighboring vertices (n1 and n2) of vertices v1 and v2, in c1 and c2 respectively, after transition, assuming c1 becomes specific to v1
	int v1_n1 = -1, v1_n2 = -1, v2_n1 = -1, v2_n2 = -1;
	float dist;
	if (c1->vertices[(v1_in_cell1 + 1) % c1->num_vertices] == v2->ind)
	{
		v1_n1 = c1->vertices[(v1_in_cell1 + 2) % c1->num_vertices];
		v1_n2 = c1->vertices[(v1_in_cell1 - 1 + c1->num_vertices) % c1->num_vertices];
	}
	else
	{
		v1_n1 = c1->vertices[(v1_in_cell1 + 1) % c1->num_vertices];
		v1_n2 = c1->vertices[(v1_in_cell1 - 2 + c1->num_vertices) % c1->num_vertices];
	}

	//get neighbours of v2 without using common_cell2
	if (s2->vertices[(v2_in_s2 + 1) % s2->num_vertices] == v1_n1 || s2->vertices[(v2_in_s2 + 1) % s2->num_vertices] == v1_n2)
	{
		v2_n2 = s2->vertices[(v2_in_s2 - 1 + s2->num_vertices) % s2->num_vertices];
	}
	else
	{
		v2_n2 = s2->vertices[(v2_in_s2 + 1) % s2->num_vertices];
	}

	if (s1->vertices[(v1_in_s1 + 1) % s1->num_vertices] == v1_n1 || s1->vertices[(v1_in_s1 + 1) % s1->num_vertices] == v1_n2)
	{
		v2_n1 = s1->vertices[(v1_in_s1 - 1 + s1->num_vertices) % s1->num_vertices];
	}
	else
	{
		v2_n1 = s1->vertices[(v1_in_s1 + 1) % s1->num_vertices];
	}

	dist = distance(v1->ind, v1_n1) + distance(v1->ind, v1_n2) + distance(v2->ind, v2_n1) + distance(v2->ind, v2_n2);

	StraightLine lv1_c1n1, lv1_c1n2, lv2_c2n1, lv2_c2n2, lv1v2, cell_edge;
	lv1_c1n1 = getLineFromEdge(v1, &this->vertices[v1_n1]);
	lv1_c1n2 = getLineFromEdge(v1, &this->vertices[v1_n2]);
	lv2_c2n1 = getLineFromEdge(v2, &this->vertices[v2_n1]);
	lv2_c2n2 = getLineFromEdge(v2, &this->vertices[v2_n2]);
	lv1v2 = getLineFromEdge(v1, v2);
	vector<StraightLine> lines = {lv1_c1n1, lv1_c1n2, lv2_c2n1, lv2_c2n2, lv1v2};
	vector<Cell *> neighbour_cells = {c1, s1, s2};
	int indices_avoid[6] = {v1->ind, v2->ind, v1_n1, v1_n2, v2_n1, v2_n2};

	Edge *ee;
	for (StraightLine l : lines)
	{
		for (Cell *cellptr : neighbour_cells)
		{
			for (int i = 0; i < cellptr->num_vertices; i++)
			{
				ee = &this->edges[cellptr->edges[i]];
				//can't be one of the edges that are going to move, and edges cannot touch
				//(contains(s1->ind, ee->cells, 2) && contains(c2->ind, ee->cells, 2)) || (contains(s2->ind, ee->cells, 2) && contains(c1->ind, ee->cells, 2))
				if (!(contains(ee->vertices[0], indices_avoid, 4) && contains(ee->vertices[1], indices_avoid, 4)) &&
					!contains(l.v1, ee->vertices, 2) && !contains(l.v2, ee->vertices, 2))
				{
					cell_edge = getLineFromEdge(ee);
					if (lines_cross(cell_edge, l))
						return -1;
				}
			}
		}
	}
	if (lines_cross(lv1_c1n1, lv2_c2n1) || lines_cross(lv1_c1n1, lv2_c2n2) || lines_cross(lv1_c1n2, lv2_c2n1) || lines_cross(lv1_c1n2, lv2_c2n2))
		return -1;
	return dist;
} ///dist t1_inwards

//c1 and c2 are cels common to v1 and v2 in original position. s1 and s2 are cells specific to v1 and v2, respectively, in original positions.
double Tissue::t1_get_dist_sum(Vertex *v1, Vertex *v2, Cell *c1, Cell *c2, Cell *s1, Cell *s2)
{

	int v1_in_cell1 = which(v1->ind, c1->vertices, c1->num_vertices);
	int v1_in_cell2 = which(v1->ind, c2->vertices, c2->num_vertices);
	int v2_in_cell1 = which(v2->ind, c1->vertices, c1->num_vertices);
	int v2_in_cell2 = which(v2->ind, c2->vertices, c2->num_vertices);
	//These integers will hold the two neighboring vertices (n1 and n2) of vertices v1 and v2, in c1 and c2 respectively, after transition, assuming c1 becomes specific to v1
	int v1_n1 = -1, v1_n2 = -1, v2_n1 = -1, v2_n2 = -1;

	float dist, dist2;
	if (c1->vertices[(v1_in_cell1 + 1) % c1->num_vertices] == v2->ind)
	{
		v1_n1 = c1->vertices[(v1_in_cell1 + 2) % c1->num_vertices];
		v1_n2 = c1->vertices[(v1_in_cell1 - 1 + c1->num_vertices) % c1->num_vertices];
	}
	else
	{
		v1_n1 = c1->vertices[(v1_in_cell1 + 1) % c1->num_vertices];
		v1_n2 = c1->vertices[(v1_in_cell1 - 2 + c1->num_vertices) % c1->num_vertices];
	}

	if (c2->vertices[(v2_in_cell2 + 1) % c2->num_vertices] == v1->ind)
	{
		v2_n1 = c2->vertices[(v2_in_cell2 + 2) % c2->num_vertices];
		v2_n2 = c2->vertices[(v2_in_cell2 - 1 + c2->num_vertices) % c2->num_vertices];
	}
	else
	{
		v2_n1 = c2->vertices[(v2_in_cell2 + 1) % c2->num_vertices];
		v2_n2 = c2->vertices[(v2_in_cell2 - 2 + c2->num_vertices) % c2->num_vertices];
	}

	dist = distance(v1->ind, v1_n1) + distance(v1->ind, v1_n2) + distance(v2->ind, v2_n1) + distance(v2->ind, v2_n2);
	dist2 = distance(v1->ind, v2_n1) + distance(v1->ind, v2_n2) + distance(v2->ind, v1_n1) + distance(v2->ind, v1_n2);
	if (dist < dist2)
	{
		return -1;
	}else{
		return 1;
	}
}

//Assumes that number of cells contacting both vertices <=3:
void Tissue::make_join_limit_edges(Rearrangement &r)
{
	int edge = r.element_index;
	int aux = 0; //Auxiliary variable used to store indices
	if (edges[edge].length > t1_transition_critical_distance)
		return; //Check that condition is still true (other rearrangements could have taken place since detection)
	if (edges[edge].cells[0] != EMPTY_CONNECTION)
		if (cells[edges[edge].cells[0]].num_vertices < 4)
			return;
	if (edges[edge].cells[1] != EMPTY_CONNECTION)
		if (cells[edges[edge].cells[1]].num_vertices < 4)
			return;

	Vertex *v1 = &vertices[edges[edge].vertices[0]];
	Vertex *v2 = &vertices[edges[edge].vertices[1]];

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
		if (v1->edges[i] != EMPTY_CONNECTION)
			aux++;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
		if (v2->edges[i] != EMPTY_CONNECTION && !contains(v2->edges[i], v1->edges, CELLS_PER_VERTEX))
			aux++;
	if (aux > 4)
		return;

	//V1 will be preserved, v2 removed
	v1->x = (v1->x + v2->x) * 0.5;
	v1->y = (v1->y + v2->y) * 0.5;

	//1) Remove old edge from v1 and v2
	v1->edges[which(edge, v1->edges, CELLS_PER_VERTEX)] = EMPTY_CONNECTION;

	//2) Remove old vertex from v1.neighbour vertices
	v1->neighbour_vertices[which(v2->ind, v1->neighbour_vertices, CELLS_PER_VERTEX)] = EMPTY_CONNECTION;

	//3) Update cells in vertices
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v2->cells[i] != EMPTY_CONNECTION && !contains(v2->cells[i], v1->cells, CELLS_PER_VERTEX))
		{ //If cell i of vertex 2 does not touch vertex 1
			v1->cells[which(EMPTY_CONNECTION, v1->cells, CELLS_PER_VERTEX)] = v2->cells[i];
		}
	}

	//4) Update edges and vertices in edges
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v2->edges[i] != EMPTY_CONNECTION && v2->edges[i] != edge)
		{ //If edge i of vertex 2 does not touch vertex 1
			v1->edges[which(EMPTY_CONNECTION, v1->edges, CELLS_PER_VERTEX)] = v2->edges[i];
			edges[v2->edges[i]].vertices[which(v2->ind, edges[v2->edges[i]].vertices, 2)] = v1->ind; //add index of v1 to cell, in the place where v2 was
		}																							 //end if edge connected to v2 is not in v1
	}																								 //end for edges in v2

	//5) Update neighbour vertices
	Vertex *vv;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v2->neighbour_vertices[i] != EMPTY_CONNECTION && v2->neighbour_vertices[i] != v1->ind && !contains(v2->neighbour_vertices[i], v1->neighbour_vertices, CELLS_PER_VERTEX))
		{
			vv = &this->vertices[v2->neighbour_vertices[i]];
			v1->neighbour_vertices[which(EMPTY_CONNECTION, v1->neighbour_vertices, CELLS_PER_VERTEX)] = vv->ind;
			vv->neighbour_vertices[which(v2->ind, vv->neighbour_vertices, CELLS_PER_VERTEX)] = v1->ind;
		}
	}

	//6) Remove old edge from cells
	for (int i = 0; i < VERTEX_PER_EDGE; i++)
	{
		if (edges[edge].cells[i] != EMPTY_CONNECTION)
		{
			aux = edges[edge].cells[i];
			removeConnectionCell(edge, cells[aux].edges, cells[aux].num_vertices);
		}
	}

	//7) Update vertices in cells
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v2->cells[i] != EMPTY_CONNECTION)
		{
			if (contains(v1->ind, cells[v2->cells[i]].vertices, cells[v2->cells[i]].num_vertices))
			{
				removeConnectionCell(v2->ind, cells[v2->cells[i]].vertices, cells[v2->cells[i]].num_vertices);
				cells[v2->cells[i]].num_vertices--;
			}
			else
			{
				cells[v2->cells[i]].vertices[which(v2->ind, cells[v2->cells[i]].vertices, cells[v2->cells[i]].num_vertices)] = v1->ind;
			}
		}
	}

	//8) Update edge lengths
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v1->edges[i] != EMPTY_CONNECTION)
			edges[v1->edges[i]].length = distance(edges[v1->edges[i]].vertices[0], edges[v1->edges[i]].vertices[1]);
	}

	//9) Update cell areas and perimeter
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v1->cells[i] != EMPTY_CONNECTION)
		{
			cells[v1->cells[i]].area = calculateCellArea(cells[v1->cells[i]]);
			cells[v1->cells[i]].perimeter = calculateCellPerimeter(cells[v1->cells[i]]);
		}
	}

	//10) Energies will be updated when moving vertices, not now

	//11) Take care of springs
	if (v2->spring != EMPTY_CONNECTION)
	{

		if (v1->spring != EMPTY_CONNECTION)
		{
			//cout << "Deleted vertex and substitute vertex have springs. v2(dead): " << v2->ind << ", v1 (substitute): " << v1->ind << ", v2 spring: " << v2->spring <<  ", v1 spring: " << v1->spring <<endl;
			springs[v2->spring].dead = true;
			//Removed after added stringEdges
/* 			int dead_vert = springs[v2->spring].vertices[0] == v2->ind ? springs[v2->spring].vertices[1] : springs[v2->spring].vertices[0];
			int other_vert = springs[v1->spring].vertices[0] == v1->ind ? springs[v1->spring].vertices[1] : springs[v1->spring].vertices[0];
			vertices[other_vert].x = (vertices[other_vert].x + vertices[dead_vert].x) * 0.5;
			vertices[other_vert].y = (vertices[other_vert].y + vertices[dead_vert].y) * 0.5;
			springs[v1->spring].length = distance(springs[v1->spring].vertices[0], springs[v1->spring].vertices[1]);
			vertices[dead_vert].dead = true;
			dead_vertices.push(dead_vert);
			this->num_vertices--; */
			this->num_springs--;
		}
		else
		{
			//cout << "Deleted vertex has springs, but substitute does not. v2(dead): " << v2->ind << ", v1 (substitute): " << v1->ind  << ", v2 spring: " << v2->spring <<  ", v1 spring: " << v1->spring << endl;
			v1->spring = v2->spring;
			springs[v2->spring].vertices[which(v2->ind, springs[v2->spring].vertices, 2)] = v1->ind;
			springs[v1->spring].length = distance(springs[v1->spring].vertices[0], springs[v1->spring].vertices[1]);
		}
	}

	v2->dead = true;
	dead_vertices.push(v2->ind);
	this->num_vertices--;
	edges[edge].dead = true;
	dead_edges.push(edge);
	this->num_edges--;

	//Energies will be updated when moving vertices, not now
	counter_edges_removed++;
	//cout << "REMOVED: " << v2->ind << " and replaced by: " << v1->ind << "; accepted moves: " << counter_moves_accepted << "; edges removed: " << counter_edges_removed << endl;
}

void Tissue::make_t2(Rearrangement &r)
{

	int cell = r.element_index;
	//cout << "A\n";
	if (cells[cell].area > t2_transition_critical_area || cells[cell].num_vertices != 3 || cells[cell].dead)
		return;

	//cout << cells[cell] << endl;
	int v1 = cells[cell].vertices[0], v2 = cells[cell].vertices[1], v3 = cells[cell].vertices[2];
	vertices[v2].dead = true;
	vertices[v3].dead = true;
	vertices[v1].x = (vertices[v1].x + vertices[v2].x + vertices[v3].x) / 3;
	vertices[v1].y = (vertices[v1].y + vertices[v2].y + vertices[v3].y) / 3;
	//cout << "B\n";
	//Remove edges between vertices of dead cell (still not remove; will be replaced later in v1 edges)
	int v2e = -1, v3e = -1, v23e = -1;
	cells[cell].dead = true;
	for (int i = 0; i < cells[cell].num_vertices; i++)
	{
		if (contains(v1, edges[cells[cell].edges[i]].vertices, 2))
		{
			if (contains(v2, edges[cells[cell].edges[i]].vertices, 2))
			{
				v2e = cells[cell].edges[i];
			}
			else if (contains(v3, edges[cells[cell].edges[i]].vertices, 2))
			{
				v3e = cells[cell].edges[i];
			}
		}
		else
		{
			v23e = cells[cell].edges[i];
		}
		edges[cells[cell].edges[i]].dead = true;
	}
	//cout << "C\n";
	//Replace neighbours
	int v2nei = -1, v3nei = -1;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		v2nei = vertices[v2].neighbour_vertices[i];
		if (v2nei == EMPTY_CONNECTION)
		{
			vertices[v1].neighbour_vertices[which(v2, vertices[v1].neighbour_vertices, CELLS_PER_VERTEX)] = EMPTY_CONNECTION;
			break;
		}
		else if (v2nei != v1 && v2nei != v3)
		{
			vertices[v1].neighbour_vertices[which(v2, vertices[v1].neighbour_vertices, CELLS_PER_VERTEX)] = vertices[v2].neighbour_vertices[i];
			vertices[v2nei].neighbour_vertices[which(v2, vertices[v2nei].neighbour_vertices, CELLS_PER_VERTEX)] = v1;
			break;
		}
	}
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		v3nei = vertices[v3].neighbour_vertices[i];
		if (v3nei == EMPTY_CONNECTION)
		{
			vertices[v1].neighbour_vertices[which(v3, vertices[v1].neighbour_vertices, CELLS_PER_VERTEX)] = EMPTY_CONNECTION;
			break;
		}
		else if (v3nei != v1 && v3nei != v2)
		{
			vertices[v1].neighbour_vertices[which(v3, vertices[v1].neighbour_vertices, CELLS_PER_VERTEX)] = vertices[v3].neighbour_vertices[i];
			vertices[v3nei].neighbour_vertices[which(v3, vertices[v3nei].neighbour_vertices, CELLS_PER_VERTEX)] = v1;
			break;
		}
	}

	//cout << "D\n";
	//change edges connecting neighbours to v2 and v3
	int echange2, echange3;
	if (v2nei != EMPTY_CONNECTION)
	{
		for (int i = 0; i < CELLS_PER_VERTEX; i++)
		{
			if (vertices[v2nei].edges[i] == EMPTY_CONNECTION)
				continue;
			if (contains(v2, edges[vertices[v2nei].edges[i]].vertices, 2))
			{
				echange2 = vertices[v2nei].edges[i];
				edges[echange2].vertices[which(v2, edges[echange2].vertices, 2)] = v1;
				edges[echange2].length = distance(v1, v2nei);
				vertices[v1].edges[which(v2e, vertices[v1].edges, CELLS_PER_VERTEX)] = echange2;
				break;
			}
		}
	}
	else
	{
		vertices[v1].edges[which(v2e, vertices[v1].edges, CELLS_PER_VERTEX)] = EMPTY_CONNECTION;
	}
	//cout << "E\n";
	if (v3nei != EMPTY_CONNECTION)
	{
		for (int i = 0; i < CELLS_PER_VERTEX; i++)
		{
			if (vertices[v3nei].edges[i] == EMPTY_CONNECTION)
				continue;
			if (contains(v3, edges[vertices[v3nei].edges[i]].vertices, 2))
			{
				echange3 = vertices[v3nei].edges[i];
				edges[echange3].vertices[which(v3, edges[echange3].vertices, 2)] = v1;
				edges[echange3].length = distance(v1, v3nei);
				vertices[v1].edges[which(v3e, vertices[v1].edges, CELLS_PER_VERTEX)] = echange3;
				break;
			}
		}
	}
	else
	{
		vertices[v1].edges[which(v3e, vertices[v1].edges, CELLS_PER_VERTEX)] = EMPTY_CONNECTION;
	}

	//Recalculate distance of v1 neighbour outside of dying cell
	int v1e;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		v1e = vertices[v1].edges[i];
		if (v1e != EMPTY_CONNECTION)
		{
			edges[v1e].length = distance(edges[v1e].vertices[0], edges[v1e].vertices[1]);
		}
	}

	int neicell2, neicell3, neicell23;
	neicell2 = edges[v2e].cells[0] == cell ? edges[v2e].cells[1] : edges[v2e].cells[0];
	neicell3 = edges[v3e].cells[0] == cell ? edges[v3e].cells[1] : edges[v3e].cells[0];
	neicell23 = edges[v23e].cells[0] == cell ? edges[v23e].cells[1] : edges[v23e].cells[0];
	//cout << "e\n";
	//Change vertices in cells

	if (neicell2 != EMPTY_CONNECTION)
	{

		//cout  << " neicell2 ind: " << neicell2 << " v2 " << v2 << " v2e: " << v2e <<endl;
		//cout << "neicell2 " << cells[neicell2]<< endl;
		removeConnectionCell(v2e, cells[neicell2].edges, cells[neicell2].num_vertices);
		removeConnectionCell(v2, cells[neicell2].vertices, cells[neicell2].num_vertices);
		cells[neicell2].num_vertices--;
	}
	//cout << "e2"<<endl;
	if (neicell3 != EMPTY_CONNECTION)
	{
		//cout<< " neicell3 ind: " << neicell3 << " v3 " << v3 << " v3e: " << v3e <<endl;
		//cout << "neicell3 " << cells[neicell3] << endl;
		removeConnectionCell(v3e, cells[neicell3].edges, cells[neicell3].num_vertices);
		removeConnectionCell(v3, cells[neicell3].vertices, cells[neicell3].num_vertices);
		cells[neicell3].num_vertices--;
	}
	//cout << "F\n";
	if (neicell23 != EMPTY_CONNECTION && neicell23 != neicell3 && neicell23 != neicell2)
	{
		removeConnectionCell(v23e, cells[neicell23].edges, cells[neicell23].num_vertices);
		removeConnectionCell(v3, cells[neicell23].vertices, cells[neicell23].num_vertices);
		cells[neicell23].num_vertices--;
		cells[neicell23].vertices[which(v2, cells[neicell23].vertices, cells[neicell23].num_vertices)] = v1;
	}
	else if (neicell23 == neicell3 && neicell23 != EMPTY_CONNECTION)
	{
		removeConnectionCell(v23e, cells[neicell23].edges, cells[neicell23].num_vertices);
		cells[neicell23].vertices[which(v2, cells[neicell23].vertices, cells[neicell23].num_vertices)] = v1;
	}
	else if (neicell23 == neicell2 && neicell23 != EMPTY_CONNECTION)
	{
		removeConnectionCell(v23e, cells[neicell23].edges, cells[neicell23].num_vertices);
		cells[neicell23].vertices[which(v3, cells[neicell23].vertices, cells[neicell23].num_vertices)] = v1;
	}
	//cout << "G\n";

	//Change cell in v1
	if (neicell23 != neicell2 && neicell23 != neicell3)
	{
		vertices[v1].cells[which(cell, vertices[v1].cells, CELLS_PER_VERTEX)] = neicell23;
	}
	else
	{
		vertices[v1].cells[which(cell, vertices[v1].cells, CELLS_PER_VERTEX)] = EMPTY_CONNECTION;
	}

	//Recalculate areas and perimeters in neighboring cells
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (vertices[v1].cells[i] != EMPTY_CONNECTION)
		{
			cells[vertices[v1].cells[i]].area = calculateCellArea(cells[vertices[v1].cells[i]]);
			cells[vertices[v1].cells[i]].perimeter = calculateCellPerimeter(cells[vertices[v1].cells[i]]);
		}
	}
	//cout << "H\n";
	//Add dead elements to dead element queues:
	dead_cells.push(cell);
	dead_vertices.push(v1);
	dead_vertices.push(v2);
	dead_edges.push(v2e);
	dead_edges.push(v3e);
	dead_edges.push(v23e);
	//cout << "I\n";
	int static_v;
	//Check if vertices were bound to springs
	if (vertices[v2].spring != EMPTY_CONNECTION)
	{
		springs[vertices[v2].spring].dead = true;
		this->num_springs--;
		static_v = edges[vertices[v2].spring].vertices[0] == v2 ? springs[vertices[v2].spring].vertices[1] : springs[vertices[v2].spring].vertices[0];
		vertices[static_v].dead = true;
		dead_vertices.push(static_v);
	}
	if (vertices[v3].spring != EMPTY_CONNECTION)
	{
		springs[vertices[v3].spring].dead = true;
		this->num_springs--;
		static_v = edges[vertices[v3].spring].vertices[0] == v3 ? springs[vertices[v3].spring].vertices[1] : springs[vertices[v3].spring].vertices[0];
		vertices[static_v].dead = true;
		dead_vertices.push(static_v);
	}
	//cout << "J\n";

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (vertices[v1].edges[i] != EMPTY_CONNECTION)
			setEdgeTension(vertices[v1].edges[i]);
	}
	this->num_cells--;
	this->num_vertices -= 2;
	this->num_edges -= 3;

	this->counter_t2++;

	if (control_cells_2sides)
	{
		int new_v, edge_to_split; //, cell_with_two_sides;
		double x, y;
		int neicell, celledge;
		for (int i = 0; i < cells[cell].num_vertices; i++)
		{
			celledge = cells[cell].edges[i];
			for (int j = 0; j < 2; j++)
			{
				neicell = edges[celledge].cells[j];
				if (neicell != EMPTY_CONNECTION && neicell != cell)
				{
					if (cells[neicell].num_vertices < 3)
					{
						if (REPORT_OUT > 1)
							cout << "  in T2: Killing cell of size 2: " << endl
								 << cells[neicell] << endl;
						killCellWith2Vertices(neicell);
					}
				}
			}
		}
	}
	if (REPORT_OUT > 1){
		printLine_t2("T2", cell, v1, vertices[v1].x, vertices[v1].y, static_cast<int>(cells[cell].type));
	}
} //End make transition T2

void Tissue::killCellWith2Vertices(int cell)
{

	if (cells[cell].dead)
	{
		//cout << "*<< 2 side cell already dead: " << cell << ">>*" << endl;
		return;
	}
	//Get neighbouring vertices
	int vnei1 = -1, vnei2 = -1, enei1 = -1, enei2 = -1; //enei1 will join vnei1 and vnei2; enei2 will die
	int auxcell1 = -1, auxcell2 = -1;
	Cell *c = &cells[cell];
	Edge *e1 = &edges[cells[cell].edges[0]];
	Edge *e2 = &edges[cells[cell].edges[1]];
	Vertex *v1 = &vertices[cells[cell].vertices[0]];
	Vertex *v2 = &vertices[cells[cell].vertices[1]];

	//Get cells
	auxcell1 = e1->cells[0] == cell ? e1->cells[1] : e1->cells[0];
	auxcell2 = e2->cells[0] == cell ? e2->cells[1] : e2->cells[0];
	//Get neighbour edges
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		enei1 = v1->edges[i];
		if (e1->ind != enei1 && e2->ind != enei1)
		{
			break;
		}
	}

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		enei2 = v2->edges[i];
		if (e1->ind != enei2 && e2->ind != enei2)
		{
			break;
		}
	}

	if (enei1 < 0 && enei2 < 0)
	{
		cout << ">>Error: When removing cell with 2 edges, did not find neighbouring vertices" << endl;
		return;
	}
	else if (enei1 < 0)
	{
		enei1 = enei2;
		enei2 = EMPTY_CONNECTION;
		v2 = &vertices[cells[cell].vertices[0]];
		v1 = &vertices[cells[cell].vertices[1]];
	}
	if (enei1 == EMPTY_CONNECTION)
	{
		c->dead = true;
		e1->dead = true;
		e2->dead = true;
		v2->dead = true;
		for (int i = 0; i < CELLS_PER_VERTEX; i++)
		{
			if (v1->neighbour_vertices[i] == v2->ind)
			{ //No which() and no break; might be more than one
				v1->neighbour_vertices[i] = EMPTY_CONNECTION;
			}
			if (v1->cells[i] == c->ind)
			{
				v1->cells[i] = EMPTY_CONNECTION;
			}
			if (v1->edges[i] == e1->ind | v1->edges[i] == e2->ind)
			{
				v1->cells[i] = EMPTY_CONNECTION;
			}
		}
		return;
	}

	//Get other vertices
	vnei1 = edges[enei1].vertices[0] == v1->ind ? edges[enei1].vertices[1] : edges[enei1].vertices[0];
	vnei2 = edges[enei2].vertices[0] == v2->ind ? edges[enei2].vertices[1] : edges[enei2].vertices[0];

	//Change neighbour vertices
	vertices[vnei1].neighbour_vertices[which(v1->ind, vertices[vnei1].neighbour_vertices, CELLS_PER_VERTEX)] = vnei2;
	vertices[vnei2].neighbour_vertices[which(v2->ind, vertices[vnei2].neighbour_vertices, CELLS_PER_VERTEX)] = vnei1;

	//Remove neighbour edge2
	vertices[vnei2].edges[which(enei2, vertices[vnei2].edges, CELLS_PER_VERTEX)] = enei1;
	edges[enei2].dead = true;
	dead_edges.push(enei2);
	num_edges--;

	//Update neighbour edge 1
	edges[enei1].vertices[which(v1->ind, edges[enei1].vertices, 2)] = vnei2;

	//Kill cell and elements
	c->dead = true;
	v1->dead = true;
	v2->dead = true;
	e1->dead = true;
	e2->dead = true;
	num_vertices -= 2;
	num_cells -= 2;
	num_edges -= 2;

	dead_cells.push(c->ind);
	dead_vertices.push(v1->ind);
	dead_vertices.push(v2->ind);
	dead_edges.push(e1->ind);
	dead_edges.push(e2->ind);

	//Remove all connections from cells
	if (auxcell1 != EMPTY_CONNECTION)
	{
		removeConnectionCell(v1->ind, cells[auxcell1].vertices, cells[auxcell1].num_vertices);
		removeConnectionCell(v2->ind, cells[auxcell1].vertices, cells[auxcell1].num_vertices - 1);
		removeConnectionCell(e1->ind, cells[auxcell1].edges, cells[auxcell1].num_vertices);
		removeConnectionCell(enei2, cells[auxcell1].edges, cells[auxcell1].num_vertices);
		cells[auxcell1].num_vertices -= 2;
	}

	if (auxcell2 != auxcell1 && auxcell2 != EMPTY_CONNECTION)
	{
		removeConnectionCell(v1->ind, cells[auxcell2].vertices, cells[auxcell2].num_vertices);
		removeConnectionCell(v2->ind, cells[auxcell2].vertices, cells[auxcell2].num_vertices - 1);
		removeConnectionCell(e2->ind, cells[auxcell2].edges, cells[auxcell2].num_vertices);
		removeConnectionCell(enei2, cells[auxcell2].edges, cells[auxcell2].num_vertices - 1);
		cells[auxcell2].num_vertices -= 2;
	}
	else if (auxcell2 == auxcell1)
	{
		removeConnectionCell(e1->ind, cells[auxcell2].edges, cells[auxcell2].num_vertices);
	}

	if (auxcell1 != EMPTY_CONNECTION)
	{
		if (cells[auxcell1].num_vertices < 3)
		{
			killCellWith2Vertices(auxcell1);
		}
	}
	if (auxcell2 != EMPTY_CONNECTION)
	{
		if (cells[auxcell2].num_vertices < 3)
		{
			killCellWith2Vertices(auxcell2);
		}
	}
	if (REPORT_OUT > 1)
		cout << "    KILLED SIZE 2 CELL: " << cell << endl
			 << "<<*" << endl;
	return;

} //Kill cell with 2 edges
//Assumes that array elements is at least 1 element longer than length
void Tissue::removeConnectionCell(int elm, int *elements, int length)
{
	int i = 0;
	while (elements[i] != elm)
		i++;
	do
	{
		elements[i] = elements[i + 1];
		i++;
	} while (i < length - 1);
	elements[length - 1] = EMPTY_CONNECTION;
} //remove element from array of cell connections

bool Tissue::check_if_edges_cross(int vertex)
{
	Vertex *v = &this->vertices[vertex];
	Edge *e1, *e2;
	Cell *c1;
	StraightLine l_e1, l_e2;

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{ //For each edge i touching the vertex of interest
		if (v->edges[i] != EMPTY_CONNECTION)
		{
			e1 = &this->edges[v->edges[i]];
			l_e1 = getLineFromEdge(e1);
			for (int j = 0; j < CELLS_PER_VERTEX; j++)
			{ //For each cell j touching the vertex of interest
				if (v->cells[j] != EMPTY_CONNECTION)
				{
					c1 = &this->cells[v->cells[j]];
					for (int k = 0; k < c1->num_vertices; k++)
					{ //For each edge k in cell j
						e2 = &this->edges[c1->edges[k]];
						if (!contains(e1->vertices[0], e2->vertices, 2) && !contains(e1->vertices[1], e2->vertices, 2))
						{ //edges i and j can't touch
							l_e2 = getLineFromEdge(e2);
							if (lines_cross(l_e1, l_e2))
								return true; //Edges cross
						}
					}
				}
			}
		}
	}

	return false; //Edges do not cross
}

inline StraightLine Tissue::getLineFromEdge(const Edge *e)
{
	Vertex *v1 = &this->vertices[e->vertices[0]];
	Vertex *v2 = &this->vertices[e->vertices[1]];
	StraightLine sl;
	sl.x1 = v1->x;
	sl.x2 = v2->x;
	sl.y1 = v1->y;
	sl.y2 = v2->y;
	sl.slope = (sl.y1 - sl.y2) / (sl.x1 - sl.x2);
	sl.intercept = sl.y2 - sl.x2 * sl.slope;
	sl.v1 = v1->ind;
	sl.v2 = v2->ind;
	sl.vertical = std::isinf(sl.slope);
	return sl;
}

inline StraightLine Tissue::getLineFromEdge(const Vertex *v1, const Vertex *v2)
{
	StraightLine sl;
	sl.x1 = v1->x;
	sl.x2 = v2->x;
	sl.y1 = v1->y;
	sl.y2 = v2->y;
	sl.slope = (sl.y1 - sl.y2) / (sl.x1 - sl.x2);
	sl.intercept = sl.y2 - sl.x2 * sl.slope;
	sl.v1 = v1->ind;
	sl.v2 = v2->ind;
	sl.vertical = std::isinf(sl.slope);
	return sl;
}

inline StraightLine Tissue::getLineFromEdge(double x1, double x2, double y1, double y2)
{
	StraightLine sl;
	sl.x1 = x1;
	sl.x2 = x2;
	sl.y1 = y1;
	sl.y2 = y2;
	sl.slope = (sl.y1 - sl.y2) / (sl.x1 - sl.x2);
	sl.intercept = sl.y2 - sl.x2 * sl.slope;
	sl.vertical = std::isinf(sl.slope);
	return sl;
}

inline bool Tissue::lines_cross(StraightLine &a, StraightLine &b)
{
	double x_intersect;
	if (a.vertical && b.vertical)
	{
		return false;
	}
	else if (a.vertical)
	{
		x_intersect = 0.5 * (a.x1 + a.x2);
	}
	else if (b.vertical)
	{
		x_intersect = 0.5 * (b.x1 + b.x2);
	}
	else
	{
		if (abs(a.slope - b.slope) < NUMERIC_THRESHOLD)
		{
			return false;
		}
		else
		{
			x_intersect = (b.intercept - a.intercept) / (a.slope - b.slope);
		}
	}
	return ((x_intersect <= a.x1 && x_intersect >= a.x2) || (x_intersect <= a.x2 && x_intersect >= a.x1)) &&
	 ((x_intersect <= b.x1 && x_intersect >= b.x2) || (x_intersect <= b.x2 && x_intersect >= b.x1));
}

std::vector<int> Tissue::getNeighbourCells(int cell)
{
	Cell *c1 = &this->cells[cell];
	Edge *e;
	std::vector<int> v;
	bool added_border = false;
	int aux;
	for (int i = 0; i < c1->num_vertices; i++)
	{
		e = &this->edges[c1->edges[i]];
		aux = e->cells[0] == c1->ind ? e->cells[1] : e->cells[0];
		if (aux != EMPTY_CONNECTION || !added_border)
		{
			v.push_back(aux);
			if (aux == EMPTY_CONNECTION)
				added_border = true;
		}
	}
	return v;
}

void Tissue::addSpringsAutomatically()
{
	bool vertices_have_spring[num_vertices];
	for (int i = 0; i < num_vertices; i++)
		vertices_have_spring[i] = false;
	float minx = vertices[0].x;
	float maxx = vertices[0].x;
	for (int i = 1; i < num_vertices; i++)
	{
		minx = vertices[i].x < minx ? vertices[i].x : minx;
		maxx = vertices[i].x > maxx ? vertices[i].x : maxx;
	}							  //min and max of the whole wing
	setHingeMinAndMaxPositions(); //min and max of hinge
	Edge e;
	for (int i = 0; i < num_edges; i++)
	{
		e = edges[i];
		if (e.type != EdgeType::tissue_boundary || e.dead)
			continue;
		if (!vertices_have_spring[e.vertices[0]])
			vertices_have_spring[e.vertices[0]] = AddSpringToVertex(e.vertices[0], minx, maxx);
		if (!vertices_have_spring[e.vertices[1]])
			vertices_have_spring[e.vertices[1]] = AddSpringToVertex(e.vertices[1], minx, maxx);
		if (!vertices_have_spring[e.vertices[0]])
			setStaticVertex(e.vertices[0]);
		if (!vertices_have_spring[e.vertices[1]])
			setStaticVertex(e.vertices[1]);
	}
	switch (spring_tension_mode)
	{
	case 1:
		setSpringTension_mode1();
		break;
	case 2:
		setSpringTension_mode2();
		break;
	case 3:
		setSpringTension_mode3();
		break;
	case 4:
		setSpringTension_mode4();
		break;
	case 0:
	default:
		break;
	}
}

void Tissue::setStaticVertex(int v)
{
	float xprop = (vertices[v].x - hinge_min_xpos) / (hinge_max_xpos - hinge_min_xpos);
	float yprop = (vertices[v].y - hinge_min_ypos) / (hinge_max_ypos - hinge_min_ypos);
	vertices[v].movable = (yprop > 0.5 && xprop < add_static_to_hinge) ? 0 : 1; //Assumes that blade is always to the right of hinge
	vertices[v].movable_x = vertices[v].movable_y = vertices[v].movable;

	//cout << "Setting vertex " << v << " as " << vertices[v].movable << ". xprop: " << xprop << ", yprop: " << yprop << ", hinge pos x: (" << hinge_min_xpos <<", "<< hinge_max_xpos<<") hinge pos y: ("<<hinge_min_ypos<<", "<<hinge_max_ypos <<")" << endl;
}

bool Tissue::AddSpringToVertex(int v, float minx, float maxx)
{
	float pos = (vertices[v].x - minx) / (maxx - minx);
	for (int minpos = NUM_SPRING_TYPES - 1; minpos >= 0; minpos--)
	{
		if (pos > spring_type_min_positions.val[minpos])
		{
			Edge e;
			e.dead = false;
			e.type = EdgeType::spring;
			e.cells[0] = EMPTY_CONNECTION;
			e.cells[1] = EMPTY_CONNECTION;
			e.tension = spring_type_constants.val[minpos];
			e.length = 0;
			e.ind = this->num_springs;
			int vnew = newVertex(vertices[v].x, vertices[v].y);
			vertices[vnew].movable = false;
			vertices[vnew].spring = e.ind;
			vertices[v].spring = e.ind;
			e.vertices[0] = v;
			e.vertices[1] = vnew;
			this->springs.push_back(e);
			this->num_springs++;
			return true;
		}
	}
	return false;
} //AddSpringToVertex

void Tissue::restoreHinge(){
	float sum_x = 0;
	int num = 0;
	// Calculate mean x coordinate
	for(Edge e : edges){
		if(e.dead || e.cells[0] == EMPTY_CONNECTION || e.cells[1] == EMPTY_CONNECTION )
			continue;
		if((cells[e.cells[0]].type == CellType::blade && cells[e.cells[1]].type == CellType::hinge) ||
		(cells[e.cells[1]].type == CellType::blade && cells[e.cells[0]].type == CellType::hinge) ||
		(cells[e.cells[0]].type == CellType::vein_blade && cells[e.cells[1]].type == CellType::vein_hinge) ||
		(cells[e.cells[1]].type == CellType::vein_blade && cells[e.cells[0]].type == CellType::vein_hinge)){
			num += 2;
			sum_x += vertices[e.vertices[0]].x + vertices[e.vertices[1]].x;
		}
	}// end calculate mean x coordinate of hinge-blade interface
	// Now change cell types accordingly (cell centroid < mean -> type=hinge, otherwise type=blade )
	sum_x /= num;
	for(Cell &c : cells){
		if(c.dead)
			continue;
		calculateCellCentroid(c);
		if(c.centroid_x <= sum_x){
			c.type = c.type == CellType::vein_blade || c.type == CellType::vein_hinge ? CellType::vein_hinge : CellType::hinge;
		}else{
			c.type = c.type == CellType::vein_blade || c.type == CellType::vein_hinge ? CellType::vein_blade : CellType::blade;
		}
	}// end change cell types depending on whether they are at one side or another from the mean x coordinate

}

void Tissue::restoreVeins(){
	std::vector<int> vnei;
	int vein_neighbours, num_neighbours;
    CellType newtypes[cells.size()];
	for(Cell c : cells){
			if(c.dead)
				continue;
			vnei = getNeighbourCells(c.ind);
			num_neighbours = vnei.size();
			vein_neighbours = 0;
			for(int nei : vnei){
				if(nei == EMPTY_CONNECTION){
					num_neighbours--;
				}else{
					if(cells[nei].type == CellType::vein_blade || cells[nei].type == CellType::vein_hinge)
						vein_neighbours++;
				}
			}//for cell neighbours
			if(vein_neighbours < num_neighbours/2.0){
				newtypes[c.ind] = c.type == CellType::hinge || c.type == CellType::vein_hinge ? CellType::hinge : CellType::blade;
			}else if(vein_neighbours > num_neighbours/2.0){
				newtypes[c.ind] = c.type == CellType::hinge || c.type == CellType::vein_hinge ? CellType::vein_hinge : CellType::vein_blade;
			}else{//Do not change type
				newtypes[c.ind] = c.type;
			}
	}//for cells
	//Now set new types
	for(Cell &c : cells){
		if(c.dead)
				continue;
		c.type = newtypes[c.ind];
	}
} //End restore veins
void Tissue::restoreShape(){
	restoreHinge();
	for(int i = 0; i < RESTORE_VEIN_ITERS; i++)
		restoreVeins();
}
void Tissue::makeVeinsThinner(int iters){
    std::vector<CellType> newtypes(cells.size(), CellType::vein_blade);
	int nei;
	for(int i = 0; i < iters; i++){
		for(Cell c: cells){
			if(c.dead || c.type == CellType::blade || c.type == CellType::hinge)
				continue;
			for(int e = 0; e < c.num_vertices; e++){
				nei = edges[c.edges[e]].cells[0] == c.ind ? edges[c.edges[e]].cells[1] : edges[c.edges[e]].cells[0];
				if(nei == EMPTY_CONNECTION)
					continue;
				if(cells[nei].type == CellType::vein_blade || cells[nei].type == CellType::vein_hinge){
					continue;
				}else{
					newtypes[c.ind] = cells[nei].type;
					break;
				}
			}
		}
		for(Cell&  c: cells){
			if(newtypes[c.ind] != CellType::vein_blade){
				c.type = newtypes[c.ind];
			}
		}
	}
}

void Tissue::emptyDivisions()
{
	divisionrecord_q aux;
	std::swap(aux, this->past_divisions);
}

bool Tissue::checkCellConsistency(int cell){
	bool consistent = true;
	Cell *c = &cells[cell];
	Edge *e;
	Vertex *v;
	int edge_in_vertices = 0, vertex_in_edges = 0;

	for(int i = 0; i < c->num_vertices; i++){
		if(c->vertices[i] == EMPTY_CONNECTION){
			consistent = false;
			cout << "__ vertices[i] empty when i < num_vertices. Cell: " << c->ind  << ", moves: " << counter_moves_accepted << endl;
		}
		if(c->edges[i] == EMPTY_CONNECTION){
			consistent = false;
			cout << "__ edges[i] empty when i < num_vertices. Cell: " << c->ind  << ", moves: " << counter_moves_accepted << endl;
		}
		//Check that all edges are in cell vertices
		for(int j = 0; j < c->num_vertices; j++){
			v = &vertices[c->vertices[j]];
			if(v->edges[0]==c->edges[i] || v->edges[1]==c->edges[i] || v->edges[2]==c->edges[i]){
				edge_in_vertices++;
			}
		}
		if(edge_in_vertices != 2){
			consistent = false;
			cout << "__ edges[i] not in 2 vertices. Cell: " << c->ind << ", moves: " << counter_moves_accepted << ", edges[i]: " << c->edges[i] << ", i: " << i << ", edge count: "<< edge_in_vertices << endl;
		}
		edge_in_vertices = 0;

		//Check that all vertices are present in 2 edges
		for(int j = 0; j < c->num_vertices; j++){
			e = &edges[c->edges[j]];
			if(e->vertices[0]==c->vertices[i] || e->vertices[1]==c->vertices[i]){
				vertex_in_edges++;
			}
		}
		if(vertex_in_edges != 2){
			consistent = false;
			cout << "__ vertices[i] not in 2 edges. Cell: " << c->ind << ", moves: " << counter_moves_accepted << ", vertices[i]: " << c->vertices[i] << ", i: " << i << ", vert count: "<< vertex_in_edges << endl;
		}
		vertex_in_edges = 0;

	}
	if(!consistent){
		produceOutputs("not_consistent" + to_string(cell));
	}
	return consistent;
} //check cell consistency

void Tissue::printLine(std::string name, int ind1, int ind2, double centroid_x, double centroid_y, int type){
	cout << ">name:" << name <<
	";v1:"<<ind1<<";v2:"<<ind2<<
	";centroid_x:"<<centroid_x<<";centroid_y:"<<
	centroid_y<<";type:"<<type<<";moves_accepted:"<<counter_moves_accepted<<endl;
}
void Tissue::printLine_t1(std::string name, int ind1, int ind2, double v1_x0, double v2_x0, double v1_y0, double v2_y0,
		double v1_x1, double v2_x1, double v1_y1, double v2_y1, int type){
			//T1_HEADER = "name\tmoves_accepted\ttype\tv1\v2\tx1_before\tx1_after\tx2_before\tx2_after\ty1_before\ty1_after\ty2_before\ty2_after\tangle"
			double vec1x = v1_x0 - v2_x0;
			double vec1y = v1_y0 - v2_y0;
			double vec2x = v1_x1 - v2_x1;
			double vec2y = v1_y1 - v2_y1;
			double dot = vec1x*vec2x + vec1y*vec2y;
			double det = vec1x*vec2y + vec1y*vec2x;
			double angle = atan2(det, dot);
			//Source for the previous: https://stackoverflow.com/questions/14066933/direct-way-of-computing-clockwise-angle-between-2-vectors
			cout << ">" + name << "\t" << counter_moves_accepted << "\t" << type << "\t" <<
			 v1_x0 << "\t" << v1_x1 << "\t" <<
			 v2_x0 << "\t" << v2_x1 << "\t" <<
			 v1_y0 << "\t" << v1_y1 << "\t" <<
			 v2_y0 << "\t" << v2_y1 << "\t" << angle << "\n";
		}
void Tissue::printLine_t2(std::string name, int ind1, int ind2, double centroid_x, double centroid_y, int type){
	cout <<  ">" + name <<
	"\t"<<ind1<<"\t"<<ind2<<
	"\t"<<centroid_x<<"\t"<<
	centroid_y<<"\t"<<type<<"\t"<<counter_moves_accepted<<endl;

		}
void Tissue::printLine_div(std::string name, int ind1, int ind2, double centroid_1x, double centroid_1y,
		double centroid_2x, double centroid_2y, int type){
	cout <<  ">" + name <<"\t"<<
	ind1<<"\t"<<ind2<<"\t"<<
	centroid_1x<<"\t"<<
	centroid_1y<<"\t"<<
	centroid_2x<<"\t"<<
	centroid_2y<<"\t"
	<<type<<"\t"<<counter_moves_accepted<<endl;

		}

std::string Tissue::getStats()
{
	if (REPORT_OUT == 0)
		return std::to_string(counter_moves_accepted);
	std::string s = "";
	s += "Move trials: " + std::to_string(counter_move_trials) + ", \n";
	s += "Moves accepted: " + std::to_string(counter_moves_accepted) + ", \n";
	s += "prop. accepted: " + std::to_string(counter_moves_accepted / float(counter_move_trials)) + ", \n";
	s += "Prop. Positive E accepted: " + std::to_string(counter_unfav_accepted / float(counter_unfav_accepted + counter_unfav_rejected)) + ", \n";
	s += "Prop. Negative E accepted: " + std::to_string(counter_favorable_accepted / float(counter_favorable_accepted + counter_favorable_rejected)) + ", \n";
	s += "T1 accepted: " + std::to_string(counter_t1) + ", \n";
	s += "T1 rejected: " + std::to_string(counter_t1_abortions) + ", \n";
	s += "T1 margin outwards: " + std::to_string(counter_t1_outwards) + ", \n";
	s += "T1 margin inwards: " + std::to_string(counter_t1_inwards) + ", \n";
	s += "T2: " + std::to_string(counter_t2) + ", \n";
	s += "Edges removed: " + std::to_string(counter_edges_removed) + ", \n";
	s += "Divisions: " + std::to_string(counter_divisions) + "\n";
	return s;
}

int Tissue::getCounterT1()
{
	return this->counter_t1;
}

void Tissue::writeCellsFile(std::string fname)
{
	fname = fname + CELLS_FILE_EXTENSION;
	ofstream fo;
	fo.open(fname);
	fo << num_cells << "\t" << MAX_SIDES_PER_CELL << "\n";
	for (Cell c : cells)
	{
		if (!c.dead)
		{
			for (int v = 0; v < MAX_SIDES_PER_CELL; v++)
			{
				fo << c.vertices[v] << "\t";
			}
			fo << int(c.type) << "\n";
		}
	}
	fo.close();
}

void Tissue::writeSpringsFile(std::string fname)
{
	fname = fname + SPRING_FILE_EXTENSION;
	ofstream fo;
	fo.open(fname);
	bool typefound = false;
	fo << num_springs << "\n";
	for (Edge s : springs)
	{
		if (!s.dead)
		{
			fo << s.vertices[0] << "\t" << s.vertices[1] << "\t";
			for (int minpos = 0; minpos < NUM_SPRING_TYPES; minpos++)
			{
				if (abs(spring_type_constants.val[minpos] - s.tension) < 0.001)
				{
					fo << minpos;
					typefound = true;
					break;
				}
			}
			if (!typefound)
			{
				fo << 0;
			}
			fo << "\n";
			typefound = false;
		}
	}
	fo.close();
}

void Tissue::writePointsFile(std::string fname)
{
	fname = fname + VERTEX_FILE_EXTENSION;
	ofstream fo;
	fo.open(fname);
	fo << num_vertices << "\n";
	for (Vertex v : vertices)
	{
		if (!v.dead)
			fo << v.x << "\t" << v.y << "\t" << v.ind << "\t" << (v.movable ? v.movable_x && v.movable_y? 1 : v.movable_x? 2 : 3 : 0) << "\n";
	}
	fo.close();
}

void Tissue::writeAllData(std::string fname)
{ //Writes tables with all data (including parameters for each element)
	ofstream of;
	of.open(fname + ".out");
	of << *this;
	of.close();
} //End  writeAll

void Tissue::writeCellDataTable(std::string fname)
{
	ofstream of;
	of.open(fname + CELLTAB_FILE_EXTENSION);
	of << CELL_HEADER;
	for (Cell c : cells)
	{
		if (!c.dead)
			of << c << "\n";
	}
	of.close();
}
//HEADER =  "ind\ttype\ttension\tlength\tstatic_vertex\tmovable_vertex\tx_static\ty_static\tx_movable\ty_movable"
void Tissue::writeSpringDataTable(std::string fname)
{
	ofstream of;
	of.open(fname + SPRINGTAB_FILE_EXTENSION);
	of << SPR_HEADER;
	int vmov, vs;
	for (Edge s : springs)
	{
		if (!s.dead)
			of << s.ind << "\t";
			of << static_cast<int>(s.type) << "\t";
			of << s.tension << "\t";
			of << s.length << "\t";
			if(vertices[s.vertices[0]].movable){
				vmov = s.vertices[0];
				vs = s.vertices[1];
			}else{
				vmov = s.vertices[1];
				vs = s.vertices[0];
			}
			of << (vertices[vs].y <= AP_compartment_limit ? "P" : "A") << "\t";
			of << vs << "\t";
			of << vmov << "\t";
			of << vertices[vs].x << "\t";
			of << vertices[vs].y << "\t";
			of << vertices[vmov].x << "\t";
			of << vertices[vmov].y << "\n";
	}
	of.close();
}

void Tissue::writeEdgeDataTable(std::string fname)
{
	ofstream of;
	of.open(fname + EDGE_FILE_EXTENSION);
	of << EDGE_HEADER;
	for (Edge e : edges)
	{
		if (!e.dead)
			of << e << "\n";
	}
	of.close();
}
void Tissue::writePointsDataTable(std::string fname)
{
	ofstream of;
	of.open(fname + POINTTAB_FILE_EXTENSION);
	of << VERTEX_HEADER;
	for (Vertex v : vertices)
	{
		if (!v.dead)
			of << v << "\n";
	}
	of.close();
}
void Tissue::printCelltypeParam(cell_type_param par, std::string name){
	cout << name << ": ";
	for(int i = 0; i < NUM_CELL_TYPES; i++) cout << i << "=" << par.val[i] << ",\t";
	cout << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//What follows is only useful to print easily
//
//Overloading of << operator for VERTICES. Useful to print
std::ostream &operator<<(std::ostream &out, const Vertex &v)
{
	int vmovable = v.movable ? v.movable_x && v.movable_y? 1 : v.movable_x? 2 : 3 : 0;
	out << v.ind << "\t" << v.x << "\t" << v.y << "\t" << v.energy << "\t" << vmovable << "\t" << v.spring << "\t" << v.moves_accepted << "\t" << v.moves_rejected << "\t";
	for (int i = 0; i < sizeof(v.cells) / sizeof(v.cells[0]); i++)
	{ // print cells touching vertex separated by ','
		out << v.cells[i] << ",";
	}
	out << "\t";
	for (int i = 0; i < sizeof(v.edges) / sizeof(v.edges[0]); i++)
	{ // print edges touching vertex separated by ','
		out << v.edges[i] << ",";
	}
	out << "\t";
	for (int i = 0; i < sizeof(v.neighbour_vertices) / sizeof(v.neighbour_vertices[0]); i++)
	{ // print neighbour vertices separated by ','
		out << v.neighbour_vertices[i] << ",";
	}

	return out;
}

//Overloading of << operator for CELLS. Useful to print
std::ostream &operator<<(std::ostream &out, const Cell &c)
{
	out << c.ind << "\t" << int(c.type) << "\t" << c.area << "\t" << c.preferred_area << "\t" << c.K << "\t" << c.perimeter << "\t" << c.perimeter_contractility << "\t";
	out << c.centroid_x << "\t" << c.centroid_y << "\t";
	out << c.base_preferred_area << "\t" << c.base_perimeter << "\t";
	out << c.division_angle_longest << "\t";
	out << c.division_angle_external << "\t";
	out << c.division_angle_random_noise << "\t";
	out << c.division_angle_external_degrees << "\t";
	out << c.max_area << "\t";
	out << c.cell_cycle_state << "\t";
	out << c.cell_cycle_limit << "\t";
	out << c.can_divide << "\t";
	out << c.num_vertices << "\t";
	for (int i = 0; i < c.num_vertices; i++)
	{ // print vertices of cell separated by ','
		out << c.vertices[i] << ",";
	}
	out << "\t";
	for (int i = 0; i < c.num_vertices; i++)
	{ // print edges of cell separated by ','
		out << c.edges[i] << ",";
	}
	out << "\t" << c.num_divisions;
	out << "\t" <<  static_cast<int>(c.vary_line_tension);
	out << "\t" <<  c.edge_angle_prop_external;
	out << "\t" << c.edge_angle_prop_uniform;
	out << "\t" << c.edge_angle_prop_maxangle;
	out << "\t" << c.edge_angle_prop_random;
	out << "\t" << c.edge_tension_external;
	out << "\t" << c.edge_maxangle;
	out << "\t" << c.edge_spatialmax_tension;
	out << "\t" << c.edge_spatialmin_tension;
	out << "\t" << c.father;
	out << "\t" << c.x_at_start;

	return out;
}

//Overloading of << operator for EDGES. Useful to print
std::ostream &operator<<(std::ostream &out, const Edge &e)
{
	out << e.ind << "\t" 
	<< int(e.type) << "\t" 
	<< e.length << "\t" 
	<< e.tension << "\t" 
	<< e.base_tension<< "\t";
	for (int i = 0; i < sizeof(e.vertices) / sizeof(e.vertices[0]); i++)
	{ //print vertices touching edge separated by ','
		out << e.vertices[i] << ",";
	}
	out << "\t";
	for (int i = 0; i < sizeof(e.cells) / sizeof(e.cells[0]); i++)
	{ //print cells touching edge separated by ','
		out << e.cells[i] << ",";
	}
	return out;
}

//Overloading of << operator for TISSUE. Useful to print
std::ostream &operator<<(std::ostream &out, const Tissue &t)
{
	out << t.simname << "\t-\t"
		<< "max. accepted moves: " << t.max_accepted_movements << ",\twrite every N moves: " << t.write_every_N_moves << "\n\n";
	std::string s = "";
	s += "Move trials: " + std::to_string(t.counter_move_trials);
	s += ", Moves accepted: " + std::to_string(t.counter_moves_accepted);
	s += ", prop. accepted: " + std::to_string(t.counter_moves_accepted / float(t.counter_move_trials));
	s += ", T1 accepted: " + std::to_string(t.counter_t1);
	s += ", T1 rejected: " + std::to_string(t.counter_t1_abortions);
	s += ", T2: " + std::to_string(t.counter_t2);
	out << s << "\n";
	//out << t.getStats() << "\n";
	out << "VERTICES:\n";
	out << VERTEX_HEADER;
	for (Vertex v : t.vertices)
	{
		if (!v.dead)
			out << v << "\n";
	}
	out << "\nCELLS:\n";
	out << CELL_HEADER;
	for (Cell c : t.cells)
	{
		if (!c.dead)
			out << c << "\n";
	}
	out << "\nEDGES:\n";
	out << EDGE_HEADER;
	for (Edge e : t.edges)
	{
		if (!e.dead)
			out << e << "\n";
	}
	out << "\nSPRINGS:\n";
	for (Edge e : t.springs)
	{
		if (!e.dead)
			out << e << "\n";
	}
	out << "****\n****\n";
	return out;
}

//Helper functions of general use.
//Returns true if integer n is present in array a, otherwise returns false.
inline bool contains(int n, const int *a, int len)
{
	for (int i = 0; i < len; i++)
	{
		if (n == a[i])
			return true;
	}
	return false;
}
//Returns -1 if integer n is not in array a, otherwise returns index of n in a.
inline int which(int n, const int *a, int len)
{
	for (int i = 0; i < len; i++)
	{
		if (n == a[i])
			return i;
	}
	return -1;
}

inline int count(int n, const int *a, int len)
{
	int cont = 0;
	for (int i = 0; i < len; i++)
	{
		if (n == a[i])
			cont++;
	}
	return cont;
}
//Returns index of FIRST element in array a that is not present in array b
inline int element_not_in(const int *a, const int *b, int len1, int len2)
{
	bool a_in_b;
	for (int i = 0; i < len1; i++)
	{
		a_in_b = false;
		for (int j = 0; j < len2; j++)
		{
			if (a[i] == b[j])
			{
				a_in_b = true;
				break;
			}
		}
		if (!a_in_b)
			return i;
	}
	return -1;
}