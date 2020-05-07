

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

//Default constructor
Tissue::Tissue() : num_cells(0), num_vertices(0), num_edges(0), counter_move_trials(0), counter_moves_accepted(0), counter_favorable_accepted(0), counter_favorable_rejected(0), counter_unfav_accepted(0), counter_unfav_rejected(0), counter_t1(0), counter_t1_abortions(0), counter_edges_removed(0), counter_divisions(0), counter_t2(0), counter_t1_outwards(0), counter_t1_inwards(0)
{

	simname = "";
	max_accepted_movements = 0;
	write_every_N_moves = 0;
	step_mode = false;

	set_default_simulation_params();
	setHingeMinAndMaxPositions();
}

//Constructor that reads vertex positions and cells from two different files, and initializes all variables using constants defined in VertexSystem.h
Tissue::Tissue(std::string starting_tissue_file, int max_accepted_movements, int write_every_N_moves, string simulname) : Tissue()
{

	this->simname = simulname == "" ? starting_tissue_file : simulname;
	this->max_accepted_movements = max_accepted_movements;
	this->write_every_N_moves = write_every_N_moves;
	step_mode = false;

	//Read file of vertices (indicates coordinates for each vertex)
	string vertexfile = starting_tissue_file + VERTEX_FILE_EXTENSION;
	std::ifstream fin_vertex;
	if (REPORT_OUT)
		cout << "reading .points file...\n";
	fin_vertex.open(vertexfile);
	initialize_vertices(fin_vertex);
	fin_vertex.close();
	if (REPORT_OUT)
		cout << ".points file read...\n";

	//Read file of cells (indicates vertices for each cell)
	try
	{
		string cellfile = starting_tissue_file + CELLS_FILE_EXTENSION;
		if (REPORT_OUT)
			cout << "reading .cells file...\n";
		ifstream fin_cells(cellfile);
		initialize_cells(fin_cells);
		fin_cells.close();
		if (REPORT_OUT)
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
		if (REPORT_OUT)
			cout << "reading .spr file...\n";
		ifstream fin_springs(springsfile);
		if (fin_springs.good())
		{
			initialize_springs(fin_springs);
			if (REPORT_OUT)
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
		if (REPORT_OUT)
			cout << msg << endl;
		if (REPORT_OUT)
			cout << "No spring file\n";
	}
	//No file with parameters, therefore use constants defined in VertexSystem.h.
	//Also initializes cell area and vertex energy
	set_default_params();
	if (REPORT_OUT)
		cout << "parameters set to cells\n";
	setHingeMinAndMaxPositions();
}

Tissue::Tissue(std::string starting_tissue_file, std::string params_file, int max_accepted_movements, int write_every_N_moves, string simulname) : num_cells(0), num_vertices(0), num_edges(0), counter_move_trials(0), counter_moves_accepted(0), counter_favorable_accepted(0), counter_favorable_rejected(0), counter_unfav_accepted(0), counter_unfav_rejected(0), counter_t1(0), counter_t1_abortions(0), counter_edges_removed(0), counter_divisions(0), counter_t2(0), counter_t1_outwards(0), counter_t1_inwards(0)
{
	cout << "reading .vp file...\n";
	this->simname = simulname == "" ? starting_tissue_file : simulname;
	this->max_accepted_movements = max_accepted_movements;
	this->write_every_N_moves = write_every_N_moves;
	step_mode = false;

	if (REPORT_OUT)
		cout << "reading .vp file...\n";
	initialize_params(params_file);
	if (REPORT_OUT)
		cout << ".vp file read...\n";
	//Read file of vertices (indicates coordinates for each vertex)
	string vertexfile = starting_tissue_file + VERTEX_FILE_EXTENSION;
	std::ifstream fin_vertex;
	if (REPORT_OUT)
		cout << "reading .points file...\n";
	fin_vertex.open(vertexfile);
	initialize_vertices(fin_vertex);
	fin_vertex.close();
	if (REPORT_OUT)
		cout << ".points file read...\n";

	//Read file of cells (indicates vertices for each cell)
	try
	{
		string cellfile = starting_tissue_file + CELLS_FILE_EXTENSION;
		if (REPORT_OUT)
			cout << "reading .cells file...\n";
		ifstream fin_cells(cellfile);
		initialize_cells(fin_cells);
		fin_cells.close();
		if (REPORT_OUT)
			cout << ".cells file read...\n";
	}
	catch (const char *msg)
	{
		cout << msg << endl;
		exit(1);
	}
	//Initialize edges from vertices and cells
	if (REPORT_OUT)
		cout << "Initializing edges...\n";
	initialize_edges();
	try
	{
		string springsfile = starting_tissue_file + SPRING_FILE_EXTENSION;
		ifstream fin_springs(springsfile);
		if (fin_springs.good())
		{
			if (REPORT_OUT)
				cout << "reading .spr file...\n";
			initialize_springs(fin_springs);
			if (REPORT_OUT)
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
		if (REPORT_OUT)
			cout << msg << endl;
		if (REPORT_OUT)
			cout << "No spring file\n";
	}
	//No file with parameters, therefore use constants defined in VertexSystem.h.
	//Also initializes cell area and vertex energy
	set_default_params();
	if (REPORT_OUT)
		cout << "parameters set to cells\n";
	setHingeMinAndMaxPositions();
}

void Tissue::set_default_simulation_params()
{
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
	xcoord_decrease_exponent = XCOORD_DECREASE_EXPONENT;

	energy_term1 = ENERGY_TERM1;
	energy_term2 = ENERGY_TERM2;
	energy_term3 = ENERGY_TERM3;

	line_tension.insert(pair<CellType, double>(CellType::blade, LINE_TENSION_BLADE));
	line_tension.insert(pair<CellType, double>(CellType::hinge, LINE_TENSION_HINGE));
	line_tension.insert(pair<CellType, double>(CellType::vein_hinge, LINE_TENSION_VEIN_HINGE));
	line_tension.insert(pair<CellType, double>(CellType::vein_blade, LINE_TENSION_VEIN_BLADE));
	line_tension_tissue_boundary.insert(pair<CellType, double>(CellType::blade, LINE_TENSION_TISSUE_BOUNDARY));
	line_tension_tissue_boundary.insert(pair<CellType, double>(CellType::hinge, LINE_TENSION_TISSUE_BOUNDARY));
	line_tension_tissue_boundary.insert(pair<CellType, double>(CellType::vein_hinge, LINE_TENSION_TISSUE_BOUNDARY));
	line_tension_tissue_boundary.insert(pair<CellType, double>(CellType::vein_blade, LINE_TENSION_TISSUE_BOUNDARY));

	perimeter_contract.insert(pair<CellType, double>(CellType::blade, PERIMETER_CONTRACT_BLADE));
	perimeter_contract.insert(pair<CellType, double>(CellType::hinge, PERIMETER_CONTRACT_HINGE));
	perimeter_contract.insert(pair<CellType, double>(CellType::vein_hinge, PERIMETER_CONTRACT_BLADE));
	perimeter_contract.insert(pair<CellType, double>(CellType::vein_hinge, PERIMETER_CONTRACT_HINGE));

	t1_transition_critical_distance = T1_TRANSITION_CRITICAL_DISTANCE;
	t2_transition_critical_area = T2_TRANSITION_CRITICAL_AREA;
	max_edge_length = MAX_EDGE_LENGTH;

	preferred_area_initial.insert(pair<CellType, double>(CellType::blade, PREFERRED_AREA_INITIAL));
	preferred_area_initial.insert(pair<CellType, double>(CellType::hinge, PREFERRED_AREA_INITIAL));
	preferred_area_initial.insert(pair<CellType, double>(CellType::vein_blade, PREFERRED_AREA_INITIAL));
	preferred_area_initial.insert(pair<CellType, double>(CellType::vein_hinge, PREFERRED_AREA_INITIAL));

	preferred_area_final.insert(pair<CellType, double>(CellType::blade, PREFERRED_AREA_FINAL));
	preferred_area_final.insert(pair<CellType, double>(CellType::hinge, PREFERRED_AREA_FINAL));
	preferred_area_final.insert(pair<CellType, double>(CellType::vein_blade, PREFERRED_AREA_FINAL));
	preferred_area_final.insert(pair<CellType, double>(CellType::vein_hinge, PREFERRED_AREA_FINAL));

	division_angle_random_noise.insert(pair<CellType, double>(CellType::blade, DIVISION_ANGLE_RANDOM_NOISE));
	division_angle_random_noise.insert(pair<CellType, double>(CellType::hinge, DIVISION_ANGLE_RANDOM_NOISE));
	division_angle_random_noise.insert(pair<CellType, double>(CellType::vein_blade, DIVISION_ANGLE_RANDOM_NOISE));
	division_angle_random_noise.insert(pair<CellType, double>(CellType::vein_hinge, DIVISION_ANGLE_RANDOM_NOISE));

	division_angle_longest_axis.insert(pair<CellType, double>(CellType::blade, DIVISION_ANGLE_LONGEST_AXIS));
	division_angle_longest_axis.insert(pair<CellType, double>(CellType::hinge, DIVISION_ANGLE_LONGEST_AXIS));
	division_angle_longest_axis.insert(pair<CellType, double>(CellType::vein_blade, DIVISION_ANGLE_LONGEST_AXIS));
	division_angle_longest_axis.insert(pair<CellType, double>(CellType::vein_hinge, DIVISION_ANGLE_LONGEST_AXIS));

	division_angle_external.insert(pair<CellType, double>(CellType::blade, DIVISION_ANGLE_EXTERNAL));
	division_angle_external.insert(pair<CellType, double>(CellType::hinge, DIVISION_ANGLE_EXTERNAL));
	division_angle_external.insert(pair<CellType, double>(CellType::vein_blade, DIVISION_ANGLE_EXTERNAL));
	division_angle_external.insert(pair<CellType, double>(CellType::vein_hinge, DIVISION_ANGLE_EXTERNAL));

	division_angle_external_degrees.insert(pair<CellType, double>(CellType::blade, DIVISION_ANGLE_EXTERNAL_DEGREES));
	division_angle_external_degrees.insert(pair<CellType, double>(CellType::hinge, DIVISION_ANGLE_EXTERNAL_DEGREES));
	division_angle_external_degrees.insert(pair<CellType, double>(CellType::vein_blade, DIVISION_ANGLE_EXTERNAL_DEGREES));
	division_angle_external_degrees.insert(pair<CellType, double>(CellType::vein_hinge, DIVISION_ANGLE_EXTERNAL_DEGREES));

	spring_type_constants.insert(pair<int, double>(0, SPRING_CONSTANT));
	spring_type_constants.insert(pair<int, double>(1, SPRING_CONSTANT));
	spring_type_constants.insert(pair<int, double>(2, SPRING_CONSTANT));
	spring_type_constants.insert(pair<int, double>(3, SPRING_CONSTANT));
	//spring_constant = SPRING_CONSTANT;

	max_cell_area.insert(pair<CellType, double>(CellType::blade, MAX_CELL_AREA));
	max_cell_area.insert(pair<CellType, double>(CellType::hinge, MAX_CELL_AREA));
	max_cell_area.insert(pair<CellType, double>(CellType::vein_blade, MAX_CELL_AREA));
	max_cell_area.insert(pair<CellType, double>(CellType::vein_hinge, MAX_CELL_AREA));
	//max_cell_area = MAX_CELL_AREA;

	cell_cycle_limit.insert(pair<CellType, double>(CellType::blade, CELL_CYCLE_LIMIT));
	cell_cycle_limit.insert(pair<CellType, double>(CellType::hinge, CELL_CYCLE_LIMIT));
	cell_cycle_limit.insert(pair<CellType, double>(CellType::vein_blade, CELL_CYCLE_LIMIT));
	cell_cycle_limit.insert(pair<CellType, double>(CellType::vein_hinge, CELL_CYCLE_LIMIT));

	vary_line_tension.insert(pair<CellType, double>(CellType::blade, -1.0));
	vary_line_tension.insert(pair<CellType, double>(CellType::hinge, -1.0));
	vary_line_tension.insert(pair<CellType, double>(CellType::vein_blade, -1.0));
	vary_line_tension.insert(pair<CellType, double>(CellType::vein_hinge, -1.0));

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

	spring_type_min_positions.insert(pair<int, double>(0, 0.5));
	spring_type_min_positions.insert(pair<int, double>(1, 0.5));
	spring_type_min_positions.insert(pair<int, double>(2, 0.5));
	spring_type_min_positions.insert(pair<int, double>(3, 0.5));
	add_static_to_hinge = -1;
}

void Tissue::readNewParameters(std::string filename)
{
	initialize_params(filename);
	set_default_params();
	setHingeMinAndMaxPositions();
}

double Tissue::read_real_par(std::vector<std::string>::iterator &it)
{
	while (it->at(0) != '>')
		it++;
	it++;
	return stod(*it);
}

cell_type_param Tissue::read_celltype_par(std::vector<std::string>::iterator &it, std::string::size_type sz)
{
	std::string s;
	CellType celltype;
	double auxd;
	cell_type_param res;

	while (it->at(0) != '>')
		it++;
	it++;

	while (it->at(0) != '<')
	{
		s = *it;
		celltype = static_cast<CellType>(stoi(s, &sz));
		s = s.substr(sz);
		s.erase(0, 1);
		auxd = stod(s, &sz);
		res.insert(pair<CellType, double>(celltype, auxd));
		it++;
	}
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
		res.insert(pair<int, double>(sprtype, auxd));
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

	if (REPORT_OUT)
		cout << "Reading .vp file...\n";
	while (getline(fin, line))
		if (!(line.empty() || line.find_first_not_of(' ') == std::string::npos))
			if (line.at(0) != '#')
				inp.push_back(line);
	std::vector<std::string>::iterator it = inp.begin();
	//READ PARAMETERS IN ORDER FROM HERE:
	integration_mode = static_cast<int>(read_real_par(it));
	min_range_vertex_movement = read_real_par(it);
	max_range_vertex_movement = read_real_par(it);
	h = read_real_par(it);
	temperature_positive_energy = read_real_par(it);
	temperature_negative_energy = read_real_par(it);
	temperature_means_proportion_of_acceptance = read_real_par(it) > 0;

	energy_term1 = read_real_par(it);
	energy_term2 = read_real_par(it);
	energy_term3 = read_real_par(it);

	//spring_constant = read_real_par(it);
	t1_transition_critical_distance = read_real_par(it);
	length_rotated_edge = read_real_par(it);
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
	xcoord_controls_size = read_real_par(it) > 0;
	xcoord_decrease_exponent = read_real_par(it);

	line_tension = read_celltype_par(it, sz);
	line_tension_tissue_boundary = read_celltype_par(it, sz);
	perimeter_contract = read_celltype_par(it, sz);
	preferred_area_initial = read_celltype_par(it, sz);
	preferred_area_final = read_celltype_par(it, sz);
	division_angle_random_noise = read_celltype_par(it, sz);
	division_angle_longest_axis = read_celltype_par(it, sz);
	division_angle_external = read_celltype_par(it, sz);
	division_angle_external_degrees = read_celltype_par(it, sz);
	spring_type_constants = read_springtype_par(it, sz);
	max_cell_area = read_celltype_par(it, sz);
	cell_cycle_limit = read_celltype_par(it, sz);

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
}

/*
Initializes vertices from file
Input: ifstream pointing to a file defining vertex coordinates
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
	e.tension = spring_type_constants[0];

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
			e.tension = this->spring_type_constants[sprtype];
		}
		catch (const std::invalid_argument &ia)
		{
			std::cerr << "Invalid argument (No third spring column): " << ia.what() << '\n';
			e.tension = spring_type_constants[0];
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
		c.cell_cycle_state = start_cell_cycle_at_random ? std::rand() % static_cast<int>(cell_cycle_limit[c.type]) : 0;
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

		cells[c].perimeter_contractility = perimeter_contract[cells[c].type];
		cells[c].preferred_area = preferred_area_initial[cells[c].type];
		cells[c].area = calculateCellArea(cells[c]);
		cells[c].perimeter = calculateCellPerimeter(cells[c]);

		cells[c].division_angle_random_noise = division_angle_random_noise[cells[c].type];
		cells[c].division_angle_longest = division_angle_longest_axis[cells[c].type];
		cells[c].division_angle_external = division_angle_external[cells[c].type];
		cells[c].division_angle_external_degrees = division_angle_external_degrees[cells[c].type];

		cells[c].cell_cycle_limit = cell_cycle_limit[cells[c].type];
		cells[c].max_area = max_cell_area[cells[c].type];

		cells[c].vary_line_tension = vary_line_tension[cells[c].type];
		cells[c].edge_angle_prop_external = edge_angle_prop_external[cells[c].type];
		cells[c].edge_angle_prop_uniform = edge_angle_prop_uniform[cells[c].type];
		cells[c].edge_angle_prop_maxangle = edge_angle_prop_maxangle[cells[c].type];
		cells[c].edge_angle_prop_random = edge_angle_prop_random[cells[c].type];
		cells[c].edge_tension_external = edge_tension_external[cells[c].type];
		cells[c].edge_maxangle = edge_maxangle[cells[c].type];
		cells[c].edge_spatialmax_tension = edge_spatialmax_tension[cells[c].type];
		cells[c].edge_spatialmin_tension = edge_spatialmin_tension[cells[c].type];
	}

	for (int e = 0; e < edges.size(); e++)
	{
		setEdgeType(e);
		setEdgeTension(e);
	}
	for (int v = 0; v < vertices.size(); v++)
	{
		vertices[v].energy = calculateEnergy(vertices[v]);
	}
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
		edges[e].tension = line_tension_tissue_boundary[cells[aux].type];
		edges[e].can_transition = true;
	}
	else if (cells[c1].type == CellType::blade && cells[c2].type == CellType::blade)
	{
		edges[e].type = EdgeType::blade;
		edges[e].tension = line_tension[CellType::blade];
		edges[e].can_transition = true;
	}
	else if (cells[c1].type == CellType::hinge && cells[c2].type == CellType::hinge)
	{
		edges[e].type = EdgeType::hinge;
		edges[e].tension = line_tension[CellType::hinge];
		edges[e].can_transition = true;
	}
	else if ((cells[c1].type == CellType::hinge && cells[c2].type == CellType::blade) || (cells[c1].type == CellType::blade && cells[c2].type == CellType::hinge))
	{
		edges[e].type = EdgeType::hinge;
		edges[e].tension = 0.5 * (line_tension[CellType::hinge] + line_tension[CellType::blade]);
		edges[e].can_transition = false;
	}
	else if ((cells[c1].type == CellType::blade && cells[c2].type == CellType::vein_blade) ||
			 (cells[c2].type == CellType::blade && cells[c1].type == CellType::vein_blade) ||
			 (cells[c1].type == CellType::blade && cells[c2].type == CellType::vein_hinge) ||
			 (cells[c2].type == CellType::blade && cells[c1].type == CellType::vein_hinge))
	{
		edges[e].type = EdgeType::vein_blade;
		edges[e].tension = line_tension[CellType::vein_blade];
		edges[e].can_transition = true; //false; true
	}
	else if ((cells[c1].type == CellType::hinge && cells[c2].type == CellType::vein_hinge) ||
			 (cells[c2].type == CellType::hinge && cells[c1].type == CellType::vein_hinge) ||
			 (cells[c1].type == CellType::hinge && cells[c2].type == CellType::vein_blade) ||
			 (cells[c2].type == CellType::hinge && cells[c1].type == CellType::vein_blade))
	{
		edges[e].type = EdgeType::vein_hinge;
		edges[e].tension = line_tension[CellType::vein_hinge];
		edges[e].can_transition = true; //false; true
	}
	else if (cells[c1].type == CellType::vein_blade && cells[c2].type == CellType::vein_blade)
	{
		edges[e].type = EdgeType::vein;
		edges[e].tension = line_tension[CellType::blade]; //vein_blade is only for vein border
		edges[e].can_transition = true;					  //false;true
	}
	else if (cells[c1].type == CellType::vein_hinge && cells[c2].type == CellType::vein_hinge)
	{
		edges[e].type = EdgeType::vein;
		edges[e].tension = line_tension[CellType::hinge]; //vein_hinge is only for vein border
		edges[e].can_transition = true;					  //false;true
	}
	else if ((cells[c1].type == CellType::vein_blade && cells[c2].type == CellType::vein_hinge) || (cells[c1].type == CellType::vein_hinge && cells[c2].type == CellType::vein_blade))
	{
		edges[e].type = EdgeType::vein;
		edges[e].tension = 0.5 * (line_tension[CellType::hinge] + line_tension[CellType::blade]); //vein_hinge and vein_blade are only for vein border
		edges[e].can_transition = false;
	}
}

void Tissue::setEdgeTension(int e)
{
	float mins, maxs, maxangle, angle, pex, pmaxan, punif, prand, tensionext, tensionrand;
	int cellvar;
	if (edges[e].type == EdgeType::tissue_boundary)
	{
		/*cellvar = edges[e].cells[0] == EMPTY_CONNECTION ? edges[e].cells[1] : edges[e].cells[0];//find out which cell
		if(! cells[cellvar].vary_line_tension) return;
		float boundary_factor = edges[e].tension/(0.5*(cells[cellvar].edge_spatialmin_tension + cells[cellvar].edge_spatialmax_tension));

		mins = cells[cellvar].edge_spatialmin_tension*boundary_factor; //Minimal tension depending on angle
		maxs = cells[cellvar].edge_spatialmax_tension*boundary_factor; //Maximal tension depending on angle
		maxangle = cells[cellvar].edge_maxangle; //Angle of max tension (in degrees)
		tensionext = cells[cellvar].edge_tension_external; //Tension set from outside (gene expression etc)
		prand = cells[cellvar].edge_angle_prop_random;
		if(vary_edge_tension_with_time){ //If proportion determined by angle varies with time
			double time_factor = static_cast<double>(counter_moves_accepted)/max_accepted_movements;
			time_factor = expAdvance(time_factor, vary_edge_tension_time_exponent);
			double mint = edge_temporal_angle_efect_min[cells[cellvar].type];
			double maxt = edge_temporal_angle_efect_max[cells[cellvar].type];
			pex = 0.0;
			pmaxan = mint + (maxt - mint)*time_factor; //Proportion is determined by exponential function of time
			punif = 1.0 - pmaxan;

		}else{ //If proportion determined by angle does not vary with time
			pex = cells[cellvar].edge_angle_prop_external; //Proportions are determined by cell params directly
			pmaxan = cells[cellvar].edge_angle_prop_maxangle;
			punif = (cells[cellvar].edge_angle_prop_uniform);
		}//end if time determines influence of angle
	*/
		return;
	}
	else
	{ //not border
		if (!cells[edges[e].cells[0]].vary_line_tension && !cells[edges[e].cells[1]].vary_line_tension)
			return;
		if (cells[edges[e].cells[0]].vary_line_tension && cells[edges[e].cells[1]].vary_line_tension)
		{ //If both cells make edge tension vary
			mins = 0.5 * (cells[edges[e].cells[0]].edge_spatialmin_tension + cells[edges[e].cells[1]].edge_spatialmin_tension);
			maxs = 0.5 * (cells[edges[e].cells[0]].edge_spatialmax_tension + cells[edges[e].cells[1]].edge_spatialmax_tension);
			maxangle = 0.5 * (cells[edges[e].cells[0]].edge_maxangle + cells[edges[e].cells[1]].edge_maxangle);
			tensionext = 0.5 * (cells[edges[e].cells[0]].edge_tension_external + cells[edges[e].cells[1]].edge_tension_external); //This is not an angle! is tension value set from gene Expression directly

			if (vary_edge_tension_with_time)
			{ //If proportion determined by angle varies with time
				double time_factor = static_cast<double>(counter_moves_accepted) / max_accepted_movements;
				time_factor = expAdvance(time_factor, vary_edge_tension_time_exponent);
				double mint = (edge_temporal_angle_efect_min[cells[edges[e].cells[0]].type] + edge_temporal_angle_efect_min[cells[edges[e].cells[1]].type]) * 0.5;
				double maxt = (edge_temporal_angle_efect_max[cells[edges[e].cells[0]].type] + edge_temporal_angle_efect_max[cells[edges[e].cells[1]].type]) * 0.5;
				pex = 0.0;
				pmaxan = mint + (maxt - mint) * time_factor;
				prand = 0.5 * (cells[edges[e].cells[0]].edge_angle_prop_random + cells[edges[e].cells[1]].edge_angle_prop_random);
				punif = 1.0 - pmaxan; //prand only used in the end; this way it is easier to avoid negative proportions
			}
			else
			{ //If proportion determined by angle does not vary with time
				pex = 0.5 * (cells[edges[e].cells[0]].edge_angle_prop_external + cells[edges[e].cells[1]].edge_angle_prop_external);
				pmaxan = 0.5 * (cells[edges[e].cells[0]].edge_angle_prop_maxangle + cells[edges[e].cells[1]].edge_angle_prop_maxangle);
				punif = 0.5 * (cells[edges[e].cells[0]].edge_angle_prop_uniform + cells[edges[e].cells[1]].edge_angle_prop_uniform);
				prand = 0.5 * (cells[edges[e].cells[0]].edge_angle_prop_random + cells[edges[e].cells[1]].edge_angle_prop_random);
			} //end if time determines influence of angle
		}
		else
		{																								  //If only one cell makes edge tension vary
			cellvar = cells[edges[e].cells[0]].vary_line_tension ? edges[e].cells[0] : edges[e].cells[1]; //find out which cell
			mins = cells[cellvar].edge_spatialmin_tension;												  //Minimal tension depending on angle
			maxs = cells[cellvar].edge_spatialmax_tension;												  //Maximal tension depending on angle
			maxangle = cells[cellvar].edge_maxangle;													  //Angle of max tension (in degrees)
			tensionext = cells[cellvar].edge_tension_external;											  //Tension set from outside (gene expression etc)
			prand = cells[cellvar].edge_angle_prop_random;
			if (vary_edge_tension_with_time)
			{ //If proportion determined by angle varies with time
				double time_factor = static_cast<double>(counter_moves_accepted) / max_accepted_movements;
				time_factor = expAdvance(time_factor, vary_edge_tension_time_exponent);
				double mint = edge_temporal_angle_efect_min[cells[cellvar].type];
				double maxt = edge_temporal_angle_efect_max[cells[cellvar].type];
				pex = 0.0;
				pmaxan = mint + (maxt - mint) * time_factor; //Proportion is determined by exponential function of time
				punif = 1.0 - pmaxan;
			}
			else
			{												   //If proportion determined by angle does not vary with time
				pex = cells[cellvar].edge_angle_prop_external; //Proportions are determined by cell params directly
				pmaxan = cells[cellvar].edge_angle_prop_maxangle;
				punif = (cells[cellvar].edge_angle_prop_uniform);
			} //end if time determines influence of angle
		}	  //End if which cells determine edge variation
	}		  //End if is border
	maxangle *= M_PI / 180;
	angle = atan2(vertices[edges[e].vertices[1]].y - vertices[edges[e].vertices[0]].y, vertices[edges[e].vertices[1]].x - vertices[edges[e].vertices[0]].x); //angle of edge
	angle = abs(sin(0.5 * M_PI + abs(angle - maxangle)));																									 //Point of sin wave (from 0 to 1)
	angle = mins + angle * (maxs - mins);																													 //Value of tension at this point of sin wave

	//Now integrate with other factors influencing tension
	tensionrand = prand > 0 ? mins + (maxs - mins) * ((double)std::rand() / (RAND_MAX)) : 0;
	edges[e].tension = (punif * edges[e].tension + pmaxan * angle + pex * tensionext + prand * tensionrand) / (punif + pmaxan + pex + prand);
	//cout << "maxt: " << maxs << ", mint: " << mins <<", tension: " << edges[e].tension << ", angle: " << 180*atan2(vertices[edges[e].vertices[1]].y - vertices[edges[e].vertices[0]].y, vertices[edges[e].vertices[1]].x - vertices[edges[e].vertices[0]].x)/M_PI << ", angle tension: " << angle << ", random: " << tensionrand << ", prand: " << prand << endl;
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
	vertices[v].spring = EMPTY_CONNECTION;
	vertices[v].dead = false;
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
	double min, max, miny, maxy;
	bool set = false;
	for (Cell c : cells)
	{
		if (c.type == CellType::hinge || c.type == CellType::vein_hinge)
		{
			calculateCellCentroid(c);
			if (set)
			{
				if (c.centroid_x < min)
					min = c.centroid_x;
				else if (c.centroid_x > max)
					max = c.centroid_x;
				if (c.centroid_y < miny)
					miny = c.centroid_y;
				else if (c.centroid_y > maxy)
					maxy = c.centroid_y;
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
	}
	hinge_min_xpos = min;
	hinge_max_xpos = max;
	hinge_min_ypos = miny;
	hinge_max_ypos = maxy;
	//cout << "Min hinge position: " << hinge_min_xpos << "; Max hinge position: " << hinge_max_xpos << endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//Area of a polygon using Shoelace formula
double Tissue::calculateCellArea(const Cell &c)
{
	double area = 0.0;
	int previous = c.num_vertices - 1;
	for (int i = 0; i < c.num_vertices; i++)
	{
		area += (vertices[c.vertices[previous]].x * vertices[c.vertices[i]].y - vertices[c.vertices[previous]].y * vertices[c.vertices[i]].x);
		previous = i;
	}
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
		a = vertices[c.vertices[previous]].x * vertices[c.vertices[i]].y - vertices[c.vertices[previous]].y * vertices[c.vertices[i]].x;
		xc += a * (vertices[c.vertices[previous]].x + vertices[c.vertices[i]].x);
		yc += a * (vertices[c.vertices[previous]].y + vertices[c.vertices[i]].y);
		previous = i;
	}
	c.centroid_x = xc / (6 * c.area);
	c.centroid_y = yc / (6 * c.area);
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
inline double Tissue::calculateEnergy(Vertex &v)
{
	double term1 = 0, term2 = 0, term3 = 0;

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.cells[i] != EMPTY_CONNECTION)
			term1 += pow(cells[v.cells[i]].area / cells[v.cells[i]].preferred_area - 1, 2);
	}
	//Second term is edge length divided by preferred area of cell.
	//Most of the time both cells will be of the same type, but sometimes one will be hinge and the other will be blade
	//Therefore, I take the mean between the preferred areas of both cells,
	//except for Edges in the border (obviously)
	double pref_area;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.edges[i] != EMPTY_CONNECTION)
		{
			pref_area = edges[v.edges[i]].type == EdgeType::tissue_boundary ? edges[v.edges[i]].cells[0] == EMPTY_CONNECTION ? cells[edges[v.edges[i]].cells[1]].preferred_area
																															 : cells[edges[v.edges[i]].cells[0]].preferred_area
																			: (cells[edges[v.edges[i]].cells[0]].preferred_area + cells[edges[v.edges[i]].cells[1]].preferred_area) * 0.5;
			term2 += edges[v.edges[i]].tension * edges[v.edges[i]].length / sqrt(pref_area);
		}
	}

	if (v.spring != EMPTY_CONNECTION)
		term2 += springs[v.spring].tension * springs[v.spring].length / sqrt(pref_area);

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.cells[i] != EMPTY_CONNECTION)
			term3 += 0.5 * cells[v.cells[i]].perimeter_contractility * pow(cells[v.cells[i]].perimeter, 2) / cells[v.cells[i]].preferred_area;
	}

	return 0.5 * term1 * energy_term1 + term2 * energy_term2 + term3 * energy_term3;
}
//
void Tissue::moveVertex(Vertex &v, float x, float y)
{
	//Set new vertex coordinates. Original object is modified because v is passed by reference
	v.x = x;
	v.y = y;

	//Loop to re-calculate Edge lengths
	//Length of edge in index v.edges[i] is updated
	//Distance function is called with indexes of vertices touching edge (edge.vertices[0] and edge.vertices[1], where edge is edges[ v.edges[i] ])
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.edges[i] != EMPTY_CONNECTION)
		{
			this->edges[v.edges[i]].length = distance(this->edges[v.edges[i]].vertices[0], this->edges[v.edges[i]].vertices[1]);
			if (UPDATE_EDGE_TENSION_EVERY_MOVE)
			{
				setEdgeType(v.edges[i]);
				setEdgeTension(v.edges[i]);
			}
		}
	}
	if (v.spring != EMPTY_CONNECTION)
		this->springs[v.spring].length = distance(this->springs[v.spring].vertices[0], this->springs[v.spring].vertices[1]);
	// Loop to re-calculate Cell areas
	// calculateCellArea is called with the Cell structure located in vector this->cells, in position v.cells[i], as only argument
	// same for calculateCellPerimeter
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{ // re-calculate cell areas
		if (v.cells[i] != EMPTY_CONNECTION)
		{
			this->cells[v.cells[i]].area = calculateCellArea(this->cells[v.cells[i]]);
			this->cells[v.cells[i]].perimeter = calculateCellPerimeter(this->cells[v.cells[i]]);
		}
	}
	//v.energy = calculateEnergy(v);
}
//counter_favorable_accepted, counter_favorable_rejected, counter_unfav_accepted, counter_unfav_rejected
bool Tissue::tryMoveVertex()
{

	int vertex_to_move;
	do
	{
		vertex_to_move = std::rand() % static_cast<int>(vertices.size());
	} while (vertices[vertex_to_move].dead || !vertices[vertex_to_move].movable); //(!vertices[vertex_to_move].movable && STATIC_PRESENT || vertices[vertex_to_move].cells[0] == EMPTY_CONNECTION)); // A or static dead vertex cannot be selected

	double old_x = vertices[vertex_to_move].x;
	double old_y = vertices[vertex_to_move].y;
	double old_energy = calculateEnergy(vertices[vertex_to_move]); //Calculate again because it is not updated every time a cell area changes etc
	double angle, radius, new_x, new_y;
	double move_prob;
	bool cell_borders_cross = false;
	int cross_trials = 0;
	do
	{
		angle = 2 * M_PI * unif(generator);
		radius = unif(generator) * max_range_vertex_movement;
		new_x = old_x + cos(angle) * radius;
		new_y = old_y + sin(angle) * radius;
		moveVertex(vertices[vertex_to_move], new_x, new_y);
		vertices[vertex_to_move].energy = calculateEnergy(vertices[vertex_to_move]);

		if (temperature_means_proportion_of_acceptance)
		{ //Prob of accepting unfavourable movement ins constant
			move_prob = vertices[vertex_to_move].energy <= old_energy ? temperature_negative_energy : temperature_positive_energy;
		}
		else
		{ //Prob of accepting unfavourable movement depends on how much unfavourable the movement is
			move_prob = vertices[vertex_to_move].energy <= old_energy ? exp((vertices[vertex_to_move].energy - old_energy) / temperature_negative_energy) : exp(-(vertices[vertex_to_move].energy - old_energy) / temperature_positive_energy);
		}
		if (CHECK_EDGES_CROSS_AFTER_MOVE)
		{
			cell_borders_cross = check_if_edges_cross(vertex_to_move);
			if (cell_borders_cross)
			{
				moveVertex(vertices[vertex_to_move], old_x, old_y);
				vertices[vertex_to_move].energy = old_energy;
			}
			cross_trials++;
			if (cross_trials > MOVE_TRIALS)
				return false;
		}
	} while (cell_borders_cross);
	double move = unif(generator);

	if (move < move_prob)
	{
		detectChangesAfterMove(vertex_to_move);
		if (autonomous_cell_cycle)
		{ //NOTE: Maybe make this before detecting changes so the effect is immediate? Otherwise cells have to wait until next round
			advanceCellCycle(vertex_to_move);
		}
		if (time_controls_size && !xcoord_controls_size)
		{
			advanceSizeWithTime(vertex_to_move);
		}
		else if (time_controls_size && xcoord_controls_size)
		{
			advanceSizeWithXcoordAndTime(vertex_to_move);
		}
		else if (xcoord_controls_size)
		{
			advanceSizeWithXcoord(vertex_to_move);
		}
		//Counters
		if (vertices[vertex_to_move].energy <= old_energy)
		{
			counter_favorable_accepted++;
		}
		else
		{
			counter_unfav_accepted++;
		} //Counters
		return true;
	}
	else
	{
		//Counters
		if (vertices[vertex_to_move].energy <= old_energy)
		{
			counter_favorable_rejected++;
		}
		else
		{
			counter_unfav_rejected++;
		} //Counters
		moveVertex(vertices[vertex_to_move], old_x, old_y);
		return false;
	}
}

inline double Tissue::expAdvance(double x, float exponent)
{
	return exponent > 0 ? (1 - exp(-pow(x, exponent))) / EXP_FACTOR : ((1 - exp(-pow(x, exponent))) - EXP_FACTOR) / (1 - EXP_FACTOR);
}

void Tissue::advanceSizeWithXcoordAndTime(int vertex_moved)
{
	int caux;
	//double auxprint;
	double time_factor = static_cast<double>(counter_moves_accepted) / max_accepted_movements;
	time_factor = expAdvance(time_factor, time_decrease_exponent); //(1 - exp(-pow(time_factor, time_decrease_exponent)))/EXP_FACTOR;
	float pos_factor = 0, aux_area;

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		cells[caux].preferred_area = preferred_area_initial[cells[caux].type] + (preferred_area_final[cells[caux].type] - preferred_area_initial[cells[caux].type]) * time_factor; //Caution: problems if step_mode is on
		//auxprint = cells[caux].preferred_area;
		if (cells[caux].type == CellType::vein_hinge)
		{
			calculateCellCentroid(cells[caux]);
			pos_factor = cells[caux].centroid_x - hinge_min_xpos;
			pos_factor = pos_factor < 0 ? 0 : pos_factor / (hinge_max_xpos - hinge_min_xpos);
			aux_area = preferred_area_initial[CellType::vein_blade] + (preferred_area_final[CellType::vein_blade] - preferred_area_initial[CellType::vein_blade]) * time_factor;
			cells[caux].preferred_area += (aux_area - cells[caux].preferred_area) * expAdvance(pos_factor, xcoord_decrease_exponent);
		}
		else if (cells[caux].type == CellType::hinge)
		{
			calculateCellCentroid(cells[caux]);
			pos_factor = cells[caux].centroid_x - hinge_min_xpos;
			pos_factor = pos_factor < 0 ? 0 : pos_factor / (hinge_max_xpos - hinge_min_xpos);
			aux_area = preferred_area_initial[CellType::blade] + (preferred_area_final[CellType::blade] - preferred_area_initial[CellType::blade]) * time_factor;
			cells[caux].preferred_area += (aux_area - cells[caux].preferred_area) * expAdvance(pos_factor, xcoord_decrease_exponent);
		}
		if (keep_area_after_division)
		{
			cells[caux].preferred_area /= pow(2, cells[caux].num_divisions);
		}
		//cout << endl;
	}
} //advanceSizeWithTimeAndXcoord

void Tissue::advanceSizeWithXcoord(int vertex_moved)
{
	int caux;
	float pos_factor = 0;

	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		if (cells[caux].type == CellType::vein_hinge)
		{
			calculateCellCentroid(cells[caux]);
			pos_factor = cells[caux].centroid_x - hinge_min_xpos;
			pos_factor = pos_factor < 0 ? 0 : pos_factor / (hinge_max_xpos - hinge_min_xpos);
			cells[caux].preferred_area = preferred_area_initial[CellType::vein_hinge] + (preferred_area_initial[CellType::vein_blade] - preferred_area_initial[CellType::vein_hinge]) * expAdvance(pos_factor, xcoord_decrease_exponent); //(1 - exp(-pow(pos_factor, xcoord_decrease_exponent)))/EXP_FACTOR;
			if (keep_area_after_division)
			{
				cells[caux].preferred_area /= pow(2, cells[caux].num_divisions);
			}
			//cout << "HingeVein cell: " << " x abs = " << cells[caux].centroid_x << ", x rel = " << pos_factor << ", pref area = " << cells[caux].preferred_area << endl;
		}
		else if (cells[caux].type == CellType::hinge)
		{
			calculateCellCentroid(cells[caux]);
			pos_factor = cells[caux].centroid_x - hinge_min_xpos;
			pos_factor = pos_factor < 0 ? 0 : pos_factor / (hinge_max_xpos - hinge_min_xpos);
			cells[caux].preferred_area = preferred_area_initial[CellType::hinge] + (preferred_area_initial[CellType::blade] - preferred_area_initial[CellType::hinge]) * expAdvance(pos_factor, xcoord_decrease_exponent); //(1 - exp(-pow(pos_factor, xcoord_decrease_exponent)))/EXP_FACTOR;
			//cout << "Hinge cell: " << " x abs = " << cells[caux].centroid_x << ", x rel = " << pos_factor << ", pref area = " << cells[caux].preferred_area << endl;
			if (keep_area_after_division)
			{
				cells[caux].preferred_area /= pow(2, cells[caux].num_divisions);
			}
		}
	}
} //advanceSizeWithXcoord

void Tissue::advanceCellCycle(int vertex_moved)
{
	int caux;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		cells[caux].cell_cycle_state += 1;
		if (cells[caux].cell_cycle_state >= cells[caux].cell_cycle_limit)
		{
			cells[caux].can_divide = true;
		}
		if (cell_cycle_controls_size && !(keep_area_after_division && (cells[caux].type == CellType::blade || cells[caux].type == CellType::vein_blade)))
		{
			cells[caux].preferred_area = preferred_area_initial[cells[caux].type] + (preferred_area_final[cells[caux].type] - preferred_area_initial[cells[caux].type]) * cells[caux].cell_cycle_state / cells[caux].cell_cycle_limit;
		}
	}
} //advanceCellCycle

void Tissue::advanceSizeWithTime(int vertex_moved)
{
	int caux;
	double time_factor = static_cast<double>(counter_moves_accepted) / max_accepted_movements;
	time_factor = expAdvance(time_factor, time_decrease_exponent); //(1 - exp(-pow(time_factor, time_decrease_exponent)))/EXP_FACTOR;
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		caux = vertices[vertex_moved].cells[i];
		if (caux == EMPTY_CONNECTION)
			continue;
		cells[caux].preferred_area = preferred_area_initial[cells[caux].type] + (preferred_area_final[cells[caux].type] - preferred_area_initial[cells[caux].type]) * time_factor; //Caution: problems if step_mode is on
		//cout << "cell="<<caux <<"; type="<<static_cast<int>(cells[caux].type) << "; time factor= " << time_factor << ", init=" << preferred_area_initial[cells[caux].type] << "; final=" << preferred_area_final[cells[caux].type] << "; area=" << cells[caux].preferred_area << endl;
		if (keep_area_after_division)
		{
			cells[caux].preferred_area /= pow(2, cells[caux].num_divisions);
		}
	}
} //advanceSizeWithTime

void Tissue::produceOutputs(std::string add_to_name)
{

std:
	string fname = simname + "_" + add_to_name + "_" + std::to_string(int(counter_moves_accepted / write_every_N_moves));
	writeCellsFile(fname);
	writePointsFile(fname);
	writeEdgeDataTable(fname);
	if (num_springs > 0)
		writeSpringsFile(fname);
	//writeAllData(fname); //THIS IS USEFUL TO DEBUG, BUT THE FORMAT IS NOT READ BY PLOTTING PROGRAM
	writeCellDataTable(fname);
	if (REPORT_OUT)
		cout << "\nWritting file: " << int(counter_moves_accepted / write_every_N_moves) << " at move " << counter_moves_accepted << endl;
		std::cout << getStats() << endl;
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
		if (v.cells[i] == EMPTY_CONNECTION)
			continue;
		c = v.cells[i];
        aux = (cells[c].area/cells[c].preferred_area- 1)/cells[c].preferred_area;
		aux_this = which(v.ind, cells[c].vertices, cells[c].num_vertices); //These three lines are repeated in TERM3 but since TERM3 is probably going to be optimized it doesnt matter
		//IMPORTANT: ASSUMES THAT VERTICES IN CELLS ARE ORDERED CLOCK-WISE
		//Calculate derivative of area according to Shoelace formula from Wikipedia
		aux_next = cells[c].vertices[(aux_this + cells[c].num_vertices - 1)%cells[c].num_vertices ];
		aux_prev = cells[c].vertices[(aux_this + 1)%cells[c].num_vertices ];
		t1x += aux*(vertices[aux_prev].y - vertices[aux_next].y);
		t1y += aux*(vertices[aux_next].x - vertices[aux_prev].x);
		cout << "    Term 1: aux " << aux << ", pref area " << cells[c].preferred_area << ", area " << cells[c].area << ", cell " << c << endl;
	
	} //Term1
	//TERM 2
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.edges[i] == EMPTY_CONNECTION)
			continue;
		c = v.edges[i];	
		pref_area = edges[c].type == EdgeType::tissue_boundary ? edges[c].cells[0] == EMPTY_CONNECTION ? cells[edges[c].cells[1]].preferred_area
																											 : cells[edges[c].cells[0]].preferred_area
																		: (cells[edges[c].cells[0]].preferred_area + cells[edges[c].cells[1]].preferred_area) * 0.5;
		aux_next = edges[c].vertices[0] == v.ind ? edges[c].vertices[1] : edges[c].vertices[0];
		if(edges[c].length > NUMERIC_THRESHOLD){//Control for division by 0
			aux = edges[c].tension / (edges[c].length * pref_area);
		}else{
			aux = edges[c].tension / (NUMERIC_THRESHOLD * pref_area);
		}
		t2x += (v.x - vertices[aux_next].x) * aux;
		t2y += (v.y - vertices[aux_next].y) * aux;
		cout << "    Term 2: aux " << aux << ", x - next: " << v.x - vertices[aux_next].x << ", y - next: " << v.y - vertices[aux_next].y << ", len: " << springs[v.spring].length << ", tens: " << springs[v.spring].tension  << ", edge " << c  << endl;	
	} //Term 2
	if (v.spring != EMPTY_CONNECTION){
		if(springs[v.spring].length > NUMERIC_THRESHOLD){//Control for division by 0
			aux = springs[v.spring].tension / (springs[v.spring].length * pref_area);
		}else{
			aux = springs[v.spring].tension / (NUMERIC_THRESHOLD * pref_area);
		}
		aux_next = springs[v.spring].vertices[0] == v.ind ? springs[v.spring].vertices[1] : springs[v.spring].vertices[0];
		t2x += (v.x - vertices[aux_next].x) * aux; // Here may be a numeric problem when distance between 2 points is already 0
		t2y += (v.y - vertices[aux_next].y) * aux;
		cout << "    Term 2 spring: aux" << aux << ", x - next: " << v.x - vertices[aux_next].x << ", y - next: " << v.y - vertices[aux_next].y << ", len: " << springs[v.spring].length << ", tens: " << springs[v.spring].tension  << endl;
	}//Spring part of term2
	//Term 3 THIS TERM MUST BE OPTIMIZED
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (v.cells[i] == EMPTY_CONNECTION)
			continue;
		c = v.cells[i];
		aux_this = which(v.ind, cells[c].vertices, cells[c].num_vertices);
		aux_prev = cells[c].vertices[(aux_this + cells[c].num_vertices - 1)%cells[c].num_vertices ];
		aux_next = cells[c].vertices[(aux_this + 1)%cells[c].num_vertices ];
		aux = cells[c].perimeter*cells[c].perimeter_contractility/cells[c].preferred_area;
        auxd1 = distance(v.ind, aux_prev);
        auxd2 = distance(v.ind, aux_next); //Maybe numeric problem
		auxd1 = auxd1 > NUMERIC_THRESHOLD ? auxd1 : NUMERIC_THRESHOLD;
		auxd2 = auxd2 > NUMERIC_THRESHOLD ? auxd2 : NUMERIC_THRESHOLD;		
		t3x += aux*( (v.x - vertices[aux_prev].x)/auxd1 + (v.x - vertices[aux_next].x)/auxd2);
		t3y += aux*( (v.y - vertices[aux_prev].y)/auxd1 + (v.y - vertices[aux_next].y)/auxd2);
		cout << "    Term 3 aux " << aux << ", d1: " << auxd1 << ",d2 " << auxd2 << endl;	
	} //Term3
	cout << "v.ind=" << v.ind << ": x1=" << t1x << ", y1=" << t1y << ", x2=" << t2x << ", y2=" << t2y << ", x3=" << t3x << ", y3=" << t3y;

	t1x = isnan(t1x) || isinf(t1x) ? 0 : t1x;
	t1y = isnan(t1y) || isinf(t1y) ? 0 : t1y;
	t2x = isnan(t2x) || isinf(t2x) ? 0 : t2x;
	t2y = isnan(t2y) || isinf(t2y) ? 0 : t2y;
	t3x = isnan(t3x) || isinf(t3x) ? 0 : t3x;
	t3y = isnan(t3y) || isinf(t3y) ? 0 : t3y;
	
	pd.x = 0.5 * t1x * energy_term1 + t2x * energy_term2 + t3x * energy_term3;
	pd.y = 0.5 * t1y * energy_term1 + t2y * energy_term2 + t3y * energy_term3;
	if (max_range_vertex_movement > 0 && (abs(pd.x) > max_range_vertex_movement || abs(pd.y) > max_range_vertex_movement))
	{
		pd.x = pd.x / (pd.x + pd.y) * max_range_vertex_movement;
		pd.y = pd.y / (pd.x + pd.y) * max_range_vertex_movement;
	}
	//Add some noise
	if(temperature_positive_energy > 0){
		double angle = 2 * M_PI * unif(generator);
		double radius = unif(generator) * temperature_positive_energy;
		pd.x += cos(angle) * radius;
		pd.y += sin(angle) * radius;
	}
	cout << ", dx=" << pd.x << ", dy=" << pd.y << endl << endl;
}//derivativeVertexPos used in Simulate Euler
void Tissue::simulateMonteCarlo()
{
	cout << "SIMULATING with MONTE CARLO method." << endl;
	if (!step_mode)
	{
		produceOutputs();
	}
	//MAIN SIMULATION LOOP
	do
	{
		if (tryMoveVertex())
		{ //Tries a random movement; if accepted:
			counter_moves_accepted++;
			performRearrangements(); // Performs any transition/rearrangement needed
			if (counter_moves_accepted % write_every_N_moves == 0 && !step_mode)
				produceOutputs(); //&& counter_moves_accepted > 70000
		}						  //End if move accepted
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
	this->generator = generator;
	this->unif = unif;
	switch (this->integration_mode)
	{
	case INTEGR_MONTECARLO:
		simulateMonteCarlo();
		break;
	case INTEGR_EULER:
		simulateEuler();
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
	for (int i = 0; i < CELLS_PER_VERTEX; i++)
	{
		if (vertices[vertex_moved].edges[i] == EMPTY_CONNECTION)
			continue;
		ee = &edges[vertices[vertex_moved].edges[i]];
		Vertex *v1 = &vertices[ee->vertices[0]];
		Vertex *v2 = &vertices[ee->vertices[1]];
		if (ee->length <= t1_transition_critical_distance && v1->movable && v2->movable && ee->can_transition)
		{ //if length < critical distance
			if (ee->cells[0] != EMPTY_CONNECTION && ee->cells[1] != EMPTY_CONNECTION)
			{ //if the edge touches two cells, and there are 4 cells involved, then T1
				int cellnum = 0;
				for (int j = 0; j < CELLS_PER_VERTEX; j++)
					if (v1->cells[j] != EMPTY_CONNECTION)
						cellnum++;
				for (int j = 0; j < CELLS_PER_VERTEX; j++)
					if (v2->cells[j] != EMPTY_CONNECTION && !contains(v2->cells[j], v1->cells, CELLS_PER_VERTEX))
						cellnum++;
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
				for (int j = 0; j < CELLS_PER_VERTEX; j++)
					if (v2->edges[j] != EMPTY_CONNECTION)
						edgenum++;
				for (int j = 0; j < CELLS_PER_VERTEX; j++)
					if (v1->edges[j] != EMPTY_CONNECTION && !contains(v1->edges[j], v2->edges, CELLS_PER_VERTEX))
						edgenum++;
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
		else if (cells[caux].area >= max_cell_area[cells[caux].type] && cells[caux].can_divide)
		{
			rearrangements_needed.push(Rearrangement{caux, RearrangementType::divide_cell});
		}
	}

} //detect rearrangements

void Tissue::performRearrangements()
{

	Rearrangement r;
	while (!rearrangements_needed.empty())
	{
		r = rearrangements_needed.front();
		rearrangements_needed.pop();
		switch (r.type)
		{
		case RearrangementType::t1:
			//cout << "Entering t1 \n";
			if (T1_ACTIVE)
				make_t1(r);
			break;
		case RearrangementType::t2:
			//cout << "Entering t2 \n";
			if (T2_ACTIVE)
				make_t2(r);
			break;
		case RearrangementType::divide_cell:
			//cout << "Entering division \n";
			if (DIVISION_ACTIVE)
				make_divide_cell(r);
			break;
		case RearrangementType::divide_edge:
			break;
		case RearrangementType::join_limit_edges:
			//cout << "Entering join edges\n";
			if (JOIN_EDGES_ACTIVE)
				make_join_limit_edges(r);
			break;
		case RearrangementType::t1_at_border_outwards:
			//cout << "Entering t1 outwards\n";
			if (T1_BORDER_OUTWARDS_ACTIVE)
				make_t1_at_border_outwards(r);
			break;
		case RearrangementType::t1_at_border_inwards:
			//cout << "Entering t1 inwards " << counter_move_trials << "\n";
			if (T1_BORDER_INWARDS_ACTIVE)
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

	if (REPORT_OUT)
		cout << "T1 transition inwards: v1=" << v1->ind << ", v2=" << v2->ind << "; mov. accepted: " << this->counter_moves_accepted << "; T1: " << counter_t1 << "; T1 abortions: " << counter_t1_abortions << ", T1 in: " << counter_t1_inwards << endl;

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

	double dist1 = t1_get_dist_sum(v1, v2, cc1, cc2, sp1, sp1); //Not very efficient: to check for edge collisions will examine twice sp1
	double dist2 = t1_get_dist_sum(v1, v2, cc2, cc1, sp1, sp1);

	//cout << "E" << endl;
	if (dist2 > 0 && (dist2 < dist1 || dist1 < 0))
	{
		//cout << "E2" << endl;
		int aux = common_cell1;
		common_cell1 = common_cell2;
		common_cell2 = aux;
		cc1 = &this->cells[common_cell1];
		cc2 = &this->cells[common_cell2];
	}
	else if (dist2 < 0 && dist1 < 0)
	{ //edges always cross
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
	if (REPORT_T1)
	{
		ff.open("T1_REPORT_outwards_" + std::to_string(counter_moves_accepted) + "_" + std::to_string(v1->ind) + "_" + std::to_string(v2->ind) + ".txt");

		ff << "\nBEFORE\n\n";
		ff << "edge: " << edge << "\n\n";
		ff << VERTEX_HEADER << *v1 << "\n"
		   << *v2 << "\n\n";
		ff << CELL_HEADER << *cc1 << "\tcommon_cell1\n"
		   << *cc2 << "\tcommon_cell2\n"
		   << *sp1 << "\tsp1\n\n";
		ff << getStats();

		ff << "\n\nDist when cell " << cc1->ind << " goes with vertex " << v1->ind << ": " << dist1 << "\n";
		ff << "Dist when cell " << cc1->ind << " goes with vertex " << v2->ind << ": " << dist2 << "\n";
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
	edges[edge].tension = line_tension_tissue_boundary[sp1->type];
	edges[edge].type = EdgeType::tissue_boundary;
	this->counter_t1_outwards++;

	if (REPORT_OUT)
		cout << "T1 transition outwards: v1=" << v1->ind << ", v2=" << v2->ind << "; mov. accepted: " << this->counter_moves_accepted << "; T1: " << counter_t1 << "; T1 abortions: " << counter_t1_abortions << ", T1 out: " << counter_t1_outwards << endl;

	if (REPORT_T1)
	{
		ff << "\n\nAFTER\n\n";
		ff << VERTEX_HEADER << *v1 << "\n"
		   << *v2 << "\n\n";
		ff << CELL_HEADER << *cc1 << "\tcommon_cell1\n"
		   << *cc2 << "\tcommon_cell2\n"
		   << *sp1 << "\tsp1\n\n";
		ff << getStats();
		ff.close();
	}

} /// End make_t1_outwards

void Tissue::make_divide_cell(Rearrangement &r)
{
	int cell = r.element_index;
	if (cells[cell].area < cells[cell].max_area || cells[cell].num_vertices < 3 || cells[cell].dead || !cells[cell].can_divide)
		return;

	double x1, x2, y1, y2;
	int e1 = -1, e2 = -1;
	if (!getDivisionPoints(cell, x1, x2, y1, y2, e1, e2))
	{
		return;
	}

	int newvind1 = newVertex(x1, y1);
	int newvind2 = newVertex(x2, y2); //create new vertices that are going to be positioned at (x1, y1) and (x2, y2) insideedges e1 and e2
	int newcind = newCell();
	cells[newcind].type = cells[cell].type;
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
	cells[newcind].cell_cycle_limit = cells[cell].cell_cycle_limit;

	cells[newcind].max_area = cells[cell].max_area;

	if (this->autonomous_cell_cycle)
	{
		//cells[newcind].cell_cycle_state = 0; //Already done in newCell()
		cells[cell].cell_cycle_state = 0;
		//cells[newcind].can_divide = 0; //Already done in newCell()
		cells[cell].can_divide = false;
	}

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

	cells[newcind].area = calculateCellArea(cells[newcind]);
	cells[newcind].perimeter = calculateCellPerimeter(cells[newcind]);
	cells[cell].area = calculateCellArea(cells[cell]);
	cells[cell].perimeter = calculateCellPerimeter(cells[cell]);

	cells[cell].num_divisions++;
	cells[newcind].num_divisions = cells[cell].num_divisions;

	if (keep_area_after_division)
	{
		//cells[newcind].preferred_area = cells[newcind].area;
		//cells[cell].preferred_area = cells[cell].area;
		cells[cell].preferred_area /= pow(2, cells[cell].num_divisions);
		cells[newcind].preferred_area = cells[cell].preferred_area;
	}
	else if (this->cell_cycle_controls_size)
	{
		cells[newcind].preferred_area = this->preferred_area_initial[cells[cell].type];
		cells[cell].preferred_area = this->preferred_area_initial[cells[cell].type];
	}
	else
	{
		cells[newcind].preferred_area = cells[cell].preferred_area;
	}

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
	if (REPORT_OUT)
		cout << ">DIVISION: moves accepted= " << counter_moves_accepted << "; divi. accepted= " << counter_divisions << "; Cell= " << cell << "; New cell=" << newcind << "; new vertex 1= " << newvind1 << "; new vertex 2=" << newvind2 << "; centroid_x= " << 0.5 * (vertices[newvind1].x + vertices[newvind2].x) << "; centroid_y= " << 0.5 * (vertices[newvind1].y + vertices[newvind2].y) << "; cell_type= " << static_cast<int>(cells[cell].type) << endl;

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
	if (REPORT_OUT && (e1 < 0 || e2 < 0))
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
	if (edges[edge].length > t1_transition_critical_distance)
		return; //Check that condition is still true (other rearrangements could have taken place since detection)

	Vertex *v1 = &vertices[edges[edge].vertices[0]];
	Vertex *v2 = &vertices[edges[edge].vertices[1]];

	if (contains(EMPTY_CONNECTION, v1->cells, CELLS_PER_VERTEX) || contains(EMPTY_CONNECTION, v2->cells, CELLS_PER_VERTEX) || contains(EMPTY_CONNECTION, v1->edges, CELLS_PER_VERTEX) || contains(EMPTY_CONNECTION, v2->edges, CELLS_PER_VERTEX))
		return;

	double old_x1 = v1->x, old_x2 = v2->x, old_y1 = v1->y, old_y2 = v2->y, old_length = edges[edge].length;

	//Get the 2 cells that are common to both vertices and cell that is exclusive to v1
	int common_cell1 = -1, common_cell2 = -1, only_v1 = -1, only_v2 = -1; //, only_v1_indv1 = -1, only_v2_indv2 = -1;
	t1_getCellRelationships(v1, v2, common_cell1, common_cell2, only_v1, only_v2);
	Cell *cc1 = &this->cells[common_cell1]; //These pointers to Cell structures will be used later
	Cell *cc2 = &this->cells[common_cell2];
	Cell *sp1 = &this->cells[only_v1];
	Cell *sp2 = &this->cells[only_v2];

	if (cc1->num_vertices < 4 || cc2->num_vertices < 4)
		return; //One of the cells can't lose any other vertex
	//t1_rotate_edge(v1, v2, edge, cc1, cc2, sp1, sp2);
	t1_rotate_edge90degrees(v1, v2, edge);

	//Determine which cell goes with which point. Compare new point locations of v1 and v2 with all the neighbours of v1 and v2.

	double dist1 = t1_get_dist_sum(v1, v2, cc1, cc2, sp1, sp2);
	double dist2 = t1_get_dist_sum(v1, v2, cc2, cc1, sp1, sp2);

	if (dist2 > 0 && (dist2 < dist1 || dist1 < 0))
	{
		int aux = common_cell1;
		common_cell1 = common_cell2;
		common_cell2 = aux;
		cc1 = &this->cells[common_cell1];
		cc2 = &this->cells[common_cell2];
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
		ff.open("T1_REPORT_" + std::to_string(counter_moves_accepted) + "_" + std::to_string(v1->ind) + "_" + std::to_string(v2->ind) + ".txt");

		ff << "\nBEFORE\n\n";
		ff << "edge: " << edge << "\n\n";
		ff << VERTEX_HEADER << *v1 << "\n"
		   << *v2 << "\n\n";
		ff << CELL_HEADER << *cc1 << "\tcommon_cell1\n"
		   << *cc2 << "\tcommon_cell2\n"
		   << *sp1 << "\tsp1\n"
		   << *sp2 << "\tsp2\n\n";
		ff << getStats();

		ff << "\n\nDist when cell " << cc1->ind << " goes with vertex " << v1->ind << ": " << dist1 << "\n";
		ff << "Dist when cell " << cc1->ind << " goes with vertex " << v2->ind << ": " << dist2 << "\n";
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

	if (REPORT_OUT)
		cout << ">T1 transition: v1=" << v1->ind << "; v2=" << v2->ind << "; mov. accepted=" << this->counter_moves_accepted << "; T1=" << counter_t1 << "; T1 abortions=" << counter_t1_abortions << "; centroid_x=" << 0.5 * (v1->x + v2->x) << "; centroid_y=" << 0.5 * (v1->y + v2->y) << "; edge type=" << static_cast<int>(edges[edge].type) << endl;
	if (REPORT_T1)
	{
		ff << "\n\nAFTER\n\n";
		ff << VERTEX_HEADER << *v1 << "\n"
		   << *v2 << "\n\n";
		ff << CELL_HEADER << *cc1 << "\tcommon_cell1\n"
		   << *cc2 << "\tcommon_cell2\n"
		   << *sp1 << "\tsp1\n"
		   << *sp2 << "\tsp2\n\n";
		ff << getStats();
		ff.close();
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
				if (REPORT_OUT)
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

	StraightLine lv1_c1n1, lv1_c1n2, lv2_c2n1, lv2_c2n2, lv1v2, cell_edge;
	lv1_c1n1 = getLineFromEdge(v1, &this->vertices[v1_n1]);
	lv1_c1n2 = getLineFromEdge(v1, &this->vertices[v1_n2]);
	lv2_c2n1 = getLineFromEdge(v2, &this->vertices[v2_n1]);
	lv2_c2n2 = getLineFromEdge(v2, &this->vertices[v2_n2]);
	lv1v2 = getLineFromEdge(v1, v2);
	vector<StraightLine> lines = {lv1_c1n1, lv1_c1n2, lv2_c2n1, lv2_c2n2, lv1v2};
	vector<Cell *> neighbour_cells = {c1, c2, s1, s2};
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
				if (!(contains(ee->vertices[0], indices_avoid, 6) && contains(ee->vertices[1], indices_avoid, 6)) &&
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
			int dead_vert = springs[v2->spring].vertices[0] == v2->ind ? springs[v2->spring].vertices[1] : springs[v2->spring].vertices[0];
			int other_vert = springs[v1->spring].vertices[0] == v1->ind ? springs[v1->spring].vertices[1] : springs[v1->spring].vertices[0];
			vertices[other_vert].x = (vertices[other_vert].x + vertices[dead_vert].x) * 0.5;
			vertices[other_vert].y = (vertices[other_vert].y + vertices[dead_vert].y) * 0.5;
			springs[v1->spring].length = distance(springs[v1->spring].vertices[0], springs[v1->spring].vertices[1]);
			vertices[dead_vert].dead = true;
			dead_vertices.push(dead_vert);
			this->num_vertices--;
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

	if (CONTROL_CELLS_2SIDES)
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
						if (REPORT_OUT)
							cout << "  in T2: Killing cell of size 2: " << endl
								 << cells[neicell] << endl;
						killCellWith2Vertices(neicell);
					}
				}
			}
		}
	}
	if (REPORT_OUT)
		cout << ">T2 transition: cell=" << cell << "; survivor vertex =" << v1 << "; v2= " << v2 << "; v3= " << v3 << "; v2nei= " << v2nei << "; v3nei= " << v3nei << "; mov. accepted= " << this->counter_moves_accepted << "; T1= " << counter_t1 << "; T1 abortions=" << counter_t1_abortions << "; T2= " << counter_t2 << "; centroid_x= " << vertices[v1].x << "; centroid_y= " << vertices[v1].y << "; cell_type= " << static_cast<int>(cells[cell].type) << endl;
	//cout << "K\n";
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
	if (REPORT_OUT)
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
	return ((x_intersect <= a.x1 && x_intersect >= a.x2) || (x_intersect <= a.x2 && x_intersect >= a.x1)) && ((x_intersect <= b.x1 && x_intersect >= b.x2) || (x_intersect <= b.x2 && x_intersect >= b.x1));
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
}

void Tissue::setStaticVertex(int v)
{
	float xprop = (vertices[v].x - hinge_min_xpos) / (hinge_max_xpos - hinge_min_xpos);
	float yprop = (vertices[v].y - hinge_min_ypos) / (hinge_max_ypos - hinge_min_ypos);
	vertices[v].movable = !(yprop > 0.5 && xprop < add_static_to_hinge); //Assumes that blade is always to the right of hinge
}

bool Tissue::AddSpringToVertex(int v, float minx, float maxx)
{
	float pos = (vertices[v].x - minx) / (maxx - minx);
	for (auto minpos = spring_type_min_positions.rbegin(); minpos != spring_type_min_positions.rend(); ++minpos)
	{
		if (pos > minpos->second)
		{
			Edge e;
			e.dead = false;
			e.type = EdgeType::spring;
			e.cells[0] = EMPTY_CONNECTION;
			e.cells[1] = EMPTY_CONNECTION;
			e.tension = spring_type_constants[minpos->first];
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

void Tissue::emptyDivisions()
{
	divisionrecord_q aux;
	std::swap(aux, this->past_divisions);
}

std::string Tissue::getStats()
{
	if (!REPORT_OUT)
		return "-";
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
			for (auto const &springtype : spring_type_constants)
			{
				if (abs(springtype.second - s.tension) < 0.001)
				{
					fo << springtype.first;
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
			fo << v.x << "\t" << v.y << "\t" << v.ind << "\t" << v.movable << "\n";
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

///////////////////////////////////////////////////////////////////////////////////////////////
//What follows is only useful to print easily, not part of the model
//
//Overloading of << operator for VERTICES. Useful to print
std::ostream &operator<<(std::ostream &out, const Vertex &v)
{
	out << v.ind << "\t" << v.x << "\t" << v.y << "\t" << v.energy << "\t" << v.movable << "\t" << v.spring << "\t";
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
	out << c.ind << "\t" << int(c.type) << "\t" << c.area << "\t" << c.preferred_area << "\t" << c.perimeter << "\t" << c.perimeter_contractility << "\t";
	out << c.centroid_x << "\t" << c.centroid_y << "\t";
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
	return out;
}

//Overloading of << operator for EDGES. Useful to print
std::ostream &operator<<(std::ostream &out, const Edge &e)
{
	out << e.ind << "\t" << int(e.type) << "\t" << e.length << "\t" << e.tension << "\t";
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
