#ifndef __BASICGRN_H
#define __BASICGRN_H

#include "VertexSystem.h"
#include "GXMatrix.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <cstdlib>
#include <random>

using namespace std;
//Files with this extension should contain:
// - Type of each gene (including cell, edge and vertex properties, which for practical purposes count as genes)
// - Parameters for each gene for each cell type (degradation rate, diffusion rate, initial expression and constant expression)
const std::string INTERACTIONS_EXTENSION = ".grn";
const std::string EXPRESSION_EXTENSION = ".expr";
const std::string OUTPUT_EXTENSION = ".grn2";

const int NOT_PRESENT = -999;

//This should be updated in other implementations
enum class GeneType{intracel = 0, diffusible = 1, cell_property = 2, edge_property = 3, vertex_property = 4};

//A map in which keys are cell types and values are vectors containing the parameter for each gene in the key cell type
typedef std::map<CellType, std::vector<double> > ctparams;

//A map in which keys are cell types and values are matrices containing the strength for each interaction in the key cell type
typedef std::map<CellType, GXMatrix<double> > ctinteractions;

//This type is used to avoid looping through all possible regulators each time (more efficient if connections are sparse, which is expected)
typedef std::vector< std::vector<int> >list_of_regulators; 
typedef std::map<CellType, list_of_regulators > map_of_regulators;


struct Gene_params{ //This should be updated in other implementations
      ctparams km, degr, diff_rate, initial_expr, constant_expr, angle_influence;
};

class basicGRN {
  private:
    double h, time;
    float grn_time_per_step; //each iteration
    int vertex_moves_per_step; //each iterat
    int final_moves_accepted;
    int write_every_N_moves; 
    int num_genes, num_cells;
    std::string name;
    GXMatrix<double> expression;
    ctinteractions interactions;
    std::vector<GXMatrix<double>> rungekutta_parts; //Store k1, k2, k3, k4
    std::vector<GeneType> gene_types;
    std::vector<bool> gene_takes_only_positive_values;
    Gene_params params;
    map_of_regulators regulators;
    Tissue *cell_grid;

    void readGeneTypes(std::vector<std::string>::iterator& it);
    void readGenesThatTakeOnlyPositiveValues(std::vector<std::string>::iterator& it);
    void readParams(std::vector<std::string>::iterator& it);
    ctparams readSingleParam(std::vector<std::string>::iterator& it);
    void readGRN(std::vector<std::string>::iterator& it);
    void get_current_cell_grid_params();
    void readExpressionFile(std::string expr_file);
    void initializeExpression();
    //Each of this function calculates increment (derivative) of a type of gene.
    void intracel_getIncrement(GXMatrix<double>& current_expr, int gene, int k);
    void diffusible_getIncrement(GXMatrix<double>& current_expr, int gene, int k); 
    void cell_property_getIncrement(GXMatrix<double>& current_expr, int gene, int k);
    void edge_property_getIncrement(GXMatrix<double>& current_expr, int gene, int k);
    void vertex_property_getIncrement(GXMatrix<double>& current_expr, int gene, int k);

    double getDegreesFromGradient(int cell);

    void changeVerticesInGrid();
    void changeEdgesInGrid();

  public:
    basicGRN();
    basicGRN(std::string sname);
    basicGRN(std::string sname, Tissue& t, std::string expr_file="");
    void simulateVertexAndNetwork(std::default_random_engine& generator, std::uniform_real_distribution<double>& unif);
    void simulate(double start, double end);
    void runge_kutta_4th();
    void activateAll(GXMatrix<double>& current_expr, int k);
    void checkNegativeValues(GXMatrix<double>& m);
    void setNewParametersToCells();
    GXMatrix<double> getExpression();
    std::string exprToString();
    std::string exprToStringConc(); //Diffusing molecules are divided by cell area
    std::string toString(bool print_current_expression);
    void addCells();  
    void produceOutputs(std::string add_to_name="moved");  


    struct { //This affects which gene index affects each cell/edge property. Should be updated in other implementations
      unsigned int cell_preferred_area = 0;
      unsigned int cell_perimeter_contractility = 1;
      unsigned int edge_tension = 2;
      unsigned int division_angle_random_noise = 3;
      unsigned int division_angle_longest = 4;
      unsigned int division_angle_external = 5;
    } property_index;

   //using increment_function = void(basicGRN::*)(GXMatrix<double>&, int, int, int);
   //std::map<GeneType, increment_function> increment_functions;
    
};


//typedef void (basicGRN::*increment_function)(GXMatrix<double>&, int, int, int);
//typedef std::map<GeneType, increment_function> function_map;
//function_map increment_functions;
/*

    increment_functions.insert(pair<GeneType, increment_function>(GeneType::intracel, basicGRN::intracel_getIncrement));
    increment_functions.insert(pair<GeneType, increment_function>(GeneType::diffusible, basicGRN::diffusible_getIncrement));
    increment_functions.insert(pair<GeneType, increment_function>(GeneType::cell_property, basicGRN::cell_property_getIncrement));
    increment_functions.insert(pair<GeneType, increment_function>(GeneType::edge_property, basicGRN::edge_property_getIncrement));
    increment_functions.insert(pair<GeneType, increment_function>(GeneType::vertex_property, basicGRN::vertex_property_getIncrement));
*/
#endif
