
#include "basicGRN.h"
#include "VertexSystem.h"
#include "GXMatrix.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <random>

basicGRN::basicGRN(){}

basicGRN::basicGRN(std::string sname){
  string grnfile = sname + INTERACTIONS_EXTENSION;
  std::ifstream fin_grn;
  fin_grn.open(grnfile);
  std::string line;
  std::vector<std::string> inp
  std::string::size_type sz;

  while(getline(fin_grn, line)) if(!(line.empty() || line.find_first_not_of(' ') == std::string::npos)) if(line.at(0) != '#') inp.push_back(line);
  std::vector<std::string>::iterator it = inp.begin(); 
  //First read value of h
  while(*it.at(0) != '>') it++;
  it++;
  h = stod(*it, &sz);

  readGeneTypes(it);
  readParams(it);
  readGRN(it);

}

basicGRN(std::string sname, Tissue& t) : basicGRN(sname){
  cell_grid = &t;
  initializeExpression();
  
}

void basicGRN::readGeneTypes( std::vector<std::string>::iterator& it){
  std::string s;
  std::string::size_type sz;

  while(*it.at(0) != '>') it++;
  it++;
  s = *it;
  
  while(s.length() != 0){
    gene_types.push_back(static_cast<GeneType>(stoi(s, &sz)));
    s = s.substr(sz);
  }
  num_genes = gene_types.size();

}//end readGeneTypes


void basicGRN::readParams( std::vector<std::string>::iterator& it){
  params.degr = readSingleParam(it);
  params.initial_expr = readSingleParam(it);
  params.constant_expr = readSingleParam(it);
}

// The return type is defined in basicGRN.h as:
//typedef std::map<CellType, std::vector<double> > ctparams;
basicGRN::ctparams basicGRN::readSingleParam( std::vector<std::string>::iterator& it){
  std::string s;
  std::string::size_type sz;
  ctparams cell_parms;
  CellType celltype;
  std::vector<double> values;

  while(*it.at(0) != '>') it++;
  it++;

  while(*it.at(0) != '<'){
    s = *it;
    celltype = static_cast<CellType>( stoi(s, &sz) );
    s = s.substr(sz);
    while(s.length() != 0){
      values.push_back(stod(s, &sz));
      s = s.substr(sz);
    }
    cell_parms.insert(pair<CellType, vector<double> >(celltype, values));
    values.clear();
    it++;
 }//while no <
 return cell_parms;
}


//Need to fill a member of type typedef std::map<CellType, GXMatrix<double> > ctinteractions;

void basicGRN::readGRN( std::vector<std::string>::iterator& it){

  std::string s = "";
  std::string::size_type sz;
  ctparams cell_parms;
  CellType celltype;
  std::vector<double> values;
  GXMatrix<double> aux(num_genes, num_genes, 0);

  list_of_regulators regs_in_cell;
  for(int i = 0; i < num_genes; i++) regs_in_cell.push_back(vector<int>());

  while(*it.at(0) != '>') it++;
  it++;
  s = *it;
  while(s.at(0)  != '<'){
    celltype = static_cast<CellType>( stoi(s, &sz) );
    //iterate over interaction matrix. Row means activated gene, col means activator
    for(int i = 0; i < num_genes; i++){
        it++;
        s = *it;
        for(int j = 0; j < num_genes; j++){
          aux(i,j) = stod(s, &sz);
          if(aux(i,j) != 0) regs_in_cell[i].push_back(j); //will make access faster
          s = s.substr(sz);
        }
    }
    interactions.insert(pair<CellType, GXMatrix<double> >(celltype, aux));
    regulators.insert(pair<CellType, list_of_regulators >(celltype, regs_in_cell));
    for(int i = 0; i < num_genes; i++) regs_in_cell[i].clear().
    it++;
    s = *it;
 }//while no <

}

int basicGRN::activateAll(GXMatrix<double> current_expr, int k);
void basicGRN::runge_kutta_4th(){
   
}



