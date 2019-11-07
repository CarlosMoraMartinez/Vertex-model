
#include "basicGRN.h"
#include "iostream"

using namespace std;

basicGRN::basicGRN(){}

basicGRN::basicGRN(std::string sname){
  name = sname;
  string grnfile = sname + INTERACTIONS_EXTENSION;
  std::ifstream fin_grn;
  fin_grn.open(grnfile);
  std::string line;
  std::vector<std::string> inp;
  std::string::size_type sz;

  cout << "Reading .grn file...\n";
  while(getline(fin_grn, line)) if(!(line.empty() || line.find_first_not_of(' ') == std::string::npos)) if(line.at(0) != '#') inp.push_back(line);

  std::vector<std::string>::iterator it = inp.begin(); 

  //First read value of h
  while(it->at(0) != '>') it++;
  it++;

  h = stod(*it, &sz);

  readGeneTypes(it);
  readParams(it);
  readGRN(it);

}

basicGRN::basicGRN(std::string sname, Tissue& t, std::string expr_file="") : basicGRN(sname){
  cell_grid = &t;
  if(expr_file == ""){
    initializeExpression();
  }
  else{
    readExpressionFile(expr_file);
  }
  get_current_cell_grid_params();
}

void basicGRN::readGeneTypes( std::vector<std::string>::iterator& it){
  std::string s;
  std::string::size_type sz;

  while(it->at(0) != '>') it++;
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
  params.diff_rate = readSingleParam(it);
  params.initial_expr = readSingleParam(it);
  params.constant_expr = readSingleParam(it);
}

// The return type is defined in basicGRN.h as:
//typedef std::map<CellType, std::vector<double> > ctparams;
ctparams basicGRN::readSingleParam( std::vector<std::string>::iterator& it){
  std::string s;
  std::string::size_type sz;
  ctparams cell_parms;
  CellType celltype;
  std::vector<double> values;
  double auxd;

  while(it->at(0) != '>') it++;
  it++;

  while(it->at(0) != '<'){
    s = *it;
    celltype = static_cast<CellType>( stoi(s, &sz) );
    s = s.substr(sz);
    while(s.length() != 0){
      std::string prueba = "a5.5";
      s.erase(0, 1);
      auxd = stod(s, &sz);
      values.push_back(auxd);
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
  CellType celltype;
  std::vector<double> values;
  GXMatrix<double> aux(num_genes, num_genes, 0);

  list_of_regulators regs_in_cell;
  for(int i = 0; i < num_genes; i++) regs_in_cell.push_back(vector<int>());

  while(it->at(0) != '>') it++;
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
          if(aux(i,j) != 0) regs_in_cell[i].push_back(j); //will make access faster during simulation
          s = s.substr(sz);
        }
    }
    interactions.insert(pair<CellType, GXMatrix<double> >(celltype, aux));
    regulators.insert(pair<CellType, list_of_regulators >(celltype, regs_in_cell));
    for(int i = 0; i < num_genes; i++) regs_in_cell[i].clear();
    it++;
    s = *it;
 }//while no '<'

}


void basicGRN::initializeExpression(){

}
void basicGRN::readExpressionFile(std::string expr_file){

}
void basicGRN::get_current_cell_grid_params(){

}

void basicGRN::activateAll(GXMatrix<double> current_expr, int k){

}
void basicGRN::runge_kutta_4th(){
   
}

std::string basicGRN::toString(){
    std::string s = name + " number of genes (+ variable properties): " + to_string(num_genes) + ".\n\n";
    s += "Gene type (intracel = 0, diffusible = 1, cell_property = 2, edge_property = 3, vertex_property = 4):\n";
    for(GeneType g : gene_types) s += to_string( static_cast<int>(g) ) + ", ";

    s += "\n\nPARAMETERS (cell types: blade = 0, hinge = 1, vein = 2, vein_hinge = 3):\n  degradation rates:\n";
    for (auto const& x : params.degr){
      s += "cell type " + to_string( static_cast<int>(x.first) ) + ": ";
      for(double y : x.second) s+= to_string(y) + ", ";
      s += "\n";
    }

    s += "\n  diffusion rates: \n";
    for (auto const& x : params.diff_rate){
      s += "cell type " + to_string( static_cast<int>(x.first) ) + ": ";
      for(double y : x.second) s+= to_string(y) + ", ";
      s += "\n";
    }

    s += "\n  Initial expression: \n";
    for (auto const& x : params.initial_expr){
      s += "cell type " + to_string( static_cast<int>(x.first) ) + ": ";
      for(double y : x.second) s+= to_string(y) + ", ";
      s += "\n";
    }

    s += "\n  Constant expression: \n";
    for (auto const& x : params.constant_expr){
      s += "cell type " + to_string( static_cast<int>(x.first) ) + ": ";
      for(double y : x.second) s+= to_string(y) + ", ";
      s += "\n";
    }

    s += "\n\nMATRIX OF INTERACTIONS BY CELL (cell types: blade = 0, hinge = 1, vein = 2, vein_hinge = 3):\n";
    GXMatrix<double> aux;
    for (auto const& x : interactions){
        aux = x.second;
        s += "cell type " + to_string( static_cast<int>(x.first) ) + ":\n" + aux.toString() + "\n";
    }

    s += "\n\nLIST OF REGULATORS BY GENE/PROPERTY BY CELL (cell types: blade = 0, hinge = 1, vein = 2, vein_hinge = 3):\n";
    for (auto const& x : regulators){
      s += "cell type " + to_string( static_cast<int>(x.first) ) + ":\n";
      for (int g = 0; g < x.second.size(); g++){
        s += "gene " + to_string(g) + ":\t";
        for (int r : x.second[g]){
          s += to_string(r) + ", ";           
        }
        s += "\n";
      }  
    }
    s += "\n*****\n";
    return s;
}


