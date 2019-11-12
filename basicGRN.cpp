
#include "basicGRN.h"
#include "iostream"

using namespace std;

basicGRN::basicGRN(){
    increment_functions.emplace(GeneType::intracel, &basicGRN::intracel_getIncrement);
    increment_functions.emplace(GeneType::diffusible, &basicGRN::diffusible_getIncrement);
    increment_functions.emplace(GeneType::cell_property, &basicGRN::cell_property_getIncrement);
    increment_functions.emplace(GeneType::edge_property, &basicGRN::edge_property_getIncrement);
    increment_functions.emplace(GeneType::vertex_property, &basicGRN::vertex_property_getIncrement);
}

basicGRN::basicGRN(std::string sname): basicGRN(){
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

    while(it->at(0) != '>') it++;
    it++;
    grn_time_per_step = stod(*it, &sz);

    while(it->at(0) != '>') it++;
    it++;
    vertex_moves_per_step = stoi(*it, &sz);

    while(it->at(0) != '>') it++;
    it++;
    final_moves_accepted = stoi(*it, &sz);

    while(it->at(0) != '>') it++;
    it++;
    write_every_N_moves = stoi(*it, &sz);

    readGeneTypes(it);
    readParams(it);
    readGRN(it);

    time = 0.0;

} //Constructor 1

basicGRN::basicGRN(std::string sname, Tissue& t, std::string expr_file) : basicGRN(sname){
    cell_grid = &t;
    cell_grid->setStepMode(true, write_every_N_moves);
    num_cells = cell_grid->num_cells;
    if(expr_file == ""){
        initializeExpression();
    }
    else{
        readExpressionFile(expr_file);
    }
    for(int k = 0; k < 5; k++)rungekutta_parts.push_back(GXMatrix<double>(num_cells, num_genes, 0.0));
    get_current_cell_grid_params();
} //Constructor 2

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
} //read params


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
} //read single parameter


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
    while(s.at(0)    != '<'){
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
    expression = GXMatrix<double>(num_cells, num_genes, 0.0);
    for(int c = 0; c < num_cells; c++){
        for(int g = 0; g < num_genes; g++){
                expression(c, g) = params.initial_expr[cell_grid->cells[c].type][g];
        }
    }
}

void basicGRN::readExpressionFile(std::string expr_file){
    //To implement
}

void basicGRN::get_current_cell_grid_params(){
    Edge e_aux;
    for(int c = 0; c < num_cells; c++){
        expression(c, property_index.cell_preferred_area) = cell_grid->cells[c].preferred_area;
        expression(c, property_index.cell_perimeter_contractility) = cell_grid->cells[c].perimeter_contractility;

        //Now guess edge tension of cells of this type
        for(int ei = 0; ei < cell_grid->cells[c].num_vertices; ei++){
            e_aux = cell_grid->edges[ cell_grid->cells[c].edges[ei] ];
            if( ( cell_grid->cells[c].type == CellType::blade && e_aux.type == EdgeType::blade ) || 
                    ( cell_grid->cells[c].type == CellType::hinge && e_aux.type == EdgeType::hinge ) ||
                    ( (cell_grid->cells[c].type == CellType::vein_blade || cell_grid->cells[c].type == CellType::vein_hinge) && e_aux.type == EdgeType::vein ) ){
                        expression(c, property_index.edge_tension) = e_aux.tension;
                        break;
            }
        } //for edges in cell

    } //for cells
} // end of get_current_cell_grid_params

void basicGRN::activateAll(GXMatrix<double>& current_expr, int k){
        increment_function f;
        for(int c = 0; c < num_cells; c++){
            for(int g = 0; g < num_genes; g++){
                if(params.constant_expr[cell_grid->cells[c].type][g] < 0){
                   // cout << "    A2, cell of type: " << static_cast<int>(cell_grid->cells[c].type) << ", gene type: " << static_cast<int>(gene_types[g]) << endl;
                    f = increment_functions[gene_types[g]];
                    (this->*f)(current_expr, c, g, k);
                }/*else{
                    rungekutta_parts[k](c, g) = 0;
                }*/ // Not necessary (runge-kutta matrices already initialized to 0) 
        }
    }
}

void basicGRN::simulateVertexAndNetwork(std::default_random_engine& generator, std::uniform_real_distribution<double>& unif){

    cell_grid->setAcceptedMovements(vertex_moves_per_step);
    cell_grid->produceOutputs();
    produceOutputs();

    while(cell_grid->counter_moves_accepted < final_moves_accepted){
        cell_grid->simulate(generator, unif);
        addCells();

        simulate(time, time + grn_time_per_step);
        setNewParametersToCells();
        if(cell_grid->counter_moves_accepted % write_every_N_moves == 0 && cell_grid->counter_moves_accepted > 0){
            cell_grid->produceOutputs();
            produceOutputs();
            //cout << expressionToString() << endl;
        }

        cell_grid->addAcceptedMovements(vertex_moves_per_step);

    } 
}


void basicGRN::simulate(double start, double end){
    for(double t = start; t < end; t += h){
        runge_kutta_4th();
        time += h;
    }
}


void basicGRN::runge_kutta_4th(){
    GXMatrix<double> aux;
    activateAll(expression, 0);
    aux = expression + rungekutta_parts[0]*(0.5*h);
    activateAll(aux, 1);
    aux = expression + rungekutta_parts[1]*(0.5*h);
    activateAll(aux, 2);
    aux = expression + rungekutta_parts[2]*h;
    activateAll(aux, 3);
    rungekutta_parts[4] = (rungekutta_parts[0] + rungekutta_parts[1]*2 + rungekutta_parts[2]*2 + rungekutta_parts[3])/6;

    expression += rungekutta_parts[4]*h;
    // for(int c = 0; c < num_cells; c++) for(int g = 0; g < num_genes; g++) if(expression(c, g) < 0) expression(c, g) = 0; //CAREFUL! don't do this for parameters
    // And maybe change here edge values. Although if it is a waste of time to iterate over all edges when several GRN steps precede a single vertex step, 
    // it is the only way to take into account edges that are between two different cell types or at border. 
    // Unless (as should be done), a different "gene" is used for each type of edge
}

void basicGRN::intracel_getIncrement(GXMatrix<double>& current_expr, int cell, int gene, int k){
    double res = 0;
    CellType ct = cell_grid->cells[cell].type;
    for(int reg : regulators[ct][gene]){
        res += current_expr(cell, reg)*interactions[ct](gene, reg);
    }
    //res = 1/(1 + exp(-1*res)) - params.degr[ct][gene]*current_expr(cell, gene);
    res = tanh(res);
    rungekutta_parts[k](cell, gene) = res < 0 ? 0 - params.degr[ct][gene]*current_expr(cell, gene) : res - params.degr[ct][gene]*current_expr(cell, gene);
}

void basicGRN::diffusible_getIncrement(GXMatrix<double>& current_expr, int cell, int gene, int k){
    intracel_getIncrement(current_expr, cell, gene, k);
    //cout << "        A2.6" << endl;
    const Cell *c = &cell_grid->cells[cell];
    //cout << "        A2.7" << endl;
    const Edge *edge;
    double diff = 0;
    int nei;
    //cout << "        A2.8" << endl;
    for(int e = 0; e < c->num_vertices; e++){
        edge = &cell_grid->edges[c->edges[e]];
        nei = edge->cells[0] == c->ind ? edge->cells[1] : edge->cells[0];
        if(nei == EMPTY_CONNECTION){
            continue;
        }else{
        //cout << "            A2.9 cell: " << cell << ", gene: " << gene << ", nei: " << nei <<endl;
        //diff += (current_expr(cell, gene)/c->area - current_expr(nei, gene)/cell_grid->cells[nei].area)*edge->length/c->perimeter; //normalizing by area?
            diff += (current_expr(cell, gene)- current_expr(nei, gene))*edge->length/c->perimeter;
        }
        //cout << "            A2.10" << endl;
    }
    //cout << "        A2.11" << endl;
    rungekutta_parts[k](cell, gene) = rungekutta_parts[k](cell, gene) - params.diff_rate[c->type][gene]*diff; //*power(diff, 2);
    //cout << "        A2.12" << endl;
} //diffussible

void basicGRN::cell_property_getIncrement(GXMatrix<double>& current_expr, int cell, int gene, int k){
    //cout << "        A2.13" << endl;
    intracel_getIncrement(current_expr, cell, gene, k);
}

void basicGRN::edge_property_getIncrement(GXMatrix<double>& current_expr, int cell, int gene, int k){
    //cout << "        A2.14" << endl;
    intracel_getIncrement(current_expr, cell, gene, k);
}
void basicGRN::vertex_property_getIncrement(GXMatrix<double>& current_expr, int cell, int gene, int k){
    //cout << "        A2.15" << endl;
    intracel_getIncrement(current_expr, cell, gene, k);
}
void basicGRN::setNewParametersToCells(){
    for(int i = 0; i < num_cells; i++){
        if(!cell_grid->cells[i].dead){
            cell_grid->cells[i].preferred_area = expression(i, property_index.cell_preferred_area);
            cell_grid->cells[i].perimeter_contractility = expression(i, property_index.cell_perimeter_contractility);
        }
    }

    //To implement: changes in vertices or edges, orientation with respect to a morphogen gradient
}

void basicGRN::addCells(){
    
    if(cell_grid->past_divisions.empty()){
        return;
    }
    for(int i = 0; i < 5; i++) rungekutta_parts[i].add_row(cell_grid->past_divisions.size(), 0.0);
    DivisionRecord d;
    vector<double> row;
    
    while(!cell_grid->past_divisions.empty()){
        d = cell_grid->past_divisions.front();
        row = expression.getRow(d.parent);
        expression.add_row(row);
        cell_grid->past_divisions.pop();
        num_cells++;
    }
}    

GXMatrix<double> basicGRN::getExpression(){
    return expression;
}

std::string basicGRN::exprToString(){
    return expression.toString();
}

void basicGRN::produceOutputs(std::string add_to_name){

    std:string fname = name + "_" + add_to_name + "_" + std::to_string(int(cell_grid->counter_moves_accepted / write_every_N_moves));
    std::string s = toString(true);
    ofstream fo;
    fo.open(fname + OUTPUT_EXTENSION);
    fo << s;
    fo.close();

    fo.open(fname + EXPRESSION_EXTENSION);
    s = expressionToString();
    fo << s;
    fo.close();
}

std::string basicGRN::expressionToString(){
    return expression.toString();
}
std::string basicGRN::toString(bool print_current_expression=false){
        std::string s = name + " time: " + to_string(time) + ", number of genes (and variable properties): " + to_string(num_genes) + ".\n\n";
        s += "Gene type (intracel = 0, diffusible = 1, cell_property = 2, edge_property = 3, vertex_property = 4):\n";
        for(GeneType g : gene_types) s += to_string( static_cast<int>(g) ) + ", ";

        s += "\n\nPARAMETERS (cell types: blade = 0, hinge = 1, vein = 2, vein_hinge = 3):\n    degradation rates:\n";
        for (auto const& x : params.degr){
            s += "cell type " + to_string( static_cast<int>(x.first) ) + ": ";
            for(double y : x.second) s+= to_string(y) + ", ";
            s += "\n";
        }

        s += "\n    diffusion rates: \n";
        for (auto const& x : params.diff_rate){
            s += "cell type " + to_string( static_cast<int>(x.first) ) + ": ";
            for(double y : x.second) s+= to_string(y) + ", ";
            s += "\n";
        }

        s += "\n    Initial expression: \n";
        for (auto const& x : params.initial_expr){
            s += "cell type " + to_string( static_cast<int>(x.first) ) + ": ";
            for(double y : x.second) s+= to_string(y) + ", ";
            s += "\n";
        }

        s += "\n    Constant expression: \n";
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
        if(print_current_expression) s += expression.toString() + "\n*****\n";
        return s;
}


