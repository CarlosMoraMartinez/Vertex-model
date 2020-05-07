

#include "GeneExpression.h"


GeneExpression::GeneExpression(){
	expr = new vector<double>();
	new_expr = new vector<double>();
}


GeneExpression::GeneExpression(int n){
	expr = new vector<double>(n, 0);
	new_expr = new vector<double>(n, 0);


}

GeneExpression::GeneExpression(int n, const double* e){
	expr = new vector<double>(n, 0);
	new_expr = new vector<double>(n, 0);
	for(int i = 0; i < n; i++) (*expr)[i] = e[i];
}

GeneExpression::GeneExpression(const GeneExpression &e){
	new_expr = new vector<double>(e.size(), 0);
	expr = new vector<double>(e.size(), 0); 
	for(int i = 0; i < e.size(); i++) (*expr)[i] = e.getExpr(i); 
}

GeneExpression::~GeneExpression(){
	(*expr).~vector<double>();
	(*new_expr).~vector<double>();
}

void GeneExpression::setNewExpr(int n, double val){
	(*new_expr)[n] = val;
}

void GeneExpression::setNewExpr(const double* val){
	for(int n = 0; n < this->size(); n++) (*new_expr)[n] = val[n];
}

void GeneExpression::update(){
	swap(expr, new_expr);
}


double& GeneExpression::operator[] (int index){
 return (*expr)[index];

}

GeneExpression& GeneExpression::operator=(const GeneExpression &gx){
	new_expr = new expression_v((*expr).size(), 0);
	expr = new expression_v(gx.size(), 0);
	for(int i = 0; i < gx.size(); i++) (*expr)[i] = gx.getExpr(i); 
}

int GeneExpression::size() const{
	return (*expr).size();
}

expression_v GeneExpression::getExpr(){
	return *expr;
}

expression_v GeneExpression::getNewExpr(){
	return *new_expr;
}

double GeneExpression::getExpr(int n) const{
	return (*expr)[n];
}

double GeneExpression::getNewExpr(int n) const{ 
	return (*new_expr)[n];
}

std::string GeneExpression::toString(){
	std::string s = "expr: ";
	for(double n : *expr) s += std::to_string(n) + ", ";
	s += "\tNewExpr: ";
	for(double n : *new_expr) s += std::to_string(n) + ", ";
	return s;
}


