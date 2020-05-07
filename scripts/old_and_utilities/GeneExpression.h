

#include <vector>
#include <string>
#include <iostream>

using namespace std;

typedef std::vector<double> expression_v;
typedef expression_v* expression_vptr;

class GeneExpression{

	public:
		GeneExpression();
		GeneExpression(int n);
		GeneExpression(int n, const double* e);
		GeneExpression(const GeneExpression &e);
		~GeneExpression();
		void update();
		double& operator[] (int index);
		GeneExpression& operator=(const GeneExpression &gx);
		int size() const;
		void setNewExpr(int n, double val);
		void setNewExpr(const double* val);
		expression_v getExpr();
		expression_v getNewExpr();
		double getExpr(int n) const;
		double getNewExpr(int n) const;
		std::string toString();
	private:		
		expression_vptr expr;
		expression_vptr new_expr;
};
