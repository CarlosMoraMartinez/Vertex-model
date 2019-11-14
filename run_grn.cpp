
#include "basicGRN.h"

using namespace std;


int main(int argc, char *argv[]){

  std::string inputfile_vertex = argv[1];
  std::string inputfile_grn = argv[2];
  cout << "a" << endl;
  Tissue t(inputfile_vertex);
  cout << "b" << endl;
  basicGRN g(inputfile_grn, t); //third argument can be an expression file
  cout << "c" << endl;

  srand(RANDOM_SEED);
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif;

  cout << "d" << endl;
  g.simulateVertexAndNetwork(generator, unif);
 /* cout << g.toString(true);
  g.simulate(0.0, 100);
  cout << "\n0 to 0.1:\n" << g.exprToString() << endl;

  g.simulate(100, 1000);
  cout << "\n0.1 to 0.2:\n" << g.exprToString() << endl;
*/
}
