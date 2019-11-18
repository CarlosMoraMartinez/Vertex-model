
#include "basicGRN.h"

using namespace std;


int main(int argc, char *argv[]){

  std::string inputfile_vertex = argv[1];
  std::string inputfile_grn = argv[2];
  std::string inputfile_vert_parms = argv[3];

  Tissue t(inputfile_vertex, inputfile_vert_parms);

  basicGRN g(inputfile_grn, t); //third argument can be an expression file


  srand(RANDOM_SEED);
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif;


  g.simulateVertexAndNetwork(generator, unif);
 /* cout << g.toString(true);
  g.simulate(0.0, 100);
  cout << "\n0 to 0.1:\n" << g.exprToString() << endl;

  g.simulate(100, 1000);
  cout << "\n0.1 to 0.2:\n" << g.exprToString() << endl;
*/
}
