
#include "basicGRN.h"

using namespace std;


int main(int argc, char *argv[]){

  cout << "A" << endl;
  Tissue t("gm0_s3.0_3x3_n0.3", 100, 20);
  cout << "B" << endl;
  basicGRN g("testgrn", t, "");
  cout << "C" << endl;
  srand(RANDOM_SEED);
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif;

  cout << g.toString(true);
  g.simulate(0.0, 100);
  cout << "\n0 to 0.1:\n" << g.exprToString() << endl;

  g.simulate(100, 1000);
  cout << "\n0.1 to 0.2:\n" << g.exprToString() << endl;

}
