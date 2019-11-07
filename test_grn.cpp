
#include "basicGRN.h"

using namespace std;


int main(int argc, char *argv[]){

  Tissue t("gm0_s3.0_20x30_n0.4", 100, 20);
  basicGRN g("testgrn", t, "");
  cout << g.toString(true);

}
