#include <TF1.h>
#include <TF2.h>

#include "constants.h"

using namespace std;


double gaus2(double* x, double* par) {
// par[0] = Scale (coefficient)
// par[1] = x center
// par[2] = y center
// par[3] = standard deviation
// x[0] = x / x[1] = y
  return par[0]*exp(-(pow((x[0]-par[1]),2)+pow((x[1]-par[2]),2))/(2*par[3]*par[3]));
}


double gaus1(double* x, double* par) {
// par[0] = Scale (coefficient)
// par[1] = x center
// par[3] = standard deviation
// x[0] = x 
  return par[0]*exp(-pow((x[0]-par[1]),2)/(2*par[3]*par[3]));
}



TF2 *f2Gaus = new TF2("f2Gaus", gaus2, 0, 1, 0, 1, 4);
// must reset range: SetRange (Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax)


TF1 *f1Gaus = new TF1("f1Gaus", gaus1, 0, 1, 3);
// must reset range: SetRange (Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax)
