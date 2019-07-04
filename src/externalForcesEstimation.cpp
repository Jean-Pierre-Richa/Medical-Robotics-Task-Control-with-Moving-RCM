#include <iostream>
#include "forces.h"

using namespace std;

void DYN::FORCE::external_forces_estimation(dynamicParams& dynParams, inertiaMat& iMat){

  // Number of joints
  int N=dynParams.N;
  // Actuators torques
  MatXd tau_act(dynParams.N,1);
  double fv=dynParams.fv, fc=dynParams.fc;

  vecxd fv_vec(7);
  fv_vec = dynParams.fv_vec;
  // Links Lengths
  vecxd lengths(7);
  lengths = dynParams.lengths;
  // Links masses
  vecxd masses(7);
  masses = dynParams.masses;
  // cout << "masses " << masses << endl;

}
