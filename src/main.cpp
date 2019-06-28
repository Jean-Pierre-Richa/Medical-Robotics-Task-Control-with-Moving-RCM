// main.cpp

#include <iostream>
#include <vector>
#include "kukaDynRCM.h"

using namespace std;

DYN::Params kinParams;

int main(int argc, char** argv){
  (void)argc;
  (void)argv;

  kinParams.l0=1;
  kinParams.l1=2;
  kinParams.l2=3;
  kinParams.l3=4;
  kinParams.l4=5;
  kinParams.l5=6;
  kinParams.l6=7;
  kinParams.l7=8;
  kinParams.off=10;
  kinParams.Xtr = 3;
  kinParams.Ytr = 6;
  kinParams.Ztr = 5;
  kinParams.target_pos << 1, 4, 7;
  kinParams.q << 1, 2, 3, 4, 5, 6, 7, 8;
  kinParams.k << 0.1,0,0,0,0,0,
                 0,0.1,0,0,0,0,
                 0,0,0.1,0,0,0,
                 0,0,0,0.1,0,0,
                 0,0,0,0,0.1,0,
                 0,0,0,0,0,0.1;

  DYN::RCM::kukaDynRCM(kinParams);


}
