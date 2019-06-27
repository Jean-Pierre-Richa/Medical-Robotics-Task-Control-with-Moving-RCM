// main.cpp

#include <iostream>
#include <vector>
#include "kukaDynRCM.h"

using namespace std;

int main(int argc, char** argv){
  (void)argc;
  (void)argv;

  vecxd q(8);
  q << 1, 2, 3, 4, 5, 6, 7, 8;

  DYN::RCM::kukaDynRCM(q);

}
