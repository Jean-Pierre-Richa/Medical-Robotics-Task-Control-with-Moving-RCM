// main.cpp

#include <vector>
#include <iostream>
#include <yaml-cpp/yaml.h>
#include "kukaDynRCM.h"

using namespace std;

DYN::Params kinParams;

// Config file dir
static const std::string CONFIG_FILE = "../config/config.yaml";
// Load the config file
YAML::Node params = YAML::LoadFile(CONFIG_FILE);
// Robot links lengths in meters
const double l0_ = params["link0"].as<double>();
// Length of link 1 in meters
const double l1_ = params["link1"].as<double>();
// Length of link 2 in meters
const double l2_ = params["link2"].as<double>();
// Length of link 3 in meters
const double l3_ = params["link3"].as<double>();
// Length of link 4 in meters
const double l4_ = params["link4"].as<double>();
// Length of link 5 in meters
const double l5_ = params["link5"].as<double>();
// Length of link 6 in meters
const double l6_ = params["link6"].as<double>();
// Length of link 7 in meters
const double l7_ = params["link7"].as<double>();
const double off_ = params["off"].as<double>();
// Target position on the X axis
const double Xtr_ = params["Xtr"].as<double>();
// Target position on the Y axis
const double Ytr_ = params["Ytr"].as<double>();
// Target position on the Z axis
const double Ztr_ = params["Ztr"].as<double>();
// Robot's joints state vector
std::vector<double> qc = params["q"].as<std::vector<double>>();
// Control gain
const double gain_ = params["gain"].as<double>();

int main(int argc, char** argv){
  (void)argc;
  (void)argv;

  vecxd q_(8);
  q_ << qc[0],qc[1],qc[2],qc[3],qc[4],qc[5],qc[6],qc[7];

  kinParams.l0=l0_;
  kinParams.l1=l1_;
  kinParams.l2=l2_;
  kinParams.l3=l3_;
  kinParams.l4=l4_;
  kinParams.l5=l5_;
  kinParams.l6=l6_;
  kinParams.l7=l7_;
  kinParams.off=off_;
  kinParams.Xtr=Xtr_;
  kinParams.Ytr=Ytr_;
  kinParams.Ztr=Ztr_;
  kinParams.target_pos << kinParams.Xtr, kinParams.Ytr, kinParams.Ztr;
  kinParams.q=q_;
  kinParams.gain=gain_;

  MatXd K_rcm = kinParams.gain*MatXd::Identity(3,3);
  MatXd K_t = kinParams.gain*MatXd::Identity(3,3);
  kinParams.k.block<3,3>(0,0)=K_rcm;
  kinParams.k.block<3,3>(3,3)=K_t;

  DYN::RCM::kukaDynRCM(kinParams);

}
