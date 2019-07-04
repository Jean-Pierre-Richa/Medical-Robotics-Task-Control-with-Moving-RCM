// main.cpp

#include <vector>
#include <iostream>
#include <yaml-cpp/yaml.h>
#include "kukaDynRCM.h"
#include "forces.h"

using namespace std;

DYN::Params kinParams;
DYN::dynamicParams dynParams;
DYN::inertiaMat iMat;

// Config file dir
static const std::string CONFIG_FILE = "../config/config.yaml";
// Load the config file
YAML::Node params = YAML::LoadFile(CONFIG_FILE);

const int N_ = params["N"].as<int>();

const double fv_ = params["fv"].as<double>();

const double fc_ = params["fc"].as<double>();

const double off_ = params["off"].as<double>();
// Target position on the X axis
const double Xtr_ = params["Xtr"].as<double>();
// Target position on the Y axis
const double Ytr_ = params["Ytr"].as<double>();
// Target position on the Z axis
const double Ztr_ = params["Ztr"].as<double>();
// Control gain
const double gain_ = params["gain"].as<double>();
// Robot's joints state vector
std::vector<double> q_ = params["q"].as<std::vector<double>>();

std::vector<double> fv_vec_ = params["fv_vec"].as<std::vector<double>>();
// Masses vector
std::vector<double> masses_ = params["masses"].as<std::vector<double>>();
// Robot links lengths in meters
std::vector<double> lengths_ = params["lengths"].as<std::vector<double>>();
// Robot Inertia matrix
std::vector<double> inertiaVecx_ = params["inertiaVecx"].as<std::vector<double>>();

int main(int argc, char** argv){
  (void)argc;
  (void)argv;

  // vecxd q_(8);
  // q_ << qc[0],qc[1],qc[2],qc[3],qc[4],qc[5],qc[6],qc[7];

  kinParams.l0=lengths_[0];
  kinParams.l1=lengths_[1];
  kinParams.l2=lengths_[2];
  kinParams.l3=lengths_[3];
  kinParams.l4=lengths_[4];
  kinParams.l5=lengths_[5];
  kinParams.l6=lengths_[6];
  kinParams.l7=lengths_[7];
  kinParams.off=off_;
  kinParams.Xtr=Xtr_;
  kinParams.Ytr=Ytr_;
  kinParams.Ztr=Ztr_;
  kinParams.q = vecxd::Map(q_.data(), q_.size());
  kinParams.target_pos << kinParams.Xtr, kinParams.Ytr, kinParams.Ztr;
  kinParams.gain=gain_;
  MatXd K_rcm = kinParams.gain*MatXd::Identity(3,3);
  MatXd K_t = kinParams.gain*MatXd::Identity(3,3);
  kinParams.k.block<3,3>(0,0)=K_rcm;
  kinParams.k.block<3,3>(3,3)=K_t;

  DYN::RCM::kukaDynRCM(kinParams);

  // External forces part
  dynParams.N=N_;
  dynParams.fv=fv_;
  dynParams.fc=fc_;
  dynParams.fv_vec = vecxd::Map(fv_vec_.data(), fv_vec_.size());
  dynParams.lengths = vecxd::Map(lengths_.data(), lengths_.size());
  dynParams.masses = vecxd::Map(masses_.data(), masses_.size());

  int element = 0;
  for (int i=0; i<7; i++){
    std::vector<double> iVec;
    for (int j=0; j<6; j++){
      iVec.push_back(inertiaVecx_[element]);
      element++;
    }
    dynParams.inertiaVectors.block<1,6>(i,0) << iVec[0], iVec[1], iVec[2], iVec[3], iVec[4], iVec[5];
  }
  // cout << dynParams.inertiaVectors << endl;

  vecEigens V(7);

  for (int i=0; i<7; i++){
    V[i] = dynParams.inertiaVectors.block<1,6>(i,0);
    // cout << "V" << "[" << i << "]" << V[i] << endl;
  }

  iMat.I1 = DYN::FORCE::calcInertiaMat(V[0]);
  iMat.I2 = DYN::FORCE::calcInertiaMat(V[1]);
  iMat.I3 = DYN::FORCE::calcInertiaMat(V[2]);
  iMat.I4 = DYN::FORCE::calcInertiaMat(V[3]);
  iMat.I5 = DYN::FORCE::calcInertiaMat(V[4]);
  iMat.I6 = DYN::FORCE::calcInertiaMat(V[5]);
  iMat.I7 = DYN::FORCE::calcInertiaMat(V[6]);


  DYN::FORCE::external_forces_estimation(dynParams, iMat);

}
