// main.cpp


#define NON_MATLAB_PARSING
#define MAX_EXT_API_CONNECTIONS 255
#define DO_NOT_USE_SHARED_MEMORY

// uncomment to recognize the sleep function on windows
// #include "windows.h"
#include <vector>
#include <iostream>
#include <yaml-cpp/yaml.h>
#include "kukaDynRCM.h"
#include "forces.h"
#include "common.h"
#include <extApi.h>
#include <extApiPlatform.h>
#include <extApi.c>
#include <extApiPlatform.c>


extern "C" {
  #include "extApi.h"
}
#define PI 3.14
using namespace std;

DYN::Params kinParams;
DYN::dynamicParams dynParams;
DYN::inertiaMat iMat;
DYN::residualParams residual;

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

  // int clientID = 0;
  int lbrJoint1=0, lbrJoint2=0, lbrJoint3=0, lbrJoint4=0, lbrJoint5=0, lbrJoint6=0, lbrJoint7=0;
  int counter = 0;
  bool WORK = true;
  simxFinish(-1);                                                     //! Close any previously unfinished business
  sleep(2);

  int clientID=simxStart((simxChar*)"127.0.0.1",19997,true,true,5000,5);
  	if (clientID==-1){
  		cout << "Could not connect to V-REP remote API server " << endl;
  		simxFinish(clientID);
  	}else{
  		cout << "Connected to remote API server" << endl;
  	}
  simxStartSimulation(clientID,simx_opmode_oneshot);
  cout << "V-REP simulation start" << endl;

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

  if (clientID != -1)
  {
      cout << " Connection status to VREP: SUCCESS" << endl;

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
      DYN::RESIDUAL::NewEul_Aux(kinParams, dynParams, residual, iMat);


      simxInt syncho = simxSynchronous(clientID, 1);
      int start = simxStartSimulation(clientID, simx_opmode_oneshot_wait);

      simxGetObjectHandle(clientID, "LBR4p_joint1", &lbrJoint1, simx_opmode_oneshot_wait);
      simxGetObjectHandle(clientID, "LBR4p_joint2", &lbrJoint2, simx_opmode_oneshot_wait);
      simxGetObjectHandle(clientID, "LBR4p_joint3", &lbrJoint3, simx_opmode_oneshot_wait);
      simxGetObjectHandle(clientID, "LBR4p_joint4", &lbrJoint4, simx_opmode_oneshot_wait);
      simxGetObjectHandle(clientID, "LBR4p_joint5", &lbrJoint5, simx_opmode_oneshot_wait);
      simxGetObjectHandle(clientID, "LBR4p_joint6", &lbrJoint6, simx_opmode_oneshot_wait);
      simxGetObjectHandle(clientID, "LBR4p_joint7", &lbrJoint7, simx_opmode_oneshot_wait);

      simxSetJointTargetPosition(clientID, lbrJoint1, 90.0* (PI / 180), simx_opmode_oneshot_wait);
      simxSetJointTargetPosition(clientID, lbrJoint2, 90.0* (PI / 180), simx_opmode_oneshot_wait);
      simxSetJointTargetPosition(clientID, lbrJoint3, 170.0* (PI / 180), simx_opmode_oneshot_wait);
      simxSetJointTargetPosition(clientID, lbrJoint4, -90.0* (PI / 180), simx_opmode_oneshot_wait);
      simxSetJointTargetPosition(clientID, lbrJoint5, 90.0* (PI / 180), simx_opmode_oneshot_wait);
      simxSetJointTargetPosition(clientID, lbrJoint6, 90.0* (PI / 180), simx_opmode_oneshot_wait);
      simxSetJointTargetPosition(clientID, lbrJoint7, 0.0* (PI / 180), simx_opmode_oneshot_wait);
      sleep(2);
      simxPauseCommunication(clientID,0);
  }

  else
  {
      cout << " Connection status to VREP: FAILED" << endl;
  }
  simxFinish(clientID);
  return clientID;
}



// }
