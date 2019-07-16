
#define NON_MATLAB_PARSING
#define MAX_EXT_API_CONNECTIONS 255
#define DO_NOT_USE_SHARED_MEMORY


// uncomment to recognize the sleep function on windows
// #include "windows.h"
#include <vector>
#include <iostream>
#include <yaml-cpp/yaml.h>
#include <sstream>
#include <fstream>
#include "kukaDynRCM.h"
#include "forces.h"
#include "common.h"
#include <extApi.h>
#include <extApiPlatform.h>
#include <extApi.c>
#include <extApiPlatform.c>
#include <boost/numeric/odeint.hpp>
// #include "matlab.hpp"

#define PI 3.1415926535897932384626433832795

extern "C" {
  #include "extApi.h"
}

using namespace std;
using namespace boost::numeric::odeint;

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
const dvec q_ = params["q"].as<dvec>();

const dvec fv_vec_ = params["fv_vec"].as<dvec>();
// Masses vector
const dvec masses_ = params["masses"].as<dvec>();
// Robot links lengths in meters
const dvec lengths_ = params["lengths"].as<dvec>();
// Robot Inertia matrix
const dvec inertiaVecx_ = params["inertiaVecx"].as<dvec>();

static const std::string time_file = "../time.txt";

dvec read_time(){

  string str;
  dvec timeVec;
  std::fstream in_file;
  in_file.open(time_file, std::fstream::in);

  if ( !in_file.is_open() ){
     std::cout << "There was a problem opening the measurements file." << std::endl;
     // return -1;
     exit(0);
   }
  // Reading line into a vector
  else{
  double num = 0.0;
  while (getline(in_file, str)){
    num = std::atof(str.c_str());
    timeVec.push_back(num);
    cout << num << endl;
  }
}
  return timeVec;

}

int main(int argc, char** argv){
  (void)argc;
  (void)argv;

  int clientID = 0;
  int lbrJoint1=0, lbrJoint2=0, lbrJoint3=0, lbrJoint4=0, lbrJoint5=0, lbrJoint6=0, lbrJoint7=0;
  int counter = 0;
  bool WORK = true;
  simxFinish(-1);                                                     //! Close any previously unfinished business
  clientID = simxStart((simxChar*)"127.0.0.1", 19997, true, true, 5000, 5);  //!< Main connection to V-REP
  sleep(2);

  kinParams.l0=lengths_[0];
  kinParams.l1=lengths_[1];
  kinParams.l2=lengths_[2];
  kinParams.l3=lengths_[3];
  kinParams.l4=lengths_[4];
  kinParams.l5=lengths_[5];
  // kinParams.l6=lengths_[6];
  kinParams.l7=lengths_[6];
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

  vecxd qrcm(8);

  double delta_t = 0.0;
  vecxd q_new(8);
  vecxd q_new_(8);
  q_new = kinParams.q;

  dvec time = read_time();
  std::vector<vecxd> positions;
  for (int i=1; i<time.size(); i++){
    cout << "i = " << i << endl;
    // cout << "kin " << kinParams.q[3] << endl;
    qrcm = DYN::RCM::kukaDynRCM(q_new, kinParams);
    cout << "q_rcm " << qrcm << "\n" << endl;
    delta_t = time[i]-time[i-1];


    cout << "delta_t " << delta_t << "\n" << endl;
    q_new += qrcm*delta_t;
    // kinParams.q = q_new;
    for (int i=0; i<q_new.size(); i++){
      q_new[i] = atan2(sin(q_new[i]),cos(q_new[i]));

      cout << "q_new_3 " << q_new_[0]*180/PI << endl;
      q_new_[i] = q_new[i];
    }
    positions.push_back(q_new_);
    cout << "q_new " << q_new << "\n" << endl;
  }
  cout << "final kinParams " << endl;
  for (int i=0; i<kinParams.q.size()-1;i++){

     cout << q_new_[i] << "    " << q_new_[i]*180/PI << endl;
  }
  cout << q_new[kinParams.q.size()-1] << endl;

  if (clientID != -1)
  {
      cout << " Connection status to VREP: SUCCESS" << endl;
      // simxInt syncho = simxSynchronous(clientID, 1);
      // int start = simxStartSimulation(clientID, simx_opmode_oneshot_wait);

      // trigger next simulation step

     simxGetObjectHandle(clientID, "LBR4p_joint1", &lbrJoint1, simx_opmode_oneshot_wait);
     simxGetObjectHandle(clientID, "LBR4p_joint2", &lbrJoint2, simx_opmode_oneshot_wait);
     simxGetObjectHandle(clientID, "LBR4p_joint3", &lbrJoint3, simx_opmode_oneshot_wait);
     simxGetObjectHandle(clientID, "LBR4p_joint4", &lbrJoint4, simx_opmode_oneshot_wait);
     simxGetObjectHandle(clientID, "LBR4p_joint5", &lbrJoint5, simx_opmode_oneshot_wait);
     simxGetObjectHandle(clientID, "LBR4p_joint6", &lbrJoint6, simx_opmode_oneshot_wait);
     simxGetObjectHandle(clientID, "LBR4p_joint7", &lbrJoint7, simx_opmode_oneshot_wait);

     float lbrJoint1Angle, lbrJoint2Angle, lbrJoint3Angle, lbrJoint4Angle, lbrJoint5Angle, lbrJoint6Angle, lbrJoint7Angle;
     simxGetJointPosition(clientID,lbrJoint1,&lbrJoint1Angle,simx_opmode_streaming);
     simxGetJointPosition(clientID,lbrJoint2,&lbrJoint2Angle,simx_opmode_streaming);
     simxGetJointPosition(clientID,lbrJoint3,&lbrJoint3Angle,simx_opmode_streaming);
     simxGetJointPosition(clientID,lbrJoint4,&lbrJoint4Angle,simx_opmode_streaming);
     simxGetJointPosition(clientID,lbrJoint5,&lbrJoint5Angle,simx_opmode_streaming);
     simxGetJointPosition(clientID,lbrJoint6,&lbrJoint6Angle,simx_opmode_streaming);
     simxGetJointPosition(clientID,lbrJoint7,&lbrJoint7Angle,simx_opmode_streaming);

     // update the joint angles:
		bool dataAvailable=(simx_error_noerror==simxGetJointPosition(clientID,lbrJoint1,&lbrJoint1Angle,simx_opmode_buffer));
		dataAvailable=dataAvailable&&(simxGetJointPosition(clientID,lbrJoint2,&lbrJoint2Angle,simx_opmode_buffer));
		dataAvailable=dataAvailable&&(simxGetJointPosition(clientID,lbrJoint3,&lbrJoint3Angle,simx_opmode_buffer));
    dataAvailable=dataAvailable&&(simxGetJointPosition(clientID,lbrJoint4,&lbrJoint4Angle,simx_opmode_buffer));
    dataAvailable=dataAvailable&&(simxGetJointPosition(clientID,lbrJoint5,&lbrJoint5Angle,simx_opmode_buffer));
    dataAvailable=dataAvailable&&(simxGetJointPosition(clientID,lbrJoint6,&lbrJoint6Angle,simx_opmode_buffer));
    dataAvailable=dataAvailable&&(simxGetJointPosition(clientID,lbrJoint7,&lbrJoint7Angle,simx_opmode_buffer));

     for (int i=0; i<positions.size(); i++){
       // simxSynchronous(clientID, 1);
       // simxStartSimulation(clientID, simx_opmode_oneshot_wait);
       simxSetJointTargetPosition(clientID, lbrJoint1, positions[i][0], simx_opmode_oneshot_wait);
       // simxSynchronousTrigger(clientID);
       simxSetJointTargetPosition(clientID, lbrJoint2, positions[i][1], simx_opmode_oneshot_wait);
       // simxSynchronousTrigger(clientID);
       simxSetJointTargetPosition(clientID, lbrJoint3, positions[i][2], simx_opmode_oneshot_wait);
       // simxSynchronousTrigger(clientID);
       simxSetJointTargetPosition(clientID, lbrJoint4, positions[i][3], simx_opmode_oneshot_wait);
       // simxSynchronousTrigger(clientID);
       simxSetJointTargetPosition(clientID, lbrJoint5, positions[i][4], simx_opmode_oneshot_wait);
       // simxSynchronousTrigger(clientID);
       simxSetJointTargetPosition(clientID, lbrJoint6, positions[i][5], simx_opmode_oneshot_wait);
       // simxSynchronousTrigger(clientID);
       simxSetJointTargetPosition(clientID, lbrJoint7, positions[i][6], simx_opmode_oneshot_wait);
       // simxSynchronousTrigger(clientID);
       // usleep(10);
       // simxSpinOnce();
       // simxSleep(50);

     }
     sleep(2);
     simxPauseCommunication(clientID,0);
 }

 else
 {
     cout << " Connection status to VREP: FAILED" << endl;
 }
 simxFinish(clientID);

      // External forces part
      dynParams.N=N_;
      dynParams.fv=fv_;
      dynParams.fc=fc_;
      dynParams.fv_vec = vecxd::Map(fv_vec_.data(), fv_vec_.size());
      dynParams.lengths = vecxd::Map(lengths_.data(), lengths_.size());
      dynParams.masses = vecxd::Map(masses_.data(), masses_.size());

      int element = 0;
      for (int i=0; i<7; i++){
        dvec iVec;
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

      DYN::FORCE::external_forces_estimation(kinParams, dynParams, residual, iMat);
      DYN::RESIDUAL::NewEul_Aux(kinParams, dynParams, residual, iMat);


  return clientID;
}
