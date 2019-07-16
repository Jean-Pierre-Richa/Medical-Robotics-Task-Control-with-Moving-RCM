
/**
 * \file      kukaDynRCM.h
 *
 * \brief     C++ Header File template
 *
 * \author    Jean Pierre Richa
 *
 * \Project:  Medical Robotics
 */

#ifndef _KUKA_DYN_RCM_
#define _KUKA_DYN_RCM_

/*******************************************************************************
 *       I N C L U D E - F I L E S
 ******************************************************************************/

#include <Eigen/Dense>
#include "common.h"


/*******************************************************************************
 *       N A M E S P A C E
 ******************************************************************************/
// static void rhs(dvec& x, dvec& dxdt, const double t );
// static void write_cout(dvec& x, const double t );

namespace DYN{
/*
 * \Struct holding the kinematic parameters
 */
struct Params{
  double gain, l0, l1, l2, l3, l4, l5, l6, l7, off, Xtr, Ytr, Ztr;
  int ee_length;
  vecxd target_pos;
  vecxd q;
  MatXd k;
  Params() : gain(0), l0(0), l1(0), l2(0), l3(0), l4(0), l5(0), l6(0), l7(0), off(0), Xtr(0), Ytr(0), Ztr(0), k(6,6), q(8), target_pos(3){}
};

/**
 * \class kukaDyn
 *
 * \Calculates the direct kinematics and jacobian matrices of the kuka LWR for
 * \the 6th and 7th joints
 *
 */

class RCM{
  public:
    static vecxd kukaDynRCM(const vecxd& q, Params& kinParams);
    static vec3Eigens dirKin(const Params& kinParams);


  }; // class kukaDyn

} //namespace
#endif // _KUKA_DYN_RCM_
