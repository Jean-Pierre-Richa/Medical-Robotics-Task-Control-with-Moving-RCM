
/**
 * \file      featuresMatcher.h
 *
 * \brief     C++ Header File template
 *
 * \author    Jean Pierre Richa
 *
 * \Project:  Medical Robotics
 */

/*******************************************************************************
 *       I N C L U D E - F I L E S
 ******************************************************************************/

#include <Eigen/Dense>

typedef Eigen::Vector3d vec3d;
typedef Eigen::RowVectorXd vecxd;
typedef Eigen::MatrixXd MatXd;

/*******************************************************************************
 *       N A M E S P A C E
 ******************************************************************************/

namespace DYN{
/*
 * \Struct holding the kinematic parameters
 */
struct Params{
  int l0, l1, l2, l3, l4, l5, l6, l7, off, Xtr, Ytr, Ztr, K;
  vecxd target_pos;
  vecxd q;
  Params() : l0(0), l1(0), l2(0), l3(0), l4(0), l5(0), l6(0), l7(0), off(0), Xtr(0), Ytr(0), Ztr(0), K(0), q(8), target_pos(3){}
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
    static void kukaDynRCM(Params& kinParams);

  }; // class kukaDyn

} //namespace
