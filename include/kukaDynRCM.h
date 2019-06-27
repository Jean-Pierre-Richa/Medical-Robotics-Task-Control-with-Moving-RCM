
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

#include <vector>
#include <Eigen/Dense>

typedef Eigen::Vector3d vec3d;
typedef Eigen::RowVectorXd vecxd;
typedef Eigen::MatrixXd MatXd;

/*******************************************************************************
 *       N A M E S P A C E
 ******************************************************************************/

namespace DYN{

/**
 * \class kukaDyn
 *
 * \Calculates the direct kinematics and jacobian matrices of the kuka LWR for
 * \the 6th and 7th joints
 *
*/

class RCM{
  public:
    static void kukaDynRCM(vecxd& q);

  }; // class kukaDyn

} //namespace
