/**
 * \file      forces.h
 *
 * \brief     C++ Header File template
 *
 * \author    Jean Pierre Richa
 *
 * \Project:  Medical Robotics
 */

#ifndef _FORCES_
#define _FORCES_

/*******************************************************************************
 *       I N C L U D E - F I L E S
 ******************************************************************************/

#include <Eigen/Dense>

typedef Eigen::Vector3d vec3d;
typedef Eigen::VectorXd vecxd;
typedef Eigen::MatrixXd MatXd;
typedef std::vector<Eigen::VectorXd> vecEigens;

/*******************************************************************************
 *       N A M E S P A C E
 ******************************************************************************/

namespace DYN{
/*
 * \Struct holding the dynamic and forces parameters
 */
struct dynamicParams{
  int N;
  double fv, fc;
  MatXd tau_act, inertiaVectors;
  vecxd fv_vec, lengths, masses;
  dynamicParams() : N(0), fv(0.0), fc(0.0), tau_act(N,1), inertiaVectors(7,6), fv_vec(7), lengths(7), masses(7){}
};
struct inertiaMat{
  MatXd I1, I2, I3, I4, I5, I6, I7;
  inertiaMat() : I1(3,3), I2(3,3), I3(3,3), I4(3,3), I5(3,3), I6(3,3), I7(3,3){}
};
/**
 * \class force
 *
 * \estimates the external forces acting on the kuka LWR
 *
 */

class FORCE{
  public:
    static void external_forces_estimation(dynamicParams& dynParams,
                                           inertiaMat& iMat);
    static MatXd calcInertiaMat(vecxd& v){
      MatXd inertiaMatrix(3,3);
      inertiaMatrix << v[0], v[3], v[4],
                       v[3], v[1], v[5],
                       v[4], v[5], v[2];
      return inertiaMatrix;
    }

  }; // class kukaDyn

} //namespace
#endif // _FORCES_
