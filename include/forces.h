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
#include "common.h"
#include "kukaDynRCM.h"

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

/*
 * \Struct holding the vectors of inertia matrices
 */
struct inertiaMat{
  MatXd I1, I2, I3, I4, I5, I6, I7;
  inertiaMat() : I1(3,3), I2(3,3), I3(3,3), I4(3,3), I5(3,3), I6(3,3), I7(3,3){}
};

/*
 * \Struct holding the residual method parameters
 * \W: angular velocity, Wd: angular acceleration, Wmd: rotor angular accelerations
 * \Pdd: linear acceleration, Pcdd: linear com acceleration
 */
struct residualParams{
  int kri, tau; // Gear Reduction Ratio and forces respectively
  double Im, mm; // Rotor inertia and Rotor mass respectively
  // vecEigen3Mat L; // vector of 3x3 eigen matrices
  vecxd m, qd, qdd, g, L; // g:1x3 vec, vectors (m & I:1x7, qd & qdd:7x1(qd and qdd should be transposed))
  vec3Eigens W, Wd, Wmd, Pdd, Pcdd, Rici, Rij; // all = 7 dimensional vector (1x7 vector) of 3x1 vectors
  vecEigen4Mat A; // 8 dimensional vector of 4x4 matrices
  vecEigen3Mat Rot, I; // 8 dimensional vector of 3x3 matrices
  residualParams() : tau(0), kri(0), L(7), Im(0.0), mm(0.0), m(7), qd(7), qdd(7), g(3), W(7), Wd(7), Wmd(7), Pdd(7), Pcdd(7), Rici(7), Rij(7){}
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
    static MatXd calcSkewSymm(vec3d& vec){
      MatXd skewSymmMat(3,3);
      skewSymmMat.setZero();
      skewSymmMat << 0,      -vec[2],  vec[1],
                     vec[2],       0, -vec[0],
                    -vec[1],  vec[0],       0;

      return skewSymmMat;
    }

  }; // class kukaDyn

class RESIDUAL{
  public:
    static void NewEulForw(Params& kinParams, dynamicParams& dynParams, residualParams& residual);
    static void NewEul_Aux(Params& kinParams, dynamicParams& dynParams, residualParams& residual);
    static MatXd KukaKinematics(Params& kinParams, dynamicParams& dynParams, residualParams& residual);
    static int NewEulBack(Params& kinParams, dynamicParams& dynParams, residualParams& residual);
};

} //namespace
#endif // _FORCES_
