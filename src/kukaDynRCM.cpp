/**
 * \file      RCM.cpp
 *
 * \brief     C++ Source File template
 *
 * \author    Jean Pierre Richa
 *
 * Project:   Medical Robotics
 */

/*******************************************************************************
 *       I N C L U D E - F I L E S
 ******************************************************************************/

#include <iostream>
#include <vector>
#include <Eigen/QR>
#include <Eigen/Dense>
#include "kukaDynRCM.h"
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

/*******************************************************************************
 *       C L A S S - M E M B E R - F U N C T I O N S
 ******************************************************************************/

// vec3Eigens DYN::RCM::dirKin(const DYN::Params& kinParams){
//
//   vecxd q(8);
//   q = kinParams.q;
//
//   double q0=q[0], q1=q[1], q2=q[2], q3=q[3], q4=q[4], q5=q[5], q6=q[6], lambda=q[7];
//   double s1=sin(q0), s2=sin(q1), s3=sin(q2), s4=sin(q3), s5=sin(q4), s6=sin(q5), s7=sin(q6);
//   double c1=cos(q0), c2=cos(q1), c3=cos(q2), c4=cos(q3), c5=cos(q4), c6=cos(q5), c7=cos(q6);
//
//   // double l0=kinParams.l0, l1=kinParams.l1,
//   //        l2=kinParams.l2, l3=kinParams.l3,
//   //        l4=kinParams.l4, l5=kinParams.l5,
     //     // l7=kinParams.l7,
//   //        off=kinParams.off;
//   // cout << "l7 " << l7 << endl;
//
//
//
//
// }


// vecxd DYN::RCM::kukaDynRCM(DYN::Params& kinParams){
vecxd DYN::RCM::kukaDynRCM(const vecxd& q, DYN::Params& kinParams){

  // cout << "kinParams.q " << kinParams.q[3] << endl;

  // vecxd q(8);
  // q = kinParams.q;


  double q1=q[0], q2=q[1], q3=q[2], q4=q[3], q5=q[4], q6=q[5], q7=q[6], lambda=q[7];
  double s1=sin(q1), s2=sin(q2), s3=sin(q3), s4=sin(q4), s5=sin(q5), s6=sin(q6), s7=sin(q7);
  double c1=cos(q1), c2=cos(q2), c3=cos(q3), c4=cos(q4), c5=cos(q5), c6=cos(q6), c7=cos(q7);
  //
  double l0=kinParams.l0, l1=kinParams.l1,
         l2=kinParams.l2, l3=kinParams.l3,
         l4=kinParams.l4, l5=kinParams.l5,
         l7=kinParams.l7,
         off=kinParams.off;

  double ee_length=0.23;
  double d1=0.0;
  double d3=0.4;
  double d5=0.39;
  double d7=0.078 + ee_length;

  // Transformations matrices
  MatXd a1(4,4);
  a1 << c1, 0,  s1, 0,
        s1, 0, -c1, 0,
        0,  1,  0,  d1,
        0,  0,  0,  1;

  MatXd a2(4,4);
  a2 << c2,  0, -s2,  0,
        s2,  0,  c2,  0,
        0,  -1,   0,  0,
        0,   0,   0,  1;

  MatXd a3(4,4);
  a3 << c3, 0,  -s3, 0,
        s3, 0, c3, 0,
        0,  -1,  0,  d3,
        0,  0,  0,  1;

  MatXd a4(4,4);
  a4 << c4, 0,  s4, 0,
        s4, 0, -c4, 0,
        0,  1,  0,  0,
        0,  0,  0,  1;

  MatXd a5(4,4);
  a5 << c5, 0,  s5, 0,
        s5, 0, -c5, 0,
        0,  1,  0,  d5,
        0,  0,  0,  1;

  MatXd a6(4,4);
  a6 << c6, 0,  -s6, 0,
        s6, 0, c6, 0,
        0,  -1,  0,  0,
        0,  0,  0,  1;

  MatXd a7(4,4);
  a7 << c7,  -s7,  0,   0,
        s7,   c7,  0,   0,
        0,     0,  1,  d7,
        0,     0,  0,   1;

  MatXd kwf1(4,4);
  MatXd kwf2(4,4);
  MatXd kwf3(4,4);
  MatXd kwf4(4,4);
  MatXd kwf5(4,4);
  MatXd kwf6(4,4);
  MatXd kwf7(4,4);

  kwf1 = a1;
  kwf2 = kwf1*a2;
  kwf3 = kwf2*a3;
  kwf4 = kwf3*a4;
  kwf5 = kwf4*a5;
  kwf6 = kwf5*a6;
  kwf7 = kwf6*a7;

  vec3d p1 = kwf1.col(3).head(3);
  vec3d p2 = kwf2.col(3).head(3);
  vec3d p3 = kwf3.col(3).head(3);
  vec3d p4 = kwf4.col(3).head(3);
  vec3d p5 = kwf5.col(3).head(3);
  vec3d p6 = kwf6.col(3).head(3);
  vec3d p7 = kwf7.col(3).head(3);

  // Rotation matrices
  MatXd r1(3,3);
  r1 = a1.block<3,3>(0,0);

  MatXd r2(3,3);
  r2 = a2.block<3,3>(0,0);

  MatXd r3(3,3);
  r3 = a3.block<3,3>(0,0);

  MatXd r4(3,3);
  r4 = a4.block<3,3>(0,0);

  MatXd r5(3,3);
  r5 = a5.block<3,3>(0,0);

  MatXd r6(3,3);
  r6 = a6.block<3,3>(0,0);

  MatXd r7(3,3);
  r7 = a7.block<3,3>(0,0);

  MatXd rwf1(3,3);
  MatXd rwf2(3,3);
  MatXd rwf3(3,3);
  MatXd rwf4(3,3);
  MatXd rwf5(3,3);
  MatXd rwf6(3,3);
  MatXd rwf7(3,3);

  rwf1 = r1;
  rwf2 = rwf1*r2;
  rwf3 = rwf2*r3;
  rwf4 = rwf3*r4;
  rwf5 = rwf4*r5;
  rwf6 = rwf5*r6;
  rwf7 = rwf6*r7;

  vec3d stat(0,0,1);
  vec3d z_0(0,0,1);
  vec3d z_1 = rwf1*stat;
  vec3d z_2 = rwf2*stat;
  vec3d z_3 = rwf3*stat;
  vec3d z_4 = rwf4*stat;
  vec3d z_5 = rwf5*stat;

  // considering that the end effector is at joint 6

  vec3d p_1 = kwf1.col(3).head(3);
  vec3d p_2 = kwf2.col(3).head(3);
  vec3d p_3 = kwf3.col(3).head(3);
  vec3d p_4 = kwf4.col(3).head(3);
  vec3d p_5 = kwf5.col(3).head(3);
  vec3d p_6_e = kwf6.col(3).head(3);

  vec3d p_6_e_0 = p_6_e;
  vec3d p_6_e_1 = p_6_e - p_1;
  vec3d p_6_e_2 = p_6_e - p_2;
  vec3d p_6_e_3 = p_6_e - p_3;
  vec3d p_6_e_4 = p_6_e - p_4;
  vec3d p_6_e_5 = p_6_e - p_5;

  vec3d linJ6vec6(0,0,0);

  MatXd linJ6(3,7);
  linJ6.setZero();
  linJ6.col(0) = z_0.cross(p_6_e_0);
  linJ6.col(1) = z_1.cross(p_6_e_1);
  linJ6.col(2) = z_2.cross(p_6_e_2);
  linJ6.col(3) = z_3.cross(p_6_e_3);
  linJ6.col(4) = z_4.cross(p_6_e_4);
  linJ6.col(5) = z_5.cross(p_6_e_5);
  linJ6.col(6) = linJ6vec6;


  vec3d z_6;
  z_6 = rwf6*stat;
  vec3d p_6 = kwf6.col(3).head(3);
  vec3d p_7_e = kwf7.col(3).head(3);

  vec3d p_7_e_0 = p_7_e;
  vec3d p_7_e_1 = p_7_e - p_1;
  vec3d p_7_e_2 = p_7_e - p_2;
  vec3d p_7_e_3 = p_7_e - p_3;
  vec3d p_7_e_4 = p_7_e - p_4;
  vec3d p_7_e_5 = p_7_e - p_5;
  vec3d p_7_e_6 = p_7_e - p_6;

  MatXd linJ7(3,7);
  linJ7.setZero();
  linJ7.col(0) = z_0.cross(p_7_e_0);
  linJ7.col(1) = z_1.cross(p_7_e_1);
  linJ7.col(2) = z_2.cross(p_7_e_2);
  linJ7.col(3) = z_3.cross(p_7_e_3);
  linJ7.col(4) = z_4.cross(p_7_e_4);
  linJ7.col(5) = z_5.cross(p_7_e_5);
  linJ7.col(6) = z_6.cross(p_7_e_5);

  // RCM Kinematics

  // Jacobian RCM
  MatXd J_lambda(3,7); // 3x7 matrix
  J_lambda = (1-lambda)*linJ7+lambda*linJ6; // 3x7 matrix
  // cout << "\n Jlambda " << J_lambda << "\n" << endl;
  vec3d g = p6-p7; // 3D vector
  MatXd J_rCM(3,8); // 3x8 matrix composed by the matrix J_lambda and vector g
  J_rCM.setZero();
  J_rCM.block<3,7>(0,0)=J_lambda; // First 3x7 block
  J_rCM.block<3,1>(0,7)=g; // 3x1 vector on the 8th column

  // Task Jacobian
  MatXd Jt(3,8); // 3x8 matrix
  Jt.setZero();
  Jt.block<3,7>(0,0)=linJ6; // First 3x7 block

  // Extention Jacobian 6x8 it includes the J_RCM and Jt (the task Jacobian)
  MatXd J_a(6,8);
  J_a.setZero();
  J_a.block<3,8>(0,0)=J_rCM;
  J_a.block<3,8>(3,0)=Jt;
  //
  // cout << "J_rCM " << J_rCM <<"\n"<< endl;
  // cout << "Jt " << Jt <<"\n"<< endl;
  // cout << "linJ6 " << linJ6 <<"\n"<< endl;
  // cout << "linJ7 " << linJ7 <<"\n"<< endl;

  // RCM and additional task errors
  double x_rcm=p7[0]+lambda*(p6[0]-p7[0]);
  double y_rcm=p7[1]+lambda*(p6[1]-p7[1]);
  double z_rcm=p7[2]+lambda*(p6[2]-p7[2]);

  vec3d err_lam(kinParams.target_pos[0]-x_rcm, kinParams.target_pos[1]-y_rcm, kinParams.target_pos[2]-z_rcm);
  vec3d err_t(kinParams.target_pos[0]-p6[0], kinParams.target_pos[1]-p6[1], kinParams.target_pos[2]+ee_length-p6[2]);

  MatXd err(6,1); // 6x1
  err.block<3,1>(0,0)=err_lam;
  err.block<3,1>(3,0)=err_t;
  // cout << "err " << err << endl;

  // Position and orientation control with RCM constraint
  MatXd identity = MatXd::Identity(8,8);
  // Jacobian pseudoinverse
  MatXd Jp(8,6);
  Jp.setZero();
  // Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cqr(J_a);
  // Jp = cqr.pseudoInverse();
  Jp = J_a.completeOrthogonalDecomposition().pseudoInverse(); // 8x6 matrix
  // cout << "JP " << Jp << endl;
  // cout << "K " << kinParams.k << endl;

  vecxd w(8);
  w << 0,0,0,0,0,0,0,lambda-0.2;
  w=-10*w;
  // cout << "W " << w << endl;
  // cout << "J_a " << J_a << endl;

  vecxd u(8);
  u.setZero();
  // // 8x6 * 6x6 * 6x1 + (8x8-(6x8*8x6)) * 8x1 // dimensions to be checked
  u=Jp*kinParams.k*err+(identity-Jp*J_a)*w;

  vecxd q_dot(8);
  q_dot.setZero();

  // cout << "u" << u << endl;

  q_dot << u[0], u[1], u[2], u[3], u[4], u[5], u[6], u[7];

  // cout << "q_dot " << q_dot << endl;

  return q_dot;

}
