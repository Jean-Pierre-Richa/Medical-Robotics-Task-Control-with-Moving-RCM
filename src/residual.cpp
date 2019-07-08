#include <iostream>
#include "common.h"
#include "forces.h"
#include "kukaDynRCM.h"

using namespace std;

DYN::inertiaMat inertiaMat;

void DYN::RESIDUAL::NewEul_Aux(Params& kinParams, dynamicParams& dynParams, residualParams& residual){

  NewEulForw(kinParams, dynParams, residual);

  residual.tau = NewEulBack(kinParams, dynParams, residual);

}

void DYN::RESIDUAL::NewEulForw(Params& kinParams, dynamicParams& dynParams, residualParams& residual){

  MatXd T(4,4);
  T.setZero();
  T = KukaKinematics(kinParams, dynParams, residual);

  vec3d prev_w(0,0,0);
  vec3d prev_wd(0,0,0);
  vec3d prev_pdd(0,0,0);

  vec3d curr_w(0,0,0);
  vec3d curr_wd(0,0,0);
  vec3d curr_pdd(0,0,0);

  vec3d curr_pcdd(0,0,0);
  vec3d curr_wmd(0,0,0);

   // Forward phase
   vec3d z0(0, 0, 1);
   for (int i=0; i<dynParams.N; i++){
     /*
      * \Vector from origin of frame i to the center of mass of link i
      * \(negative quantity siciliano's book)
      */
     vec3d rici(-residual.L(i)/2, 0, 0);
     /*
      * \Vector from the origin of frame j-1 (frame of link i) to the origin of
      * \frame j (j=i+1=i) from siciliano's book
      */
     MatXd a(4,4);
     a.setZero();
     a = residual.A[i];
     // cout << "a" << a << endl;
     vec3d rij(0,0,0);
     rij = a.block<3,1>(0,3);
     // cout << "rij" << rij << endl;
     /*
      * \R from frame i to frame i-1 (1->0, 2->1,...) we'll use its inverse
      */
     MatXd R(3,3);
     R.setZero();
     R = residual.Rot[i];
     // cout << "R" << R << endl;
     /*
      * \Angular velocity of link i w.r.t frame i. (Equation 7.107 (Siciliano))
      */
     MatXd Rt(3,3);
     Rt.setZero();
     Rt = R.transpose();
     // angular velocity of link i expressed in frame i
     curr_w = Rt*(prev_w+residual.qd(i)*z0);
     /*
      * \Angular acceleration w.r.t frame i. (Equation 7.108 (Siciliano))
      */
     // angular acceleration of link i expressed in frame i
     curr_wd = Rt*(prev_wd+residual.qdd(i)*z0+residual.qd(i)*prev_w.cross(z0));
     /*
      * \Linear acceleration w.r.t frame i (Equation 7.109 (Siciliano))
      */
     curr_pdd = Rt*prev_pdd+curr_wd.cross(rij)+curr_w.cross(curr_w.cross(rij));
     /*
      * \linear acceleration of link i expressed in frame i (Equation 7.110(Siciliano))
      */
     curr_pcdd = curr_pdd+curr_w.cross(rici)+curr_w.cross(curr_w.cross(rici));
     /*
      * \Angular acceleration of the rotor i w.r.t frame i-1 (Equation 7.111 (Siciliano))
      */
     curr_wmd = prev_wd + residual.kri*residual.qdd(i)*z0+residual.kri*residual.qd(i)*prev_w.cross(z0);
     /*
      * \List of velocities and accelerations to pass to the backward step
      */
     residual.Rot[i] = R;
     residual.W[i] = curr_w;
     residual.Wd[i] = curr_wd;
     residual.Pdd[i] = curr_pdd;
     residual.Pcdd[i] = curr_pcdd;
     residual.Wmd[i] = curr_wmd;
     residual.Rij[i] = rij;
     residual.Rici[i] = rici;
     prev_w = curr_w;
     prev_wd = curr_wd;
     prev_pdd = curr_pdd;

   }

}

int DYN::RESIDUAL::NewEulBack(Params& kinParams, dynamicParams& dynParams, residualParams& residual){


  int ta = 0;
  // Backwards
  vec3d z0(0,0,1);
  vec3d prev_f(0,0,0);
  vec3d prev_u(0,0,0);
  vec3d curr_f(0,0,0);
  vec3d curr_u(0,0,0);

  vec3d zm;
  zm.setZero();
  zm=z0;
  cout << "zm " << zm << endl;

  MatXd tau(dynParams.N,1);
  tau.setZero();

  cout << "tau " << tau << endl;

  // MatXd constructS = DYN::FORCE::calcSkewSymm(vec);


  vecxd qd(dynParams.N+1);
  qd.setZero();
  vecxd qdd(dynParams.N+1);
  qdd.setZero();


  for (int i=dynParams.N; i>1; i--){
    // Rotation matrix
    MatXd R(3,3);
    R.setZero();
    R = residual.Rot[i];
    // Rotation matrix transpose
    MatXd Rt(3,3);
    Rt.setZero();
    Rt = R.transpose();
    // Linear acceleration
    vec3d pcdd(0,0,0);
    pcdd = residual.Pcdd[i];

    vec3d rij(0,0,0);
    rij = residual.Rij[i];

    vec3d rici(0,0,0);
    rici = residual.Rici[i];
    // Angular velocity
    vec3d w(0,0,0);
    w = residual.W[i];
    // Angular acceleration
    vec3d wd(0,0,0);
    wd = residual.Wd[i];
    // Rotor angular acceleration
    vec3d wmd(0,0,0);
    wmd = residual.Wmd[i];
    /*
     * \Force Equation (siciliano 7.112)
     */
    // curr_f = R*prev_f + residual.m[i]  *pcdd;
    /*
     * \get the skew symmetric matrix
     */
    MatXd Skew(3,3);
    Skew.setZero();
    Skew = DYN::FORCE::calcSkewSymm(rici);
    /*
     * \Link inertia tensor. Steiner Theorem, equation 7.60 (siciliano)
     */
    MatXd I_tensor(3,3);
    I_tensor.setZero();
    // I_tensor = MatXd::Identity(3,3)*residual.I[i] + m[i] * ((Skew.transpose())*Skew);

    // cout << I_tensor << endl;

  }
  return ta;
}

MatXd DYN::RESIDUAL::KukaKinematics(Params& kinParams, dynamicParams& dynParams, residualParams& residual){

  vecxd q(8);
  q = kinParams.q;

  double q1=q(0), q2=q(1), q3=q(2), q4=q(3), q5=q(4), q6=q(5), q7=q(6), lambda=q(7);
  double s1=sin(q1), s2=sin(q2), s3=sin(q3), s4=sin(q4), s5=sin(q5), s6=sin(q6), s7=sin(q7);
  double c1=cos(q1), c2=cos(q2), c3=cos(q3), c4=cos(q4), c5=cos(q5), c6=cos(q6), c7=cos(q7);

  double d1, d3, d5, d7;
  d1 = residual.L(0) + residual.L(1);
  d3 = residual.L(2) + residual.L(3);
  d5 = residual.L(4) + residual.L(5);
  d7 = residual.L(6) + kinParams.ee_length;

  MatXd a1(4,4);
  a1 << c1, 0,  s1, 0,
        s1, 0, -c1, 0,
        0,  1,  0,  d1,
        0,  0,  0,  1;

  MatXd a2(4,4);
  a2 << c2, 0,  s2, 0,
        s2, 0, -c2, 0,
        0,  1,  0,  0,
        0,  0,  0,  1;

  MatXd a3(4,4);
  a3 << c3, 0,  s3, 0,
        s3, 0, -c3, 0,
        0,  1,  0,  d3,
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
  a6 << c6, 0,  s6, 0,
        s6, 0, -c6, 0,
        0,  1,  0,  0,
        0,  0,  0,  1;

  MatXd a7(4,4);
  a7 << c7,  s7, 0,  0,
        s7, -c7, 0,  0,
        0,   1,  0,  d7,
        0,   0,  0,  1;

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

  MatXd r1(3,3);
  r1 << c1, 0,  s1,
        s1, 0, -c1,
        0,  1,   0;

  MatXd r2(3,3);
  r2 << c2, 0,  s2,
        s2, 0, -c2,
        0,  1,   0;

  MatXd r3(3,3);
  r3 << c3, 0,  s3,
        s3, 0, -c3,
        0,  1,   0;

  MatXd r4(3,3);
  r4 << c4, 0,  s4,
        s4, 0, -c4,
        0,  1,   0;

  MatXd r5(3,3);
  r5 << c5, 0,  s5,
        s5, 0, -c5,
        0,  1,   0;

  MatXd r6(3,3);
  r6 << c6, 0,  s6,
        s6, 0, -c6,
        0,  1,   0;

  MatXd r7(3,3);
  r7 << c7,  s7, 0,
        s7, -c7, 0,
        0,    1, 0;

  residual.Rot.push_back(r1);
  residual.Rot.push_back(r2);
  residual.Rot.push_back(r3);
  residual.Rot.push_back(r4);
  residual.Rot.push_back(r5);
  residual.Rot.push_back(r6);
  residual.Rot.push_back(r7);

  residual.A.push_back(a1);
  residual.A.push_back(a2);
  residual.A.push_back(a3);
  residual.A.push_back(a4);
  residual.A.push_back(a5);
  residual.A.push_back(a6);
  residual.A.push_back(a7);

  cout << "r1 " << r1 << endl;
  // cout << residual.Rot << endl;

  return kwf7;

}
