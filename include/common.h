/**
 * \file      forces.h
 *
 * \brief     C++ Header File template
 *
 * \author    Jean Pierre Richa
 *
 * \Project:  Medical Robotics
 */

#ifndef _COMMON_H_
#define _COMMON_H_

/*******************************************************************************
 *       I N C L U D E - F I L E S
 ******************************************************************************/

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <boost/numeric/odeint.hpp>

typedef Eigen::Vector3d vec3d;
typedef Eigen::VectorXd vecxd;
typedef Eigen::MatrixXd MatXd;
typedef std::vector<double> dvec;
typedef std::vector<Eigen::Vector3d> vec3Eigens;
typedef std::vector<Eigen::VectorXd> vecEigens;
typedef std::vector<Eigen::MatrixXd> vecEigenMats;
typedef std::vector<Eigen::Matrix<double,4,4> > vecEigen4Mat;
typedef std::vector<Eigen::Matrix<double,3,3> > vecEigen3Mat;

template <typename T>
int sign (const T &val){
  return (val>0)-(val<0);
}
template <typename T>
std::vector<int> sign (const std::vector<T> &v){
  std::vector<int> r(v.size());
  std::transform(v.begin(), v.end(), r.begin(), (int(*)(const T&))sign);
  return r;
}
#endif // _COMMON_H_
