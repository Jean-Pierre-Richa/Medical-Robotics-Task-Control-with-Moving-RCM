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

typedef Eigen::Vector3d vec3d;
typedef Eigen::VectorXd vecxd;
typedef Eigen::MatrixXd MatXd;
typedef std::vector<Eigen::Matrix<double,4,4> > vecEigen4Mat;
typedef std::vector<Eigen::Matrix<double,3,3> > vecEigen3Mat;
typedef std::vector<Eigen::Vector3d> vec3Eigens;
typedef std::vector<Eigen::MatrixXd> vecEigenMats;
typedef std::vector<Eigen::VectorXd> vecEigens;

#endif // _COMMON_H_
