/**
 * TSAI 1999 Hand-eye calibration
 * Emanuele Ruffaldi 2016
 *
 * Copyright Scuola Superiore Sant'Anna (2016) 
 * Emanuele Ruffaldi, e.ruffaldi@sssup.it
 */
#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>


template <class T>
 using Vector3 = Eigen::Matrix<T,3,1>;

template <class T>
 using Matrix3 = Eigen::Matrix<T,3,3>;

template <class T>
 using Matrix4 = Eigen::Matrix<T,4,4>;


/**
 * Computes Tsai calibration: four frames
 *	 R = robot root
 * 	 E = robot end-effector
 *   C = camera
 *   M = marker
 *
 * The objective is EC
 * We move E wrt R and keep R and M fixed each other
 */

float calibrationTsai(const std::vector<Eigen::Matrix4f> & cMm,
                     const std::vector<Eigen::Matrix4f> & rMe,
                     Eigen::Matrix4f &eMc, bool computeResidual = false);

double calibrationTsai(const std::vector<Eigen::Matrix4d> & cMm,
                     const std::vector<Eigen::Matrix4d> & rMe,
                     Eigen::Matrix4d &eMc, bool computeResidual = false);

/*
 * These are internal functions used for testing the Tsai
 * They mainly deal with the specific Tsai parametrization of the vector that is more robust to small values than the regular one
 */
template <class T>
Eigen::Matrix<T,3,1> paratsaiprime2paratsai(const Eigen::Matrix<T,3,1> & a);
template <class T>
Eigen::Matrix<T,3,1> paratsai2paratsaiprime(const Eigen::Matrix<T,3,1> & a);
template <class T>
Eigen::Matrix<T,3,3> paratsai2rot(const Eigen::Matrix<T,3,1> & a);
template <class T>
Eigen::Quaternion<T> paratsai2quat(const Eigen::Matrix<T,3,1> & a);
template <class T>
Eigen::Matrix<T,3,1> quat2paratsai(const Eigen::Quaternion<T> & q);
template <class T>
T paratsaiprime2theta(const Eigen::Matrix<T,3,1> & x);
template <class T>
T paratsai2theta(const Eigen::Matrix<T,3,1> & x);


/// returns the skew of the vector
template <class T>
inline Matrix3<T> skewtsai(const Vector3<T> & a)
{
  Matrix3<T> r;
  r  << 0, -a(2), a(1),  a(2),0,-a(0),-a(1),a(0),0;
  return r;
}


template <class T>
Vector3<T> paratsaiprime2paratsai(const Vector3<T> & a)
{
    return 2*a/sqrt(1+a.squaredNorm()); // equation 14
}

template <class T>
Vector3<T> paratsai2paratsaiprime(const Vector3<T> & a)
{
  return a/sqrt(4-a.squaredNorm());
}

/// converts the parametric in tsai form to rotation
/// NO NEED OF TRIGONOMETRICS
/// equation 10
template <class T>
Matrix3<T> paratsai2rot(const Vector3<T> & a)
{
  T sn = a.squaredNorm();
  T alpha = sqrt(4-sn);
  return (1-sn/2)*Matrix3<T>::Identity() + 0.5*(a*a.transpose() + alpha * skewtsai(a));
}

/// converts the parametric in tsai form to rotation
/// NO NEED OF TRIGONOMETRICS
/// equation 9 modified NOT USED
template <class T>
Eigen::Quaternion<T> paratsai2quat(const Vector3<T> & a)
{
  double sinthetahalf = a.norm()/2; // 2 * sin(theta)/2 => sin(theta)/2
  return Eigen::Quaternion<T>(sqrt(1-sinthetahalf*sinthetahalf),a.x()/2,a.y()/2,a.z()/2); // cos(theta/2)
}

template <class T>
T paratsai2theta(const Vector3<T> & x)
{
  return 2*asin(x.norm()/2); // eq. 9    
}
