
/* 
 * Tsai algorithm 1989, initially ported to Eigen from VISP and then manually modified for efficiency, clearness and more Eigenability
 * 
 * Copyright Scuola Superiore Sant'Anna (2016) 
 * Emanuele Ruffaldi, e.ruffaldi@sssup.it
 *
 * TODO: exclude collinear pairs
 * TODO: compute residuals for providing a measure of the error 
 * 
 * Structure optimization (double the memory but faster)
 *    A1,A2,b1,b2,tmp
 *    loop all pairs
 *      exclude bad pairs
 *      write A1 and A2
 *      write b1 and partial b2 the part that doesn't require Rcg
 *    solve Rcg using A1,b1
 *    update b2 using eRcg*cijo.matrix().template block<3,1>(0,3) (so we need to store cijo translation part into tmp)
 *    solve translation for A2,b2
 */
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include "tsai.hpp"



/// extracts the translation part
template <class T>
static Vector3<T> gettx(const Matrix4<T> & x)
{
  return x.template block<3,1>(0,3);
}

/// extracts the translation part
template <class T>
static Vector3<T> gettx(const Eigen::Transform<T,3,Eigen::Affine> & x)
{
  return x.matrix().template block<3,1>(0,3);
}

/// extracts rotation
template <class T>
static Matrix3<T> getrot(const Eigen::Transform<T,3,Eigen::Affine> & x)
{
  return x.matrix().template block<3,3>(0,0);
}

/// extracts rotation
template <class T>
static Matrix3<T> getrot(const Matrix4<T> & x)
{
  return x.template block<3,3>(0,0);
}

/// extracts quat and transform from matrix 4x4
template <class T>
static void extract(const Matrix4<T> & x,Eigen::Quaternion<T> & q, Vector3<T> & t)
{
  q = x.template block<3,3>(0,0);
  t = x.template block<3,1>(0,3);
}

/// extract rotation
template <class T>
static void extract(const Matrix4<T> & x,Eigen::Quaternion<T> & q)
{
  q = x.template block<3,3>(0,0);
}

/// makes affine3f from matrix using extraction
template <class T>
static void extract(const Matrix4<T> & x,Eigen::Transform<T,3,Eigen::Affine> & a)
{
  auto q = x.template block<3,3>(0,0);
  auto t = x.template block<3,1>(0,3);
  a.fromPositionOrientationScale(t,q,Vector3<T>(1,1,1));
}




template <class T>
static inline T sinc(T sx,T x)
{
  return fabs(x) < 1e-8 ? 1 : sx/x;
}

template <class T>
static inline T sinc(T x)
{
  return fabs(x) < 1e-8 ? 1 : sin(x)/x;
}

// modified tsai is: 2 sin(theta/2) axis
// BUT quaternion real part is:  sin(theta/2) * axis
// equation 9
// THIS IS FROM VISP
template <class T>
Vector3<T> quat2paratsai(const Eigen::Quaternion<T> & q)
{
  //Eigen::AngleAxisf aa(q);
  //return 2*sin(aa.angle()/2)*aa.axis();
  //return 2*Vector3<T>(q.x(),q.y(),q.z()); 
  const double minimum = 0.0001;
  Matrix3<T> R(q);
  Vector3<T> r;

  double s = (R(1,0)-R(0,1))*(R(1,0)-R(0,1))
    + (R(2,0)-R(0,2))*(R(2,0)-R(0,2))
    + (R(2,1)-R(1,2))*(R(2,1)-R(1,2));
   s = sqrt(s)/2.0;
   double c = (R(0,0)+R(1,1)+R(2,2)-1.0)/2.0;
   double theta=atan2(s,c);  /* theta in [0, PI] since s > 0 */

  // General case when theta != pi. If theta=pi, c=-1
  if ( (1+c) > minimum) // Since -1 <= c <= 1, no fabs(1+c) is required
  {
    double si = sinc(s,theta);

    r[0] = (R(2,1)-R(1,2))/(2*si);
    r[1] = (R(0,2)-R(2,0))/(2*si);
    r[2] = (R(1,0)-R(0,1))/(2*si);
  }
  else /* theta near PI */
  {
    if ( (R(0,0)-c) < std::numeric_limits<double>::epsilon() )
      r[0] = 0.;
    else
      r[0] = theta*(sqrt((R(0,0)-c)/(1-c)));
    if ((R(2,1)-R(1,2)) < 0) r[0] = -r[0];

    if ( (R(1,1)-c) < std::numeric_limits<double>::epsilon() )
      r[1] = 0.;
    else
      r[1] = theta*(sqrt((R(1,1)-c)/(1-c)));

    if ((R(0,2)-R(2,0)) < 0) r[1] = -r[1];

    if ( (R(2,2)-c) < std::numeric_limits<double>::epsilon() )
      r[2] = 0.;
    else
      r[2] = theta*(sqrt((R(2,2)-c)/(1-c)));

    if ((R(1,0)-R(0,1)) < 0) r[2] = -r[2];
  }
  return r;

}


template <class T>
T paratsaiprime2theta(const Vector3<T> & x)
{
  return 2*atan(x.norm()); // eq. 13
}


// port and improve: https://github.com/thomas-moulard/visp-deb/blob/master/src/camera/calibration/vpCalibrationTools.cpp
template<class T>
T calibrationTsaiT(const std::vector<Matrix4<T> > & cMo,const std::vector<Matrix4<T> > & rMe, Matrix4<T> &eMc, bool computeResidual)
{
  assert(cMo.size() == rMe.size() && !cMo.empty() && "not empty calibration matrices");
  const unsigned int nbPose = cMo.size();
  Vector3<T> x; // tsaipara result
  Eigen::Matrix<T,Eigen::Dynamic,3> A((nbPose*(nbPose-1)/2)*3,3); // reused between the two loops over all pairs
  Eigen::Matrix<T,Eigen::Dynamic,1> B(A.rows(),1);
  {
    unsigned int k = 0 ;
    // for all couples ij
    for (unsigned int i=0 ; i < nbPose ; i++)
    {
      Eigen::Quaternion<T> eiR,ioR;
      extract(rMe[i],eiR);
      extract(cMo[i],ioR);

      for (unsigned int j=i+1; j < nbPose ; j++, k+= 3)
      {
          Eigen::Quaternion<T> ejR,joR;
          extract(rMe[j],ejR);
          extract(cMo[j],joR);

          Eigen::Quaternion<T> rRgij = ejR.conjugate() * eiR;
          Eigen::Quaternion<T>  cRijo = joR * ioR.conjugate();

          // CODE: theta = sqrt(|rotationalpart of rRgij|) then multiply rotation part by sinc(theta/2)
          // going from norm theta to scaled by sinc

          // rotation axis (eq 11a,11b,11c)
          Vector3<T>  rPgij = quat2paratsai<T>(rRgij);
          Vector3<T>  cPijo = quat2paratsai<T>(cRijo);

          // TODO: in paper we can remove pair if not good

          // skewtsai(Pgij+Pcij) Pcg' = Pcij - Pgij  equation 12
          A.template block<3,3>(k,0) = skewtsai(Vector3<T>(rPgij + cPijo)) ;
          B.template block<3,1>(k,0) = cPijo - rPgij;       

#if 0
          if(k == 0)
          {            
            std::cout << "rPeij0 R\n" << Matrix3<T>(rRgij) << std::endl;
            std::cout << "rPeij0\n" << rPgij.transpose() << std::endl;
            std::cout << "rPeij0\n" << Eigen::AngleAxisf(rRgij).angle() << "|" << Eigen::AngleAxisf(rRgij).axis() << std::endl;
            std::cout << "cijPo\n" << cPijo.transpose() << std::endl;

            std::cout << "As0\n " << A.template block<3,3>(0,0) << " \n b0" << B.template block<3,1>(0,0).transpose() << std::  endl;
          }  
#endif
      }
    }

    // the output is in the MODIFIED form
    x = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B);

    // convert output to PARAMETRIC eq. 14
    // Pcg' = 1/(4-|Pcg|^2) Pcg
    // Pcg  = 2 Pcg' /(1+|Pcg|^2)
    x = paratsaiprime2paratsai(x);

    // verification of the quality of the vector
    const double theta = paratsai2theta(x); 

    if (std::fabs(theta) < std::numeric_limits<double>::epsilon())
      x.setZero();
    // x *= theta/(2*sin(theta/2));     // paratsai => para NOT NEEDED we work in paratsai
  }

  // Building of the rotation matrix eRc using eq. 10
  Matrix3<T>  eRcg = paratsai2rot(x);

  {
    // Building of the system for the translation estimation
    // for all couples ij
    unsigned int k = 0 ;
    for (unsigned int i = 0 ; i < nbPose ; i++)
    {
      /*
      Eigen::Affine3f eiA,ioA;
      extract(rMe[i],eiA);
      extract(cMo[i],ioA);
      */

      Vector3<T>  eiT = gettx(rMe[i]);
      Vector3<T>  ioT = gettx(cMo[i]);
      Matrix3<T>  eiR = getrot(rMe[i]);


      for (unsigned int j = i+1 ; j < nbPose ; j++, k+= 3)
      {
          /*
          Eigen::Affine3f ejA,joA;
          extract(rMe[j],rejA);
          extract(cMo[j],joA);
          Eigen::Affine3f jiA = ejA.inverse() * eiA; // gij of paper
          */
          Vector3<T>  ejT = gettx(rMe[j]);
          Vector3<T> joT = gettx(cMo[j]);
          Matrix3<T> ejR = getrot(rMe[j]);
          Matrix3<T> jiR = ejR.transpose()*eiR;

          // TODO: in paper we can remove pair if not good COLLINEARITY ISSUE
        
          // equation 15: (Rgij-I) T = Rcg Tcij - Tgij
          A.template block<3,3>(k,0) = jiR-Matrix3<T>::Identity();
          B.template block<3,1>(k,0) = eRcg*joT - jiR*eRcg*ioT + ejR.transpose()*(ejT-eiT);
      }
    }

    eMc.setIdentity();
    eMc.template block<3,3>(0,0) = eRcg;
    eMc.template block<3,1>(0,3) = A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B);

    // A (n*3,3) x(3,1) - B(n*3,1)
    return computeResidual ? (A*eMc.template block<3,1>(0,3)-B).squaredNorm()/A.rows() : 0; // this is the variance of result
 }
}

float calibrationTsai(const std::vector<Eigen::Matrix4f> & cMo,
                      const std::vector<Eigen::Matrix4f> & rMe,
                      Eigen::Matrix4f &eMc,
                      bool computeResidual)
{
    return calibrationTsaiT<float>(cMo,rMe,eMc,computeResidual);
}

double calibrationTsai(const std::vector<Eigen::Matrix4d> & cMo,
                       const std::vector<Eigen::Matrix4d> & rMe,
                       Eigen::Matrix4d &eMc,bool computeResidual)
{
    return calibrationTsaiT<double>(cMo,rMe,eMc,computeResidual);
}
