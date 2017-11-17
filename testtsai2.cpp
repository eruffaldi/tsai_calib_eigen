#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <fstream>
#include "tsai.hpp"

#ifdef USEVISP
#include <visp/vpCalibration.h>
#include <visp/vpMath.h>
#include <visp/vpPose.h>
#include <visp/vpPixelMeterConversion.h>
#endif
//using LTD = SE3dist<float>;
using LTD = float;

#ifdef USEVISP


void assign(Eigen::Matrix4f & a, const 	vpHomogeneousMatrix & b)
{
	vpQuaternionVector q;
	vpTranslationVector t;
	b.extract(q);
	b.extract(t);
	Eigen::Vector3f et(t[0],t[1],t[2]);
	Eigen::Quaternionf eq(q.w(),q.x(),q.y(),q.z());
	a.setIdentity();
	a.block<3,3>(0,0) = Eigen::Matrix3f(eq);
	a.block<3,1>(0,3) = et;
}

void assign(vpHomogeneousMatrix & b,const Eigen::Matrix4f & a)
{
	Eigen::Quaternionf eq(Eigen::Matrix3f(a.block<3,3>(0,0)));
	Eigen::Vector3f et = a.block<3,1>(0,3);

	vpQuaternionVector q(eq.x(),eq.y(),eq.z(),eq.w());
	vpTranslationVector t(et.x(),et.y(),et.z());
	b.buildFrom(t,q);
}



void calibrationTsai(std::vector<vpHomogeneousMatrix>& cMo,
						std::vector<vpHomogeneousMatrix>& rMe,
						vpHomogeneousMatrix &eMc){

  vpColVector x ;
  unsigned int nbPose = (unsigned int)cMo.size();
  if(cMo.size()!=rMe.size()) throw vpCalibrationException(vpCalibrationException::dimensionError,"cMo and rMe have different sizes");
  {
    vpMatrix A ;
    vpColVector B ;
    unsigned int k = 0 ;	
    // for all couples ij
    for (unsigned int i=0 ; i < nbPose ; i++)
    {
      vpRotationMatrix rRei, ciRo ;
      rMe[i].extract(rRei) ;
      cMo[i].extract(ciRo) ;
      //std::cout << "rMei: " << std::endl << rMe[i] << std::endl;

      for (unsigned int j=i+1; j < nbPose ; j++)
      {
        {
          vpRotationMatrix rRej, cjRo ;
          rMe[j].extract(rRej) ;
          cMo[j].extract(cjRo) ;
	  //std::cout << "rMej: " << std::endl << rMe[j] << std::endl;

          vpRotationMatrix rReij = rRej.t() * rRei;

          vpRotationMatrix cijRo = cjRo * ciRo.t();

          vpThetaUVector rPeij(rReij);

          double theta = sqrt(rPeij[0]*rPeij[0] + rPeij[1]*rPeij[1]
                              + rPeij[2]*rPeij[2]);

          for (unsigned int m=0;m<3;m++) rPeij[m] = rPeij[m] * vpMath::sinc(theta/2);

          vpThetaUVector cijPo(cijRo) ;
          theta = sqrt(cijPo[0]*cijPo[0] + cijPo[1]*cijPo[1]
                       + cijPo[2]*cijPo[2]);
          for (unsigned int m=0;m<3;m++) cijPo[m] = cijPo[m] * vpMath::sinc(theta/2);

          vpMatrix As;
          vpColVector b(3) ;

          As = vpColVector::skew(vpColVector(rPeij) + vpColVector(cijPo)) ;

          b =  (vpColVector)cijPo - (vpColVector)rPeij ;           // A.40

          if (k==0)
          {
            A = As ;
            B = b ;
#if 0           
            auto w = rReij.t();
            std::cout << "rPeij0 R\n";
            std::cout << w.getCol(0) << "\n" << w.getCol(1) << "\n" << w.getCol(2) << std::endl;
            std::cout << "rPeij0 is " << rPeij[0] << " " << rPeij[1] << " " << rPeij[2] << std::endl;
            std::cout << "cijPo0 is " << cijPo[0] << " " << cijPo[1] << " " << cijPo[2] << std::endl;
            std::cout << "As0 is \n "; 
            As.csvPrint (std::cout);
            std::cout<< "b0 " << b[0] << " " << b[1] << " " << b[2] << std::endl;
#endif
          }
          else
          {
            A = vpMatrix::stack(A,As) ;
            B = vpColVector::stack(B,b) ;
          }
          k++ ;
        }
      }
    }
	
    // the linear system is defined
    // x = AtA^-1AtB is solved
    vpMatrix AtA = A.AtA() ;

    vpMatrix Ap ;
    AtA.pseudoInverse(Ap, 1e-6) ; // rank 3
    x = Ap*A.t()*B ;

//     {
//       // Residual
//       vpColVector residual;
//       residual = A*x-B;
//       std::cout << "Residual: " << std::endl << residual << std::endl;

//       double res = 0;
//       for (int i=0; i < residual.getRows(); i++)
// 	res += residual[i]*residual[i];
//       res = sqrt(res/residual.getRows());
//       printf("Mean residual = %lf\n",res);
//     }

    // extraction of theta and U
    double theta ;
    double   d=x.sumSquare() ;
    for (unsigned int i=0 ; i < 3 ; i++) x[i] = 2*x[i]/sqrt(1+d) ;
    theta = sqrt(x.sumSquare())/2 ;
    theta = 2*asin(theta) ;
    //if (theta !=0)
    if (std::fabs(theta) > std::numeric_limits<double>::epsilon())
    {
      for (unsigned int i=0 ; i < 3 ; i++) x[i] *= theta/(2*sin(theta/2)) ;
    }
    else
      x = 0 ;
  }

  // Building of the rotation matrix eRc
  vpThetaUVector xP(x[0],x[1],x[2]);
  vpRotationMatrix eRc(xP);

  {
    vpMatrix A ;
    vpColVector B ;
    // Building of the system for the translation estimation
    // for all couples ij
    vpRotationMatrix I3 ;
    I3.eye() ;
    int k = 0 ;
    for (unsigned int i=0 ; i < nbPose ; i++)
    {
      vpRotationMatrix rRei, ciRo ;
      vpTranslationVector rTei, ciTo ;
      rMe[i].extract(rRei) ;
      cMo[i].extract(ciRo) ;
      rMe[i].extract(rTei) ;
      cMo[i].extract(ciTo) ;


      for (unsigned int j=i+1 ; j < nbPose ; j++)
      {
        {

          vpRotationMatrix rRej, cjRo ;
          rMe[j].extract(rRej) ;
          cMo[j].extract(cjRo) ;

          vpTranslationVector rTej, cjTo ;
          rMe[j].extract(rTej) ;
          cMo[j].extract(cjTo) ;

          vpRotationMatrix rReij = rRej.t() * rRei ;

          vpTranslationVector rTeij = rTej+ (-rTei);

          rTeij = rRej.t()*rTeij ;

          vpMatrix a = vpMatrix(rReij) - vpMatrix(I3);

          vpTranslationVector b ;
          b = eRc*cjTo - rReij*eRc*ciTo + rTeij ;

          if (k==0)
          {
            A = a ;
            B = b ;
          }
          else
          {
            A = vpMatrix::stack(A,a) ;
            B = vpColVector::stack(B,b) ;
          }
          k++ ;

        }
      }
    }

    // the linear system is solved
    // x = AtA^-1AtB is solved
    vpMatrix AtA = A.AtA() ;
    vpMatrix Ap ;
    vpColVector AeTc ;
    AtA.pseudoInverse(Ap, 1e-6) ;
    AeTc = Ap*A.t()*B ;

//     {
//       // residual
//       vpColVector residual;
//       residual = A*AeTc-B;
//       std::cout << "Residual: " << std::endl << residual << std::endl;
//       double res = 0;
//       for (int i=0; i < residual.getRows(); i++)
// 	res += residual[i]*residual[i];
//       res = sqrt(res/residual.getRows());
//       printf("mean residual = %lf\n",res);
//     }

    vpTranslationVector eTc(AeTc[0],AeTc[1],AeTc[2]);

    eMc.insert(eTc) ;
    eMc.insert(eRc) ;
  }
}


#endif
void out(std::ostream & onf, const Eigen::Matrix4f  & q)
{
	for(int i = 0; i < 16; i++)
		onf << q.data()[i] << ' ';
}

Eigen::Matrix4f makerandom(const LTD &ld)
{
  return Eigen::Matrix4f::Identity(); //ld.sample().asMatrix();
}

Eigen::Matrix<float,6,1> distance(Eigen::Matrix4f a, Eigen::Matrix4f b)
{
  return Eigen::Matrix<float,6,1> ::Zero();
/*	SE3group<float> ag(a);
	SE3group<float> bg(b);
	return ag.distance(bg).get();*/
}

// pos and ang
Eigen::Matrix<float,2,1> distance2(Eigen::Matrix4f a, Eigen::Matrix4f b)
{
    //auto q = SE3group<float>(a).distance(SE3group<float>(b)).get();
    return Eigen::Matrix<float,2,1>::Zero(); //{q.segment<3>(0).norm(),q.segment<3>(3).norm()};
}

// pos and ang
template <class T>
T distanceT(Eigen::Matrix<T,4,4> a, Eigen::Matrix<T,4,4> b)
{
    return (a*b.inverse()).template block<3,1>(0,3).norm();
}

template <class T>
Vector3<T> normalize(Vector3<T>  q)
{
	return q / q.norm();
}

int main(int argc, char const *argv[])
{
#if 0
	Eigen::AngleAxisf aa(0.001*M_PI,normalize(Eigen::Vector3f(0.2,0.3,0.4)));
	Eigen::Quaternionf q(aa);
	Eigen::Vector3f p = quat2paratsai(q);
	Eigen::Quaternionf Q = paratsai2quat(p);
	Eigen::Matrix3f R = paratsai2rot(p);
	Eigen::Quaternionf QR(R);
	Eigen::Vector3f pp = paratsai2paratsaiprime(p);
	Eigen::Vector3f ppp = paratsaiprime2paratsai(pp);

	float pangle = paratsai2theta(p);
	std::cout << "angle original  " << aa.angle() << std::endl;
	std::cout << "angle from para " << pangle << std::endl;
	std::cout << "angle from quat " << Eigen::AngleAxisf(Q).angle() << std::endl;
	std::cout << "angle from rot  " << Eigen::AngleAxisf(R).angle() << std::endl;
	std::cout << q << std::endl;
	std::cout << Q << std::endl;
	std::cout << QR << std::endl;
	std::cout << ppp.transpose() << std::endl;
	std::cout << p.transpose() << std::endl;
#endif	

    if(argc < 2)
        return -1;
    using FT = double;
    std::ifstream inf(argv[1],std::ios::binary);
    std::vector<Eigen::Matrix<FT,4,4> > cMm,rMe;
    Eigen::Matrix<FT,4,4> m1fp,m2fp,m1f,m2f;
    std::cout << "assuming file of pairs of matrices: cMm and cMe" << std::endl;
    while(inf)
    {
        Eigen::Matrix<double,4,4,Eigen::RowMajor> m1,m2;
        inf.read((char*)m1.data(),16*sizeof(double));
        if(!inf)
            break;
        inf.read((char*)m2.data(),16*sizeof(double));
        if(!inf)
            break;
        // data is stored rowmajor
        m1f = m1.cast<FT>();
        m2f = m2.cast<FT>();
        cMm.push_back(m1f);
        rMe.push_back(m2f);
        if(cMm.size() == 1)
        {
            std::cout << "m1[0] is\n" << m1f << std::endl;
            std::cout << "m2[0] is\n" << m2f << std::endl;
        }
        else
        {
            // previous is valid

            // use 2D (pos,rot) or 1D (pos) distance
            //std::cout << "Diff m1: " << distance2(m1f,m1fp).transpose() << " Diff m2: " << distance2(m2f,m2fp).transpose() << std::endl;
            std::cout << "Diff m1: " << distanceT(m1f,m1fp) << " Diff m2: " << distanceT(m2f,m2fp) << std::endl;
        }
        m1fp = m1f;
        m2fp = m2f;

    }
    Eigen::Matrix<FT,4,4> r;

        FT res = calibrationTsai(cMm,rMe,r);


        std::cout << "Output:   \n" << r << std::endl;
        std::cout << "Residual:   \n" << res << std::endl;
        /*std::cout << "Expected: \n" << ec << std::endl;
	std::cout << "THIS diff:  " << distance(r,ec).transpose() << std::endl;
	std::cout << "!THIS error: " << distance(r,ec).transpose().norm() << std::endl;
    */
#ifdef USEVISP
	vpCalibration c;
	std::vector<vpHomogeneousMatrix> vcMo(cMm.size());
	std::vector<vpHomogeneousMatrix> vrMe(cMm.size());
	vpHomogeneousMatrix veMc;
    Eigen::Matrix4f eMc;
	for(int i = 0; i < cMm.size(); i++)
	{
		assign(vcMo[i],cMm[i].cast<float>());
		assign(vrMe[i],rMe[i].cast<float>());
	}
         calibrationTsai(vcMo,vrMe,veMc);
	assign(eMc,veMc);

        std::cout << "MY output res:\n" << res << std::endl;
        std::cout << "VISP output:\n" << eMc << std::endl;
  std::ofstream onf("tsaiout.bin",std::ios::binary);
  Matrix4<double> out(eMc.cast<FT>());
  onf.write((char*)out.data(),sizeof(double)*16);
/*	std::cout << "Expected:   \n" << ec << std::endl;
	std::cout << "VISP diff:  " << distance(eMc,ec).transpose() << std::endl;
	std::cout << "!VISP error:" << distance(eMc,ec).transpose().norm() << std::endl;
*/
#endif
#if 0
	if(argc == 2)
	{
		std::ofstream onf(argv[1]);
		out(onf,ec);		onf << std::endl;
		out(onf,r);		onf << std::endl;
		onf << cMm.size() << std::endl;
		for(int i = 0; i < cMm.size(); i++)
		{
			out(onf,cMm[i]);
			onf << std::endl;
		}
		for(int i = 0; i < rMe.size(); i++)
		{
			out(onf,rMe[i]);
			onf << std::endl;
		}
	}
#endif
	return 0;
}

/*
 * 0.000  1.000  0.000  0.038
-1.000 -0.000  0.000  0.012
0.000  0.000  1.000  0.015
0.000  0.000  0.000  1.000
*/
