// Author    : Aniruddha Shembekar, Research Engineer, University of Southern California

#include "mex.h"
#include "matrix.h"
#include <iostream>
#include "string.h"
#include </usr/local/include/eigen3/Eigen/Eigen>
#include <stdio.h>
#include <vector>
#include <cmath>

Eigen::MatrixXd eul2bxbybz(Eigen::MatrixXd);
Eigen::Matrix3d eul2rot(Eigen::MatrixXd);	 // fixed ZYX seq

void mexFunction (int _OutArgs, mxArray *MatlabOut[], int _InArgs, const mxArray *MatlabIn[] )
{
    // Define Input
    int cba_cols = 3;
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> cba (mxGetPr(MatlabIn[0]), mxGetNumberOfElements(MatlabIn[0])/cba_cols, cba_cols); 
     
    // Method 
    Eigen::MatrixXd bxbybz = eul2bxbybz(cba);
    Eigen::MatrixXd bx(bxbybz.rows(),3);
    bx << bxbybz.block(0,0,bxbybz.rows(),3);
    Eigen::MatrixXd by(bxbybz.rows(),3);
    by << bxbybz.block(0,3,bxbybz.rows(),3);
    Eigen::MatrixXd bz(bxbybz.rows(),3);
    bz << bxbybz.block(0,6,bxbybz.rows(),3);

    // Define Output
    MatlabOut[0] = mxCreateDoubleMatrix(bx.rows(), 3, mxREAL);
    MatlabOut[1] = mxCreateDoubleMatrix(by.rows(), 3, mxREAL);
    MatlabOut[2] = mxCreateDoubleMatrix(bz.rows(), 3, mxREAL);
    
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M0 ( mxGetPr(MatlabOut[0]), bx.rows(), 3 );
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M1 ( mxGetPr(MatlabOut[1]), by.rows(), 3 );
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M2 ( mxGetPr(MatlabOut[2]), bz.rows(), 3 );
    
    M0 = bx.array();  
	M1 = by.array();  
	M2 = bz.array();  
}

Eigen::MatrixXd eul2bxbybz(Eigen::MatrixXd eul_angles)
{
	// input euler alpha,beta,gamma for ZYX (in radians)
	Eigen::MatrixXd bxbybz = Eigen::MatrixXd::Constant(eul_angles.rows(),9,0);
	for (unsigned int i=0;i<eul_angles.rows();++i)
	{
		Eigen::Matrix3d R = eul2rot(eul_angles.row(i));
		bxbybz.row(i) << R(0,0),R(1,0),R(2,0),R(0,1),R(1,1),R(2,1),R(0,2),R(1,2),R(2,2);
	}
	return bxbybz;
}

Eigen::Matrix3d eul2rot(Eigen::MatrixXd eul_angles)
{
	std::string seq = "ZYX";
	Eigen::Matrix3d rot_mat = Eigen::Matrix3d::Identity();
	for (int i=0; i<3; ++i)
	{
		if(seq[i]=='X')
			rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitX());
		else if(seq[i]=='Y')
			rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitY());			
		else if(seq[i]=='Z')
			rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitZ());					
	}
	return rot_mat; 
}