// Author    : Aniruddha Shembekar, Research Engineer, University of Southern California

#include "mex.h"
#include "matrix.h"
#include <iostream>
#include "string.h"
#include </usr/local/include/eigen3/Eigen/Eigen>
#include <stdio.h>
#include <vector>

Eigen::Matrix3d qt2rot(Eigen::MatrixXd);

void mexFunction (int _OutArgs, mxArray *MatlabOut[], int _InArgs, const mxArray *MatlabIn[] )
{
    // Define Input
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> quat (mxGetPr(MatlabIn[0]), 1, 4); 
    
    // Method 
    Eigen::MatrixXd rot_mat = qt2rot(quat);
    
    // Define Output
    MatlabOut[0] = mxCreateDoubleMatrix(3,3,mxREAL);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M0 (mxGetPr(MatlabOut[0]),3,3);
    M0 = rot_mat.array(); 
}

Eigen::Matrix3d qt2rot(Eigen::MatrixXd quat)
{
    Eigen::Quaterniond q;
    q.x() = quat(0,0);
    q.y() = quat(0,1);
    q.z() = quat(0,2);
    q.w() = quat(0,3);  
    return q.normalized().toRotationMatrix();
}
