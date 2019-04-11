#include "mex.h"
#include "matrix.h"
#include <iostream>
#include "string.h"
#include </usr/local/include/eigen3/Eigen/Eigen>
#include <stdio.h>
#include <vector>

Eigen::MatrixXd rot2qt(Eigen::Matrix3d);

void mexFunction (int _OutArgs, mxArray *MatlabOut[], int _InArgs, const mxArray *MatlabIn[] )
{
    // Define Input
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> rot_mat (mxGetPr(MatlabIn[0]), 3, 3); 
    
    // Method 
    Eigen::MatrixXd quat = rot2qt(rot_mat);
    
    // Define Output
    MatlabOut[0] = mxCreateDoubleMatrix(1,4,mxREAL);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M0 (mxGetPr(MatlabOut[0]),1,4);
    M0 = quat.array(); 
}

Eigen::MatrixXd rot2qt(Eigen::Matrix3d rot_mat)
{
    Eigen::MatrixXd quat(1,4);
    Eigen::Quaterniond q(rot_mat);
    quat(0,0) = q.x();
    quat(0,1) = q.y();
    quat(0,2) = q.z();
    quat(0,3) = q.w();
    return quat;
}
