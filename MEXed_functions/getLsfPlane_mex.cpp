// Author    : Aniruddha Shembekar, Research Engineer, University of Southern California

#include "mex.h"
#include "matrix.h"
#include <iostream>
#include "string.h"
#include </usr/local/include/eigen3/Eigen/Eigen>
#include <stdio.h>
#include <vector>

Eigen::MatrixXd getLsfPlane(const Eigen::MatrixXd&);

void mexFunction (int _OutArgs, mxArray *MatlabOut[], int _InArgs, const mxArray *MatlabIn[] )
{
    // Define Input
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> pts (mxGetPr(MatlabIn[0]), mxGetNumberOfElements(MatlabIn[0])/3, 3); 

    // Method 
    Eigen::MatrixXd T = getLsfPlane(pts);

    // Define Output
    MatlabOut[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    MatlabOut[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    MatlabOut[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    MatlabOut[3] = mxCreateDoubleMatrix(1,1,mxREAL);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M0 (mxGetPr(MatlabOut[0]),1,1);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M1 (mxGetPr(MatlabOut[1]),1,1);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M2 (mxGetPr(MatlabOut[2]),1,1);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M3 (mxGetPr(MatlabOut[3]),1,1);
    M0 = T(0,0);
    M1 = T(0,1);
    M2 = T(0,2);
    M3 = T(0,3); 
}

Eigen::MatrixXd getLsfPlane(const Eigen::MatrixXd& pts_for_plane)
{
    Eigen::MatrixXd x_vec(pts_for_plane.rows(),1);
    Eigen::MatrixXd y_vec(pts_for_plane.rows(),1);
    Eigen::MatrixXd z_vec(pts_for_plane.rows(),1);
    double x_avg, y_avg, z_avg, L00, L11, L01, R0, R1;   
    Eigen::MatrixXd ret(1,4);
    x_vec = pts_for_plane.block(0,0,pts_for_plane.rows(),1);
    y_vec = pts_for_plane.block(0,1,pts_for_plane.rows(),1);
    z_vec = pts_for_plane.block(0,2,pts_for_plane.rows(),1); 
    x_avg = x_vec.sum()/pts_for_plane.rows();
    y_avg = y_vec.sum()/pts_for_plane.rows();
    z_avg = z_vec.sum()/pts_for_plane.rows();
    L00 = ((x_vec.array() - x_avg).array().pow(2)).sum();
    L01 = ((x_vec.array() - x_avg).array()*(y_vec.array() - y_avg).array()).sum(); 
    L11 = ((y_vec.array() - y_avg).array().pow(2)).sum();
    R0 = ((z_vec.array() - z_avg).array()*(x_vec.array() - x_avg).array()).sum(); 
    R1 = ((z_vec.array() - z_avg).array()*(y_vec.array() - y_avg).array()).sum(); 
    ret(0,0) = -((L11*R0-L01*R1)/(L00*L11-L01*L01));
    ret(0,1) = -((L00*R1-L01*R0)/(L00*L11-L01*L01));
    ret(0,2) = 1;
    ret(0,3) = -(z_avg+ret(0,0)*x_avg+ret(0,1)*y_avg);
    return ret;
}