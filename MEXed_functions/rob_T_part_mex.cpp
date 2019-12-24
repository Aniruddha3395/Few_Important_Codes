// Author    : Aniruddha Shembekar, Research Engineer, University of Southern California

#include "mex.h"
#include "matrix.h"
#include <iostream>
#include "string.h"
#include </usr/local/include/eigen3/Eigen/Eigen>
#include <stdio.h>
#include <vector>
#include <cmath>

Eigen::Matrix4d get_rob_T_part(Eigen::MatrixXd, Eigen::MatrixXd);
Eigen::MatrixXd mean(Eigen::MatrixXd);

void mexFunction (int _OutArgs, mxArray *MatlabOut[], int _InArgs, const mxArray *MatlabIn[] )
{
    // Define Input
    int data_cols = 3;
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> part_pts (mxGetPr(MatlabIn[0]), mxGetNumberOfElements(MatlabIn[0])/data_cols, data_cols); 
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> rob_pts (mxGetPr(MatlabIn[1]), mxGetNumberOfElements(MatlabIn[1])/data_cols, data_cols); 
     
    // Method 
    Eigen::Matrix4d transform_mat = get_rob_T_part(part_pts, rob_pts);
    
    // Define Output
    MatlabOut[0] = mxCreateDoubleMatrix(4, 4, mxREAL);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M0 ( mxGetPr(MatlabOut[0]),4,4);
    M0 = transform_mat.array();  
}

Eigen::Matrix4d get_rob_T_part(Eigen::MatrixXd part_pts, Eigen::MatrixXd rob_pts)
{
// part_pts: co-ordinates with respect to CAD part frame
// rob_pts: points with respect to robot base frame (or world frame if it matches with robot base frame)
// input: part_pts = [x1, y1, z1;
//                    x2, y2, z2;
//                         :
//                         :
//                    xn, yn, zn]
// input: rob_pts = [x1, y1, z1;
//                   x2, y2, z2;
//                         :
//                         :
//                   xn, yn, zn]
    Eigen::MatrixXd centroid_part_pts(1,part_pts.cols());
    Eigen::MatrixXd centroid_rob_pts(1,rob_pts.cols());
    Eigen::MatrixXd shifted_part_pts(part_pts.rows(),part_pts.cols());
    Eigen::MatrixXd shifted_rob_pts(rob_pts.rows(),rob_pts.cols());
    Eigen::MatrixXd cros_cov_mat(part_pts.cols(),rob_pts.cols());
    Eigen::Matrix3d R;
    Eigen::Matrix3d U_T;
    Eigen::Matrix3d V;
    Eigen::Vector3d T;
    Eigen::Matrix3d M = Eigen::Matrix3d::Identity();
    Eigen::Matrix4d transform_mat = Eigen::Matrix4d::Constant(0);
    if (part_pts.rows()==rob_pts.rows())
    { 
        centroid_part_pts = mean(part_pts);
        centroid_rob_pts = mean(rob_pts);
        shifted_part_pts.col(0) = part_pts.col(0).array() - centroid_part_pts(0,0);
        shifted_part_pts.col(1) = part_pts.col(1).array() - centroid_part_pts(0,1);
        shifted_part_pts.col(2) = part_pts.col(2).array() - centroid_part_pts(0,2);
        shifted_rob_pts.col(0) = rob_pts.col(0).array() - centroid_rob_pts(0,0);
        shifted_rob_pts.col(1) = rob_pts.col(1).array() - centroid_rob_pts(0,1);
        shifted_rob_pts.col(2) = rob_pts.col(2).array() - centroid_rob_pts(0,2);
        cros_cov_mat = shifted_part_pts.transpose()*shifted_rob_pts;
    
        // Singular Value Decomposition
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(cros_cov_mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
        
        // Take care of reflection case due to negative eigen vectors
        U_T = svd.matrixU().transpose();    V = svd.matrixV();
        M(2,2) = (V*U_T).determinant();
        R = V*M*U_T;
        if (R.determinant()>0)
        {
            T = -R*centroid_part_pts.transpose() + centroid_rob_pts.transpose();
            transform_mat.block(0,0,3,3) = R;
            transform_mat.block(0,3,3,1) = T;
            transform_mat(3,3) = 1; 
        }
        else
        {
            std::cerr << "ERROR: Determinant of rotation matrix is negative..." << std::endl;
        }   
    }
    else
    {
        std::cerr << "ERROR: FUNCTION ERROR: For correspondance, number of rows of both matrices should be same..." << std::endl;
    }
    return transform_mat;
}

Eigen::MatrixXd mean(Eigen::MatrixXd mat)
{
    Eigen::VectorXd vec(mat.cols());
    for (int i=0;i<mat.cols();++i)
    {
        vec(i) =  mat.block(0,i,mat.rows(),1).sum()/mat.rows();
    }
    return vec.transpose();
}