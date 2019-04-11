#include "mex.h"
#include "matrix.h"
#include <iostream>
#include "string.h"
#include </usr/local/include/eigen3/Eigen/Eigen>
#include <stdio.h>
#include <vector>

Eigen::MatrixXd eul2rot(Eigen::MatrixXd, std::string);
std::string validate_seq(std::string);

void mexFunction (int _OutArgs, mxArray *MatlabOut[], int _InArgs, const mxArray *MatlabIn[] )
{
    if (_InArgs==2)
    {
        // Define Input
        std::string seq_in = "";
        Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> eul_angles (mxGetPr(MatlabIn[0]), 1, 3); 
        char* pCharArray = mxArrayToString(MatlabIn[1]); 
        std::string str(pCharArray);
        seq_in = str;    
            
        // Method 
        Eigen::MatrixXd m = eul2rot(eul_angles, seq_in);
        
        // Define Output
        MatlabOut[0] = mxCreateDoubleMatrix(3,3,mxREAL);
        Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M0 (mxGetPr(MatlabOut[0]),3,3);
        M0 = m.array(); 
    }
    else
    {
        // Define Input
        Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> eul_angles (mxGetPr(MatlabIn[0]), 1, 3); 
            
        // Method 
        std::string seq_in = "";
        Eigen::MatrixXd m = eul2rot(eul_angles, seq_in);
        
        // Define Output
        MatlabOut[0] = mxCreateDoubleMatrix(3,3,mxREAL);
        Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M0 (mxGetPr(MatlabOut[0]),3,3);
        M0 = m.array();
    }
}

Eigen::MatrixXd eul2rot(Eigen::MatrixXd eul_angles, std::string seq)
{
    seq = validate_seq(seq);
    Eigen::MatrixXd rot_mat = Eigen::MatrixXd::Identity(3,3);
    for (int i=0; i<3; ++i)
    {
        if(seq[i]=='X' || seq[i]=='x')
            rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitX());
        else if(seq[i]=='Y' || seq[i]=='y')
            rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitY());           
        else if(seq[i]=='Z' || seq[i]=='z')
            rot_mat = rot_mat * Eigen::AngleAxisd(eul_angles(0,i), Eigen::Vector3d::UnitZ());                   
    }
    return rot_mat; 
}

std::string validate_seq(std::string seq)
{
    if(seq =="")
        seq = "ZYX";    
    bool invalid_flag = false;
    if(seq.size()!=3)
    {
        invalid_flag = true;
    }
    for (int i =0;i<3;++i)
    {
        if(seq[i]!='X' && seq[i]!='Y' && seq[i]!='Z' && seq[i]!='x' && seq[i]!='y' && seq[i]!='z')
        {
            invalid_flag = true; 
            break;
        }
    }
    if(invalid_flag)
    {
        std::cerr << "ERROR: Invalid Rotations Sequence: returning the value in 'ZYX' sequence " << seq << std::endl;
        return "ZYX";       
    }
    return seq;
}