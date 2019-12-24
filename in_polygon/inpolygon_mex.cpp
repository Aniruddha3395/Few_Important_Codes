//****************************************************************************************
//
// Author : Aniruddha Shembekar, University of Southern California
//
//****************************************************************************************
        
#include "mex.h"
#include "matrix.h"
#include <iostream>
#include "string.h"
#include </usr/local/include/eigen3/Eigen/Eigen>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <chrono>

Eigen::MatrixXd in_poly(Eigen::MatrixXd , Eigen::MatrixXd );
bool lines_intersect(double l1[2][2], double l2[2][2]);

void mexFunction (int _OutArgs, mxArray *MatlabOut[], int _InArgs, const mxArray *MatlabIn[] )
{
    // Define Input
    // q : querry points
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> q (mxGetPr(MatlabIn[0]), mxGetNumberOfElements(MatlabIn[0])/2, 2); 
    // p : polygon points
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> p (mxGetPr(MatlabIn[1]), mxGetNumberOfElements(MatlabIn[1])/2, 2); 
    
    // Method 
    Eigen::MatrixXd in = in_poly(q,p);

    // Define Output
    MatlabOut[0] = mxCreateDoubleMatrix(in.rows(), 1, mxREAL);
    Eigen::Map<Eigen::ArrayXXd,Eigen::Aligned> M0 ( mxGetPr(MatlabOut[0]), in.rows(), 1 );
    M0 = in.array();  
}

Eigen::MatrixXd in_poly(Eigen::MatrixXd q, Eigen::MatrixXd p)
{
	double l1[2][2];
	double l2[2][2];

	Eigen::MatrixXd in = Eigen::MatrixXd::Constant(q.rows(),1,0);

	double xmin = p.col(0).minCoeff();
	double xmax = p.col(0).maxCoeff();
	double ymin = p.col(1).minCoeff();
	double ymax = p.col(1).maxCoeff();

	for (long i=0;i<q.rows();++i)
	{
		// bounding box test
		if (q(i,0)<xmin || q(i,0)>xmax || q(i,1)<ymin || q(i,1)>ymax)
		{
			continue;
		}
		int intersection_count = 0;
		Eigen::MatrixXd cont_lines = Eigen::MatrixXd::Constant(p.rows(),1,0);
		for (int j=0;j<p.rows();++j)
		{
			if (j==0)
			{
				l1[0][0] = q(i,0);l1[0][1] = q(i,1);
				l1[1][0] = xmax;l1[1][1] = q(i,1);
				l2[0][0] = p(p.rows()-1,0);l2[0][1] = p(p.rows()-1,1);
				l2[1][0] = p(j,0);l2[1][1] = p(j,1);
				if (lines_intersect(l1,l2))
				{
					intersection_count++;
					cont_lines(j,0) = 1;
				}	
			}
			else
			{
				l1[0][0] = q(i,0);l1[0][1] = q(i,1);
				l1[1][0] = xmax;l1[1][1] = q(i,1);
				l2[0][0] = p(j,0);l2[0][1] = p(j,1);
				l2[1][0] = p(j-1,0);l2[1][1] = p(j-1,1);
				if (lines_intersect(l1,l2))
				{
					intersection_count++;
					cont_lines(j,0) = 1;
					if (cont_lines(j-1,0)==1)
					{
						if (p(j-1,1)==q(i,1))
						{
							if (j-1==0)
							{
								if (!((p(p.rows()-1,1)<p(j-1,1) && p(j,1)<p(j-1,1)) || (p(p.rows()-1,1)>p(j-1,1) && p(j,1)>p(j-1,1))))
								{
									intersection_count--;
								}
							}
							else
							{
								if (!((p(j-2,1)<p(j-1,1) && p(j,1)<p(j-1,1)) || (p(j-2,1)>p(j-1,1) && p(j,1)>p(j-1,1))))
								{
									intersection_count--;
								}
							}
						}
					}
				}
			}
		}
		if (intersection_count%2==1)
		{
			in(i,0) = 1;
		}
	}
	return in;
}

bool lines_intersect(double l1[2][2], double l2[2][2])
{
	// l1 for horizontal ray line...slope is always zero

	// checking if other slope is zero
	if (l2[0][1]==l2[1][1])
	{
    	return false;
	}
	else
	{
		// checking both pts of second line above first line
		if ((l2[0][1]>l1[0][1] && l2[1][1]>l1[0][1]) || (l2[0][1]<l1[0][1] && l2[1][1]<l1[0][1]))
		{
			return false;
		}
		else
		{
			// checking both pts of second line either on right or on left of fist line
			if ((l2[0][0]<l1[0][0] && l2[1][0]<l1[0][0]) || (l2[0][0]>l1[1][0] && l2[1][0]>l1[1][0]))
			{
				return false;
			}
			else
			{
				// checking if other line is vertical
				if (l2[0][0]== l2[1][0])
				{
					return true;
				}
				else
				{
					// getting intersection point
					double m2 = (l2[1][1]-l2[0][1])/(l2[1][0]-l2[0][0]);		
					double x = (l1[0][1]+m2*l2[0][0]-l2[0][1])/m2;
					// checking if intersection point lies on the first line
					if ((x>l1[0][0] || std::abs(x-l1[0][0])<1e-9) && (x<l1[1][0] || std::abs(x-l1[1][0])<1e-9))
					{
						return true;
					}
					else
					{
						return false;
					}
				}
			}
		}
	} 
	return false;
}