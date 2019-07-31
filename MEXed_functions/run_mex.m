% run this file to mex all the .cpp codes.
% NOTE: need to have the eigen c++ libarary installed (http://eigen.tuxfamily.org/index.php?title=Main_Page)

mex -R2018a apply_transformation_mex.cpp;
mex -R2018a bxbybz_to_euler_mex.cpp;
mex -R2018a euler_to_bxbybz_mex.cpp;
mex -R2018a rob_T_part_mex.cpp;
mex -R2018a Rotx_mex.cpp;
mex -R2018a Roty_mex.cpp;
mex -R2018a Rotz_mex.cpp;
mex -R2018a eul2rotm_mex.cpp;  
mex -R2018a rotm2eul_mex.cpp;  
mex -R2018a quat2rotm_mex.cpp;
mex -R2018a rotm2quat_mex.cpp;
mex -R2018a quat2eul_mex.cpp;
mex -R2018a eul2quat_mex.cpp;
mex -R2018a getLsfPlane_mex.cpp;

