%****************************************************************************************
%
% Author : Aniruddha Shembekar, University of Southern California
%
%****************************************************************************************

function [Transformation] = rob_T_part(part_pts,rob_pts)

% part with respect to robot transformation 

centroid_part_pts = mean(part_pts);
centroid_rob_pts = mean(rob_pts);
part_pts = part_pts - centroid_part_pts;
rob_pts = rob_pts - centroid_rob_pts;

CrossCovariance_Mat = part_pts'*rob_pts;
[U,~,V] = svd(CrossCovariance_Mat);         % singular value decomposition
R = V*[1 0 0;0 1 0;0 0 det(V*U')]*U';       % to take care of reflection case due to negative eigen vectors

if det(R)>0
    T = -R*centroid_part_pts' + centroid_rob_pts';
    Transformation = [[R,T];0 0 0 1];
else
    fprintf('Determinant of rotation matrix is negative...')
end

end




