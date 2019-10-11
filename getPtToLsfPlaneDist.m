function [dist_p_plane] = getPtToLsfPlaneDist(p,pts_for_plane)

[A,B,C,D] = getLsfPlane(pts_for_plane);
dist_p_plane = abs(A*p(1)+B*p(2)+C*p(3)+D)/(sqrt(A^2+B^2+C^2));

end