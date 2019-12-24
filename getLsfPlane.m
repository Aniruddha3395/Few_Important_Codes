% Author    : Aniruddha Shembekar, Research Engineer, University of
% Southern California

function [A,B,C,D] = getLsfPlane(pts_for_plane)

x_avg = sum(pts_for_plane(:,1))/size(pts_for_plane,1);
y_avg = sum(pts_for_plane(:,2))/size(pts_for_plane,1);
z_avg = sum(pts_for_plane(:,3))/size(pts_for_plane,1);

L00 = sum((pts_for_plane(:,1)-x_avg).^2);
L01 = sum((pts_for_plane(:,1)-x_avg).*(pts_for_plane(:,2)-y_avg));
L11 = sum((pts_for_plane(:,2)-y_avg).^2);
R0 = sum((pts_for_plane(:,3)-z_avg).*(pts_for_plane(:,1)-x_avg));
R1 = sum((pts_for_plane(:,3)-z_avg).*(pts_for_plane(:,2)-y_avg));

A = -((L11*R0-L01*R1)/(L00*L11-L01^2));
B = -((L00*R1-L01*R0)/(L00*L11-L01^2));
C = 1;
D = -(z_avg+A*x_avg+B*y_avg);

end