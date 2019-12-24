%****************************************************************************************
%
% Author : Aniruddha Shembekar, University of Southern California
%
%****************************************************************************************

function showTransformFrame(T)

% code to show origin in the figure

% use this for meters
dir_vec_len = 0.050;
% use this for mm
% dir_vec_len = 50;

hold on;
scatter3(T(1,4),T(2,4),T(3,4),50,'d','filled','k');
hold on;
quiver3(T(1,4),T(2,4),T(3,4),T(1,1),T(2,1),T(3,1),dir_vec_len,'r');
hold on;
quiver3(T(1,4),T(2,4),T(3,4),T(1,2),T(2,2),T(3,2),dir_vec_len,'g');
hold on;
quiver3(T(1,4),T(2,4),T(3,4),T(1,3),T(2,3),T(3,3),dir_vec_len,'b');
hold on;

end