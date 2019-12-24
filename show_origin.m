% Author    : Aniruddha Shembekar, Research Engineer, University of
% Southern California

function show_origin()

% code to show origin in the figure

% use this for meters
dir_vec_len = 0.050;
% use this for mm
% dir_vec_len = 50;

hold on;
scatter3(0,0,0,50,'d','filled','k');
hold on;
quiver3(0,0,0,1,0,0,dir_vec_len,'r');
hold on;
quiver3(0,0,0,0,1,0,dir_vec_len,'g');
hold on;
quiver3(0,0,0,0,0,1,dir_vec_len,'b');
hold on;
end