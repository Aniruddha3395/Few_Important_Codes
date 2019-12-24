% Author    : Aniruddha Shembekar, Research Engineer, University of
% Southern California

function [xyz_n_transformed] = apply_transformation_to_xyz_normals(xyz_n,transformation_matrix)
% code for applying transformation on xyz points as well as normals
% vectors bx, by and bz

xyz_n_transformed = zeros(size(xyz_n));
xyz_n_transformed(:,1:3) = apply_transformation(xyz_n(:,1:3),transformation_matrix);

% making translation component of transformation matrix to zero
transformation_matrix(1:3,4) = 0;
xyz_n_transformed(:,4:6) = apply_transformation(xyz_n(:,4:6),transformation_matrix);

end