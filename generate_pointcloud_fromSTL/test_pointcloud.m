% Author    : Aniruddha Shembekar, Research Engineer, University of
% Southern California

clc;
clear;
close all;

% test pointcloud generation

clc;
clear;
close all;
set(0, 'DefaultFigureRenderer', 'opengl');

%% adding woking directory and all dependancies
[working_dir,~,~] = fileparts(mfilename('fullpath'));
addpath(genpath(working_dir),'-end');

%% import the STL part and its corresponding pointcloud

part_name = 'test9';

stl_file_name = strcat('CAD_stl/',part_name,'.STL');
ptcloud_file_name = strcat('data_files/',part_name,'_ptcloud_mm_gap1mm.csv');

[v,f,n, ~] = stlRead(stl_file_name);

gap_x = 0.5;
gap_y = 0.5;

% tic;
% [pts_set,normals] = generate_pointcloud(v,f,n,gap_x,gap_y);
% toc;
tic;
[pts_set,normals] = generate_pointcloud_mex(v,f,n,gap_x,gap_y);
toc;

figure;
show_origin();
hold on;
view_stl_with_VF(v,f,[],[],[],[]);
hold on;
scatter3d(pts_set,'.');
hold on;
% quiver3(pts_set(:,1),pts_set(:,2),pts_set(:,3),normals(:,1),normals(:,2),normals(:,3),5,'r');