%****************************************************************************************
%
% Author : Aniruddha Shembekar, University of Southern California
%
%****************************************************************************************

clc;
clear;
close all;

use_fixed_polygon = true;   % false for random 4 sided polygon
use_int_querries = false;    % false for float querries

if use_fixed_polygon 
    % fixed polygon
    p = [    -1    -8;
         9     8;
        -7     2;
        -5     1;
        -7    -7];
else
    %variable polygon
    p = -10 + rand(4,2)*20;
end

if use_int_querries
    % generate integer querries
    q = randi([-15,15],10000,2);
else
    % generate float querries
    q = -15 + rand(10000,2)*30;
end

%% Custom MATLAB implementation
tic;
in = inPolygon(q,p);
toc;

%% MEX implementation
tic;
in2 = inpolygon_mex(q,p);
toc;

%% MATLAB Default function
tic;
in3 = inpolygon(q(:,1),q(:,2),p(:,1),p(:,2));
toc;

%% plotting MEX function output data
figure;
xlim([-15,15]);ylim([-15,15]);
hold on;
p = [p;p(1,:)];
plot(p(:,1),p(:,2));
hold on;
qr = [];qg = [];
for i=1:size(in2,1)
    if in2(i,1)==0
        qr = [qr;q(i,:)];
    else
       qg = [qg;q(i,:)];
    end
end
scatter(qr(:,1),qr(:,2),'.r');
scatter(qg(:,1),qg(:,2),'.g');
