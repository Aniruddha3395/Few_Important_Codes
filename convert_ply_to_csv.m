clc;
clear;
close all;

pathdir = '/home/aniruddha/Downloads/';

filename= 'CapturedPLY';

pcdata = pcread(strcat(pathdir,filename,'.ply'));
pts = zeros(pcdata.Count,3);

for i=1:pcdata.Count
    pts(i,:) = pcdata.Location(i,1:3);
end

pts = pts.*1000;
scatter3d(pts,'.');
dlmwrite(strcat(pathdir,filename,'.csv'),pts);