clc;
clear;
close all;

%% addddd graph

%give pairs of connected nodes
% x = [1 1 1 2 3 2 2 3 3 4 5 5 7 6];
% y = [2 3 4 3 4 5 7 5 6 6 6 8 8 8];
% wt = [2 5 4 2 1 7 12 4 3 4 1 5 3 7];
% names = {'O' 'A' 'B' 'C' 'D' 'E' 'F' 'T'};


% fileID = fopen('file.txt','r');
% format = '%d %d %d';
% a = fscanf(fileID,format,[3 Inf])';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test case generator
total_nodes = 5;
edges = 10;
wll = 1;
wul = 100;
a = zeros(edges,3);
for i=1:edges
    flag = 0;
    while flag==0
        a(i,1:2) = randi([1,total_nodes],1,2);
        if a(i,1)~=a(i,2)
            flag = 1;
        end
    end
end
a = unique(a,'rows');
a(:,3) = randi([wll,wul],size(a,1),1);
a = [a,(a(:,1)).^2+(a(:,2)).^2];
[~, rows] = unique(a(:, 4));
b = a(rows, :);
a = b(:,1:3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x = a(:,1)';
y = a(:,2)';
wt = a(:,3)';
names = [1:total_nodes];
link = [x',y'];
G = graph(x,y,wt);
figure;
graph_plot = plot(G,'EdgeLabel',G.Edges.Weight);

%% Dijkstra

%input
start_node = input('input start node: ');       %start node
end_node = input('input end node: ');       %end node

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
fprintf('----------------------------------------------------------\n');
[P,d] = shortestpath(G,start_node,end_node);       %Matlab Answer for reference
MATLAB_Answer = d
toc;
fprintf('----------------------------------------------------------\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
tic;
inf = 10000;        % large number given as a weight to all other nodes at the start
save_start_node = start_node;
save_end_node = end_node;
visit = [];         % visited nodes
N_visit = [1:size(names,2)];        % unvisited nodes
node = [1:size(names,2)]';          % nodes range
dist = [inf*ones(1,size(names,2))]';        % distnace from the start node
parent = [zeros(1,size(names,2))]';         % parent node
tab = [node,dist,parent];            % data containing node, weight and parent
tab(find(tab(:,1)==start_node),2)=0;
Initial_table = tab;

%algorithm
while any(visit==end_node)==0
% looking at the first column of link
    indx = find(link(:,1)==start_node);
    if isempty(indx)==0
        for i=1:size(indx)
            if any(visit==link(indx(i),2))==0
                new_wt = wt(indx(i))+ tab(find(tab(:,1)==start_node),2);
                wt_loc = find(tab(:,1)==link(indx(i),2));
                if new_wt < tab(wt_loc,2)
                    tab(wt_loc,2) = new_wt;
                    tab(wt_loc,3) = start_node;
                end
            end
        end
    end
    
% looking at the second column of link
    indx = find(link(:,2)==start_node);
    if isempty(indx)==0
        for i=1:size(indx)
            if any(visit==link(indx(i),1))==0
                new_wt = wt(indx(i))+ tab(find(tab(:,1)==start_node),2);
                wt_loc = find(tab(:,1)==link(indx(i),1));
                if new_wt < tab(wt_loc,2)
                    tab(wt_loc,2) = new_wt;
                    tab(wt_loc,3) = start_node;
                end
            end
        end
    end
    
    % looking in the table for the node with lowest weight (but not zero) and which is not
    % visited yet.
    visit = [visit,start_node];
    N_visit(find(N_visit)==start_node) = [];
    store = inf;
    for i=1:size(tab,1)
        if any(visit==tab(i,1))==0 && tab(i,2)~=0 && tab(i,2)<store
            store = tab(i,2);
            store_idx = i;
        end
    end
    start_node = tab(store_idx,1);
end
Node = names';
Weight = tab(:,2);
Parent_Node = tab(:,3);
Final_table = table(Node,Weight, Parent_Node);

% visualization

track_path = [end_node];
counter = 0;
while any(track_path==save_start_node)==0
    counter = counter+1;
    end_node_loc = find(tab(:,1)==end_node);
    track_path = [track_path,tab(end_node_loc,3)];
    end_node = track_path(end);
end

hold on;
highlight(graph_plot,track_path,'EdgeColor','r','LineWidth',1.5)
fprintf('----------------------------------------------------------\n');
Shortest_Distance = tab(find(tab(:,1)==save_end_node),2)
toc;
fprintf('----------------------------------------------------------\n');



