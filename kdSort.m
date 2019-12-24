%****************************************************************************************
%
% Author : Aniruddha Shembekar, University of Southern California
%
%****************************************************************************************

function points = kdSort(points,m,n)
       sorted_data = [];
       for i = 1:m-1
           origin = [min(points(:,1))-0.1,min(points(:,2))-0.1];
           dist = vecnorm(points(:,1:2)-origin,2,2);
           [~,idx] = min(dist);
           o = points(idx,:);
           sorted_data = [sorted_data;o];
           for i=1:n-1
               points = setdiff(points,o,'rows');
               kdtree = KDTreeSearcher(points);
               if i==1
                   Idx = knnsearch(kdtree,o,'K',3);
                    np = points(Idx,:);
                    [~,id] = min(np(:,1));
               elseif i==n-1
                    Idx = knnsearch(kdtree,o,'K',4);
                    np = points(Idx,:);
                    [~,id] = min(np(:,1));   
               else
                   Idx = knnsearch(kdtree,o,'K',6);
                   np = points(Idx,:);
                   [~,id] = mink(np(:,1),2);
                   if norm(np(id(1),:)-o)<norm(np(id(2),:)-o)
                       id = id(1);
                   else
                       id = id(2);
                   end
               end
               sorted_data = [sorted_data;np(id,:)];
               o = np(id,:);
           end
           points = setdiff(points,o,'rows');
       end
       points = sortrows(points,2);
       sorted_data = [sorted_data;points];
       points = sorted_data;
      
   end

