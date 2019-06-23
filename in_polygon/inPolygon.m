function in = inPolygon(q,p)

%     p : polygon points
%     q : querry points
    

if p(1,1)==p(end,1) && p(1,2)==p(end,2)
    p = p(1:end-1,:);
end

in = zeros(size(q,1),1);
xmin = min(p(:,1));
xmax = max(p(:,1));
ymin = min(p(:,2));
ymax = max(p(:,2));

for i=1:size(q,1)
    % bounding box test
    if q(i,1)<xmin || q(i,1)>xmax || q(i,2)<ymin || q(i,2)>ymax
        continue;
    end
    intersection_count = 0;
    cont_lines = zeros(size(p,1),1);
    for j=1:size(p,1)
       if j==1
           if lines_intersect([q(i,:);xmax,q(i,2)],[p(end,:);p(j,:)])==1
               intersection_count = intersection_count+1;
               cont_lines(j,1) = 1;
           end
       else
           if lines_intersect([q(i,:);xmax,q(i,2)],[p(j,:);p(j-1,:)])==1
               intersection_count = intersection_count+1;
               cont_lines(j,1) = 1;
               if cont_lines(j-1,1)==1
                  if p(j-1,2)== q(i,2)
                      if j-1==1
                          if (p(end,2)<p(j-1,2) && p(j,2)<p(j-1,2)) || (p(end,2)>p(j-1,2) && p(j,2)>p(j-1,2))
                          else
                              intersection_count = intersection_count -1;
                          end
                      else
                          if (p(j-2,2)<p(j-1,2) && p(j,2)<p(j-1,2)) || (p(j-2,2)>p(j-1,2) && p(j,2)>p(j-1,2))
                          else
                              intersection_count = intersection_count -1;
                          end
                      end
                          
                  end
               end
           end
       end
    end
    if rem(intersection_count,2)==1
        in(i,1) = 1;
    end
end
end