function r = lines_intersect(l1,l2)

r = -1;
% l1 for horizontal ray line...slope is always zero

if l2(1,2)==l2(2,2) % checking if other slope is zero
    r = 0;
else
    % checking both pts of second line above first line
    if (l2(1,2)>l1(1,2) && l2(2,2)>l1(1,2)) || (l2(1,2)<l1(1,2) && l2(2,2)<l1(1,2))
        r = 0;
    else
        % checking both pts of second line either on right or on left of fist line
        if (l2(1,1)<l1(1,1) && l2(2,1)<l1(1,1)) || (l2(1,1)>l1(2,1) && l2(2,1)>l1(2,1))
            r = 0;
        else
            % checking if other line is vertical
            if l2(1,1)== l2(2,1)
               r=1; 
            else
                % getting intersection point
                m2 = (l2(2,2)-l2(1,2))/(l2(2,1)-l2(1,1));
                x = round(((l1(1,2)+m2*l2(1,1)-l2(1,2))/m2),6);
                % checking if intersection point lies on the first line
                if x>=l1(1,1) && x<=l1(2,1)
                    r = 1;
                else
                    r=0;
                end
            end
        end
    end
end

end