clc;
clear;

n = input('number of masses in the system = ');

for i = 1:n
    m{i,1} = input('give name to mass = ','s');
    x{i,1} = input('give name to its displacement = ','s');
    fprintf('in which direction does this mass go?\n');
    d{i,1} = input('1.right(enter r)\n2.left(enter l)\n=','s');
end

for i = 1:n
    fprintf('is there any spring attached to left side of %s ?\n',m{i,1})
    L{i,1} = input('1.yes(enter y)\n2.no(enter n) \n=','s');
    if L{i,1}=='y'
        Lk{i,1} = input('enter spring constant = \n','s');
        fprintf('is there any mass attached to left side of this spring?\n');
        LLk{i,1} = input('1.yes(enter y)\n2.no(enter n) \n=','s');
        if LLk{i,1} == 'y'
            LLx{i,1} = input('enter its displacement = ','s');
            fprintf('enter its direction = \n');
            LLd{i,1} = input('1.right(enter r)\n2.left(enter l)\n=','s');
        end
    end 
end

for i = 1:n
    fprintf('is there any spring attached to right side of %s ?\n',m{i,1})
    R{i,1} = input('1.yes(enter y)\n2.no(enter n) \n=','s');
    if R{i,1}=='y'
        Rk{i,1} = input('enter spring constant = \n','s');
        fprintf('is there any mass attached to right side of this spring?\n');
        RRk{i,1} = input('1.yes(enter y)\n2.no(enter n) \n=','s');
        if RRk{i,1} == 'y'
            RRx{i,1} = input('enter its displacement = ','s');
            fprintf('enter its direction = \n');
            RRd{i,1} = input('1.right(enter r)\n2.left(enter l)\n=','s');
        end
    end 
end

for i = 1:n
    if d{i,1}=='r'
        Lkx(i,1) = {'+'};
        Rkx(i,1) = {'-'};
    else
        Lkx(i,1) = {'-'};
        Rkx(i,1) = {'+'};
    end
    if LLd{i,1}=='r'
        LLxLk(i,1) = {'-'};
    else
        LLxLk(i,1) = {'+'};
    end
    if RRd{i,1} == 'r'
        RRxRk(i,1) = {'+'};
    else
        RRxRk(i,1) = {'-'};
    end
   fprintf('equation on LHS of mass %s pointing leftwards  = \n',m{i,1});
   fprintf('(%s %s*%s %s %s*%s)\n',Lkx{i,1},Lk{i,1},x{i,1},LLxLk{i,1},LLx{i,1},Lk{i,1});
   fprintf('equation on RHS of mass %s pointing leftwards  = \n',m{i,1});
   fprintf('(%s %s*%s %s %s*%s)\n',Rkx{i,1},Rk{i,1},x{i,1},RRxRk{i,1},RRx{i,1},Rk{i,1});
end