%% original code

%%% Function for factorial expt
% % n are number of factors
% function F = factorial_expt(n)
% lp = 1;
% for i = 1:n
%     fprintf('enter number of levels for factor %d = \n',i);
%     L(1,i) = input('');
%     lp = lp*L(1,i);
% end
% F = zeros(lp,n);
% Lp = lp;
% lp =1;
% for i = 1:n
%         fprintf('value for levels for factor %d\n',i);
%     for j = 1:L(1,i)
%         val(1,j) = input('value = ');
%     end
%     v = repmat(val,[lp 1]);
%     v1 = (v(:)')';
%     lp = lp*L(1,i);
%     
%     Dp = Lp/lp;
%     F(:,i) = repmat(v1',[1 Dp])';
%     v = [];
%     val = [];
% end
% F;
% end


%% modified version
% Function for factorial expt
% n are number of factors
function F = factorial_expt(n)
lp = 16*16*16*16; 
L = [16,16,16,16];
F = zeros(lp,n);
Lp = lp;
lp =1;
for i = 1:n
%         fprintf('value for levels for factor %d\n',i);
    val = [1:16];
    v = repmat(val,[lp 1]);
    v1 = (v(:)')';
    lp = lp*L(1,i);
    
    Dp = Lp/lp;
    F(:,i) = repmat(v1',[1 Dp])';
    v = [];
    val = [];
end
F;
end