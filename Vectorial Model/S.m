
% Function to form skew symmetric matrix

function [s] = S(x1,x2,x3)
        s = [0 -x3 x2;
             x3 0 -x1;
             -x2 x1 0];
end