function [q] = mysignal(x)
%%%%%%%%%%%%%%
%% Input
%  x: some x values
%% Output
%  q: this will be the function, in this case a mollifier
%%%%%%%%%%%%%%
Index = x.^2<0.999999;
xn = x.*Index;
q =  1*exp(-1./(1-xn.^2)).*Index;
end