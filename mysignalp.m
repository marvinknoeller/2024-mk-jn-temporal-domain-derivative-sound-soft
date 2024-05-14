function [qp] = mysignalp(x)
%%%%%%%%%%%%%%
%% Input
%  x: some x values
%% Output
%  qp: this will be the function, in this case the derivative of a mollifier
%%%%%%%%%%%%%%
Index = x.^2<0.999999;
xn = x.*Index;
qp = -2*xn./((1-xn.^2).^2) .* exp(-1./(1-xn.^2)).*Index;
end