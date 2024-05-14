function [x,y] = F_op(freq,S,V,b,d1,dist)
%%%%%%%%%%%%%%
%% Input
%  freq: this is s
%  S: function handle of the single layer potential
%  V: function handle of the single layer operator, i.e. V=S|_{\partial D}
%  b: right hand side 
%  d1: number of evaluation points in space
%  dist: the minimum distance to the measurement points
%% Output
%  x: this is SV^{-1}b
%  y: this is V^{-1}b (interesting, when requiring Neumann data of u+u^i)
%%%%%%%%%%%%%%
% little abuse of the function here: 
if length(b) == 1
    x = size(V,2);
    return
end

%%
if nargout == 2
    if exp(-real(freq)*dist)*norm(freq^(5/2)*b)<10^(-12)
        x = zeros(d1,size(b,2));
        y = zeros(size(S,2),1);
    else
        y = V\b;
        x = S*y;
    end
elseif nargout == 1
    if exp(-real(freq)*dist)*norm(freq^(5/2)*b)<10^(-12)
        x = zeros(d1,size(b,2));
    else
        x = S*(V\b);
    end
end

end