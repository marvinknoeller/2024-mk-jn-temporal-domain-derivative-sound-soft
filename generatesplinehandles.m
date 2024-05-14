function perspline = generatesplinehandles(vals)
%%%%%%%%%%%%%%
%% Input
%  vals: the values of the points to be interpolated by a periodic cubic
%  spline 
%% Output
%  perspline: the spline
%%%%%%%%%%%%%%
assert(norm(vals(:,1) - vals(:,end))<1e-14 )
x = linspace(0,2*pi,length(vals));
perspline = csape(x,vals,'periodic');
assert(norm(perspline.breaks(end)-2*pi)<1e-14)
end