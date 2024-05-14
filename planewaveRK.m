function [u,dnu]=planewaveRK(f,fp,d,tlag,x,normal,T,M,option,c)
%%%%%%%%%%%%%%
%% Input
%  f: function handle of 1d smooth function
%  fp: derivative of f
%  d: direction of propagation
%  tlag: a potential lag (shift to not start at time t=0)
%  x: points in space, where you want to have the wave
%  normal: the exterior unit normals to the boundary
%  T: final time 
%  M: number of steps in time
%  option: either 1 or 0 (do you want the derivative?)
%  c: coefficient from the RK method
%% Output
%  u: the plane wave at the stages needed for the CQ at the points x
%  dnu: see u but for the derivative
%%%%%%%%%%%%%%
tau = T/M;
m = length(c);
t = get_time_points(T,tau, m, c);
dir = x*d.';
u = f(dir - (t-tlag)); 
% important: u is nothing but:
% u(:,1) discribes incoming wave on the boundary at the first time point
% u(:,2) "-" at the second time point, and so on...
if option==0
    dnu=[];
    return
end
nor = normal*d.';
dnuh = fp(dir-(t-tlag));
dnu = nor .* dnuh;


return

